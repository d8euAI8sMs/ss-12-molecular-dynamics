#pragma once

#include <util/common/geom/point.h>
#include <util/common/math/vec.h>
#include <util/common/plot/plot.h>
#include <util/common/math/fuzzy.h>

#include <vector>
#include <map>
#include <array>

#include <omp.h>

namespace model
{

    /*****************************************************/
    /*                     params                        */
    /*****************************************************/

    namespace consts
    {
        static const double eps    = 2.17;
        static const double alpha  = 2.59;
        static const double rho0   = 1.12;
        static const double rho1   = 1.8;
        static const double lambda = 21.0;
        static const double gamma  = 1.20;
        static const double A      = 7.05;
        static const double B      = 0.602;
        static const double m      = 28.086;
        static math::v3<> r[2][4] = {
            {
                math::v3<>{  1, -1, -1 } * alpha / 2 / 2,
                math::v3<>{ -1, -1,  1 } * alpha / 2 / 2,
                math::v3<>{ -1,  1, -1 } * alpha / 2 / 2,
                math::v3<>{  1,  1,  1 } * alpha / 2 / 2
            },
            {
                math::v3<>{  1,  1, -1 } * alpha / 2 / 2,
                math::v3<>{ -1,  1,  1 } * alpha / 2 / 2,
                math::v3<>{  1, -1,  1 } * alpha / 2 / 2,
                math::v3<>{ -1, -1, -1 } * alpha / 2 / 2
            }
        };
    };

    struct parameters
    {
        // system params
        size_t n;

        // other params
        double dt, dx, eps;
        bool freebc, wipeke, hardwipe;
    };

    inline parameters make_default_parameters()
    {
        parameters p =
        {
            // system params
            10,

            // other params
            1e-1, consts::rho0 / 100, 1e-9,
            false, true, false
        };
        return p;
    }

    /*****************************************************/
    /*                     data                          */
    /*****************************************************/

    struct particle
    {
        math::v3<> x, v0, v;
        std::vector < particle * > neighbors;
        bool outer = false;
    };

    using particle_fuzzy_t = math::fuzzy < math::fuzzy_weak_double_traits > ;

    inline bool sameloc(particle & p1, particle & p2)
    {
        return particle_fuzzy_t::eq(p1.x.x, p2.x.x) &&
               particle_fuzzy_t::eq(p1.x.y, p2.x.y) &&
               particle_fuzzy_t::eq(p1.x.z, p2.x.z);
    }

    class particles
    {
    public:
        std::vector < particle > all;
        std::vector < size_t > outer;
        std::vector < std::pair < size_t, size_t > > edges;
        double radius;
        bool skip_0 = false;
    public:
        void init(const parameters & p)
        {
            all.clear();
            outer.clear();
            edges.clear();
            
            std::vector < std::pair < particle, size_t > > prev_layer;
            std::vector < std::pair < particle, size_t > > cur_layer;
            std::vector < std::pair < particle, size_t > > next_layer;
            std::vector < std::pair < size_t, size_t > > cur_edges;
            std::vector < std::pair < size_t, size_t > > next_edges;

            cur_layer.emplace_back();

            radius = 0;

            for (size_t n = 0; n < p.n; ++n)
            {
                for (size_t k = 0; k < cur_layer.size(); ++k)
                {
                    auto & pt = cur_layer[k];
                    radius = std::fmax(radius, math::sqnorm(pt.first.x));
                    for (size_t i = 0; i < 4; ++i)
                    {
                        particle pt0 = pt.first; pt0.x = pt0.x + consts::r[n & 1][i];
                        bool skip = false;
                        #pragma omp parallel for reduction(|:skip)
                        for (int j = 0; j < prev_layer.size(); ++j)
                            if (sameloc(pt0, prev_layer[j].first)) { skip |= true; }
                        #pragma omp parallel for reduction(|:skip)
                        for (int j = 0; j < next_layer.size(); ++j)
                        {
                            if (sameloc(pt0, next_layer[j].first))
                            {
                                skip |= true;
                                next_edges.emplace_back(all.size(), j);
                            }
                        }
                        if (!skip)
                            next_layer.emplace_back(pt0, all.size());
                    }
                    all.push_back(pt.first);
                    if (n + 1 == p.n)
                    {
                        all.back().outer = true;
                        outer.push_back(all.size() - 1);
                    }
                    edges.emplace_back(all.size() - 1, pt.second);
                }
                for (size_t i = 0; i < cur_edges.size(); ++i)
                {
                    edges.emplace_back(cur_edges[i].first,
                                       all.size() - cur_layer.size() + cur_edges[i].second);
                }
                cur_edges.swap(next_edges);
                next_edges.clear();
                prev_layer.swap(cur_layer);
                cur_layer.swap(next_layer);
                next_layer.clear();
            }

            double thres2 = 1.5 * consts::rho1; thres2 *= thres2;

            #pragma omp parallel for
            for (int i = 0; i < all.size(); ++i)
            for (int j = 0; j < i; ++j)
            {
                if (i == j) continue;
                if (math::sqnorm(all[i].x - all[j].x) < thres2)
                {
                    #pragma omp critical
                    {
                        all[i].neighbors.push_back(&all[j]);
                        all[j].neighbors.push_back(&all[i]);
                    }
                }
            }

            radius = std::sqrt(radius) + consts::rho1;
        }

        double stillinger_weber_potential(
            size_t a, math::v3<> ax)
        {
            double e = 0;
            
            auto & ap = all[a];

            #pragma omp parallel for reduction(+:e)
            for (int i = 0; i < ap.neighbors.size(); ++i)
            {
                auto & ip = *ap.neighbors[i];
                if (skip_0 && (&all[0] == &ip)) continue;

                auto v_ai = ip.x - ax;
                double r_ai = math::norm(v_ai);

                double e0 = 0;

                e0 += f_ij(r_ai);
                
                for (size_t j = 0; j < i; ++j)
                {
                    auto & jp = *ap.neighbors[j];
                    auto v_ij = jp.x - ip.x;
                    auto v_ja = ax - jp.x;
                    double r_aj = math::norm(v_ja);
                    double r_ij = math::norm(v_ij);
                    double s_iaj = - (v_ai * v_ja);
                    double s_aij = - (v_ai * v_ij);
                    double s_ija = - (v_ij * v_ja);
                    e0 += h_ijk(r_ai, r_aj, s_iaj / r_ai / r_aj) +
                          h_ijk(r_ai, r_ij, s_aij / r_ai / r_ij) +
                          h_ijk(r_ij, r_aj, s_ija / r_ij / r_aj);
                }

                e += e0;
            }

            e *= consts::eps;

            return e;
        }

        double penergy()
        {
            double e = 0;
            for (size_t i = skip_0 ? 1 : 0; i < all.size(); ++i)
            {
                e += stillinger_weber_potential(i, all[i].x);
            }
            return e / 2;
        }

        double kenergy()
        {
            double e = 0;
            for (size_t i = skip_0 ? 1 : 0; i < all.size(); ++i)
            {
                e += math::sqnorm(all[i].v);
            }
            return consts::m * e / 2;
        }

        double f_ij(double r_ij)
        {
            if (r_ij > consts::rho1) return 0;
            return consts::A * (consts::B / r_ij / r_ij / r_ij / r_ij - 1) *
                        std::exp(1 / (r_ij - consts::rho1));
        }

        double h_ijk(double r_ij, double r_ik, double cos_jik)
        {
            if (r_ij > consts::rho1) return 0;
            if (r_ik > consts::rho1) return 0;
            return consts::lambda * std::exp(consts::gamma / (r_ij - consts::rho1) +
                                             consts::gamma / (r_ik - consts::rho1)) *
                    (cos_jik + 1 / 3.) * (cos_jik + 1 / 3.);
        }

        void next(const parameters & p)
        {
            for (size_t i = skip_0 ? 1 : 0; i < all.size(); ++i)
            {
                auto & a = all[i];
                if (!p.freebc & a.outer) continue;
                auto f0 = get_f(p, i);
                a.x = a.x + a.v * p.dt + f0 / 2 / consts::m * p.dt * p.dt;
                auto f1 = get_f(p, i);
                a.v = a.v + (f0 + f1) / 2 / consts::m * p.dt;
                if (p.wipeke && (math::sqnorm(a.v) < math::sqnorm(a.v0)))
                {
                    if (p.hardwipe)
                    {
                        #pragma omp parallel for
                        for (int j = skip_0 ? 1 : 0; j < all.size(); ++j)
                            all[j].v = {};
                    }
                    else
                    {
                        a.v = {};
                    }
                }
                a.v0 = a.v;
            }
        }

        math::v3<> get_f(const parameters & p, size_t i)
        {
            auto & a = all[i];
            auto U = [&] (const math::v3<> & d) { return stillinger_weber_potential(i, a.x + d); };
            return {
                - (U({ p.dx, 0, 0 }) - U({ -p.dx, 0, 0 })) / 2 / p.dx,
                - (U({ 0, p.dx, 0 }) - U({ 0, -p.dx, 0 })) / 2 / p.dx,
                - (U({ 0, 0, p.dx }) - U({ 0, 0, -p.dx })) / 2 / p.dx,
            };
        }
    };

    /*****************************************************/
    /*                     drawing                       */
    /*****************************************************/

    using points_t = std::vector < geom::point2d_t > ;

    struct plot_data
    {
        util::ptr_t < points_t > data;
        plot::list_drawable < points_t > :: ptr_t plot;
    };

    struct plot_config
    {
        plot::world_t::ptr_t world;
        plot::auto_viewport < points_t > :: ptr_t autoworld;
    };

    struct model_data
    {
        util::ptr_t < parameters > params;
        plot_config config;
        plot_data   penergy_data;
        plot_data   kenergy_data;
        plot_data   senergy_data;
        plot_data   denergy_data;
        particles   system_data;
    };

    inline static plot_data make_plot_data
    (
        plot::palette::pen_ptr pen = plot::palette::pen(0xffffff),
        plot::list_data_format data_format = plot::list_data_format::chain
    )
    {
        plot_data pd;
        pd.data = util::create < points_t > ();
        pd.plot = plot::list_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr, // no point painter
            pen
        );
        pd.plot->data_format = data_format;
        return pd;
    }

    inline static plot::drawable::ptr_t make_root_drawable
    (
        const plot_config & p,
        std::vector < plot::drawable::ptr_t > layers
    )
    {
        using namespace plot;

        return viewporter::create(
            tick_drawable::create(
                layer_drawable::create(layers),
                const_n_tick_factory<axe::x>::create(
                    make_simple_tick_formatter(6, 8),
                    0,
                    5
                ),
                const_n_tick_factory<axe::y>::create(
                    make_simple_tick_formatter(6, 8),
                    0,
                    5
                ),
                palette::pen(RGB(80, 80, 80)),
                RGB(200, 200, 200)
            ),
            make_viewport_mapper(make_world_mapper < points_t > (p.autoworld))
        );
    }

    inline plot_config make_plot_config()
    {
        plot_config cfg;
        cfg.world = plot::world_t::create();
        cfg.autoworld = plot::min_max_auto_viewport < points_t > :: create();
        return cfg;
    }

    inline model_data make_model_data(const parameters & p = make_default_parameters())
    {
        model_data md;
        md.config = make_plot_config();
        md.params = util::create < parameters > (p);
        md.kenergy_data = make_plot_data(plot::palette::pen(0x0000ff, 2));
        md.penergy_data = make_plot_data(plot::palette::pen(0x0000ff, 2));
        md.senergy_data = make_plot_data(plot::palette::pen(0x0000ff, 2));
        md.denergy_data = make_plot_data(plot::palette::pen(0x0000ff, 2));
        return md;
    }
}