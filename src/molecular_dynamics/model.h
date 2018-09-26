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
        static const double m     = 28.086;
        static const double H0    = 4.67;
        static const double alpha = 0.543;
        static const double rho0  = 0.235;
        static const double b     = 0.0956;
        static const double d     = 0.0239;
        static const double c     = 0.192;
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
            1e-2, consts::rho0 / 100, 1e-9,
            false, true, true
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
        std::map < size_t, std::vector < size_t > > edge_neighbors;
        double radius;
    public:
        void init(const parameters & p, bool skip_0)
        {
            all.clear();
            outer.clear();
            edges.clear();
            edge_neighbors.clear();
            
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

            for (auto it = edges.begin(); it != edges.end();)
            {
                if ((it->first == it->second) || skip_0 && (it->first == 0 || it->second == 0))
                {
                    it = edges.erase(it);
                }
                else
                {
                    ++it;
                }
            }

            if (skip_0)
            {
                edges.emplace_back(1, 2);
                edges.emplace_back(3, 4);
            }

            for (int i = 0; i < edges.size(); ++i)
            {
                auto & e = edges[i];
                all[e.first].neighbors.push_back(&all[e.second]);
                all[e.second].neighbors.push_back(&all[e.first]);
            }
            
            std::vector < std::pair < particle *, particle * > > coincidents;
            for (size_t i = 0; i < edges.size(); ++i)
            {
                auto & e = edges[i];
                auto & p1 = all[e.first];
                auto & p2 = all[e.second];
                coincidents.clear();
                for (size_t j = 0; j < p1.neighbors.size(); ++j)
                {
                    if (p1.neighbors[j] != &p2)
                    {
                        coincidents.emplace_back(p1.neighbors[j], &p1);
                    }
                }
                for (size_t j = 0; j < p2.neighbors.size(); ++j)
                {
                    if (p2.neighbors[j] != &p1)
                    {
                        coincidents.emplace_back(p2.neighbors[j], &p2);
                    }
                }
                for (size_t j = 0; j < coincidents.size(); ++j)
                {
                    auto & c = coincidents[j];
                    size_t i1 = 0, i2 = 0;
                    for (size_t k = 0; k < all.size(); ++k)
                    {
                        if (&all[k] == c.first) { i1 = k; }
                        if (&all[k] == c.second) { i2 = k; }
                    }
                    for (size_t k = 0; k < edges.size(); ++k)
                    {
                        if (i1 == edges[k].first && i2 == edges[k].second ||
                            i2 == edges[k].first && i1 == edges[k].second)
                        {
                            edge_neighbors[i].push_back(k);
                            break;
                        }
                    }
                }
            }

            radius = std::sqrt(radius);
        }

        double vukcevic_potential(size_t a, size_t ae, math::v3<> ax)
        {
            double e = 0;

            if (all[a].outer || all[edges[ae].first].outer || all[edges[ae].second].outer) return 0;
            
            math::v3<> c1, c2;
            if (a == edges[ae].first) c1 = ax - all[edges[ae].second].x;
            else                      c1 = all[edges[ae].first].x - ax;
            double c1n2 = math::sqnorm(c1);
            double c1n = std::sqrt(c1n2);

            for (size_t i = 0; i < edge_neighbors[ae].size(); ++i)
            {
                if (a == edges[edge_neighbors[ae][i]].first)
                    c2 = ax - all[edges[edge_neighbors[ae][i]].second].x;
                else if (a == edges[edge_neighbors[ae][i]].second)
                    c2 = all[edges[edge_neighbors[ae][i]].first].x - ax;
                else
                    c2 = (all[edges[edge_neighbors[ae][i]].first].x - all[edges[edge_neighbors[ae][i]].second].x);
                auto cos2 = (c1 * c2) / c1n2 / math::sqnorm(c2) * (c1 * c2);
                e += std::exp(-(cos2 - 1/9.) * (cos2 - 1/9.) / consts::c);
            }

            e /= edge_neighbors[ae].size();
            e *= consts::H0 / 2 / (consts::b - consts::d);
            e *= (consts::d * std::exp((consts::rho0 - c1n) / consts::d) - consts::b * std::exp((consts::rho0 - c1n) / consts::b));
            e += consts::H0 / 2;

            return e;
        }

        double vukcevic_potential(size_t a, math::v3<> ax)
        {
            double e = 0;

            for (size_t j = 0; j < edges.size(); ++j)
            {
                if (edges[j].first == a || edges[j].second == a)
                    e += vukcevic_potential(a, j, ax);
            }

            return e;
        }

        double penergy()
        {
            double e = 0;
            for (size_t i = 0; i < all.size(); ++i)
            {
                e += vukcevic_potential(i, all[i].x);
            }
            return e;
        }

        double penergy_eq()
        {
            double e = 0;
            for (size_t j = 0; j < edges.size(); ++j)
            {
                if (all[edges[j].first].outer || all[edges[j].second].outer) continue;
                e -= consts::H0;
            }
            return e;
        }

        double kenergy()
        {
            double e = 0;
            for (size_t i = 0; i < all.size(); ++i)
            {
                e += math::sqnorm(all[i].v);
            }
            return consts::m * e / 2;
        }

        void next(const parameters & p)
        {
            for (size_t i = 0; i < all.size(); ++i)
            {
                auto & a = all[i];
                if (!p.freebc & a.outer) continue;
                bool skip = false;
                if (!p.freebc)
                {
                    for (size_t j = 0; j < a.neighbors.size(); ++j)
                    {
                        if (a.neighbors[j]->outer) { skip = true; break; }
                    }
                    if (skip) continue;
                }
                auto f0 = get_f(p, i);
                a.x = a.x + a.v * p.dt + f0 / 2 / consts::m * p.dt * p.dt;
                auto f1 = get_f(p, i);
                a.v = a.v + (f0 + f1) / 2 / consts::m * p.dt;
                if (p.wipeke && (math::sqnorm(a.v) < math::sqnorm(a.v0)))
                {
                    if (p.hardwipe)
                    {
                        #pragma omp parallel for
                        for (int j = 0; j < all.size(); ++j)
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
            auto U = [&] (const math::v3<> & d) { return vukcevic_potential(i, a.x + d); };
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