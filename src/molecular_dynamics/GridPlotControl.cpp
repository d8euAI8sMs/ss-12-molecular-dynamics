// SurfacePlotControl.cpp : implementation file
//

#include "stdafx.h"
#include "GridPlotControl.h"

#include <gl/GLU.h>
#include <gl/GL.h>

// CSurfacePlotControl

IMPLEMENT_DYNAMIC(CGridPlotControl, COglControl)

CGridPlotControl::CGridPlotControl()
    : m_data(nullptr)
{
}

CGridPlotControl::~CGridPlotControl()
{
}


BEGIN_MESSAGE_MAP(CGridPlotControl, COglControl)
END_MESSAGE_MAP()



// CSurfacePlotControl message handlers

void CGridPlotControl::OnDrawItemOGL()
{
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glLineWidth(5);

    m_zoomMultiplier = model::consts::rho0 / 10;
    m_translateXMultiplier = model::consts::alpha;
    m_translateYMultiplier = model::consts::alpha;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScaled(1 / m_data->system_data.radius,
             1 / m_data->system_data.radius,
             1 / m_data->system_data.radius);

	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-2, 2, -2, 2, -2, 2);
    glTranslated(GetTranslateX(), GetTranslateY(), -model::consts::rho0);
    glRotated(GetRotAngleH() / M_PI * 180, 1, 0, 0);
    glRotated(GetRotAngleV() / M_PI * 180, 0, 1, 0);
    glScaled(GetZoom(), GetZoom(), GetZoom());

    glMatrixMode(GL_MODELVIEW);

    glClearColor(0, 0, 0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (m_data->system_data.all.empty()) return;

    GLUquadric * q = gluNewQuadric();
    gluQuadricDrawStyle(q, GLU_FILL);

    for (size_t i = 1; i < m_data->system_data.all.size(); ++i)
    {
        auto & e = m_data->system_data.all[i];
        if (e.outer) continue;
        glTranslated(e.x.x, e.x.y, e.x.z);
        glColor3d(1, 0, 0); gluSphere(q, model::consts::rho0 / 4, 10, 10);
        glTranslated(-e.x.x, -e.x.y, -e.x.z);
    }
    
    {
        auto & e = m_data->system_data.all[0];
        glTranslated(e.x.x, e.x.y, e.x.z);
        glColor3d(0, 0, 1); gluSphere(q, model::consts::rho0 / 4, 10, 10);
        glTranslated(-e.x.x, -e.x.y, -e.x.z);
    }

    for (size_t i = 0; i < m_data->system_data.outer.size(); ++i)
    {
        auto & e = m_data->system_data.all[m_data->system_data.outer[i]];
        glTranslated(e.x.x, e.x.y, e.x.z);
        glColor3d(0, 0, 1); gluSphere(q, model::consts::rho0 / 4, 10, 10);
        glTranslated(-e.x.x, -e.x.y, -e.x.z);
    }

    gluDeleteQuadric(q);
    
    glBegin(GL_LINES);
    for (size_t i = 0; i < m_data->system_data.edges.size(); ++i)
    {
        auto & e = m_data->system_data.edges[i];
        auto & v1 = m_data->system_data.all[e.first].x;
        auto & v2 = m_data->system_data.all[e.second].x;
        glColor3d(0, 1, 0); glVertex3d(v1.x, v1.y, v1.z); glVertex3d(v2.x, v2.y, v2.z);
    }
    glEnd();

    glFinish();
}
