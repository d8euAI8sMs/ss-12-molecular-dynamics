#pragma once

#include <vector>

#include <util/common/geom/point.h>
#include <util/common/plot/shape.h>
#include <util/common/gui/OglControl.h>

#include "model.h"

class CGridPlotControl : public COglControl
{
    DECLARE_DYNAMIC(CGridPlotControl)

public:
    model::model_data * m_data;
public:
    CGridPlotControl();
    virtual ~CGridPlotControl();
    virtual void OnDrawItemOGL();

protected:
    DECLARE_MESSAGE_MAP()
};
