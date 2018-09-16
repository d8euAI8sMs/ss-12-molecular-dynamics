
// molecular_dynamicsDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/gui/PlotControl.h>
#include "afxwin.h"

#include "model.h"
#include "GridPlotControl.h"

// CMolecularDynamicsDlg dialog
class CMolecularDynamicsDlg : public CSimulationDialog
{
// Construction
public:
    CMolecularDynamicsDlg(CWnd* pParent = NULL);    // standard constructor

// Dialog Data
    enum { IDD = IDD_MOLECULAR_DYNAMICS_DIALOG };

    protected:
    virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support


// Implementation
protected:
    HICON m_hIcon;

    // Generated message map functions
    virtual BOOL OnInitDialog();
    afx_msg void OnPaint();
    afx_msg HCURSOR OnQueryDragIcon();
    DECLARE_MESSAGE_MAP()
public:
    virtual void OnSimulation();
    CGridPlotControl m_gridCtrl;
    CPlotControl m_energyCtrl;
    model::model_data m_data;
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton2();
    CEdit m_initialEnergy;
    CEdit m_vacancyEnergy;
    BOOL m_withVacancy;
    CButton m_radioEp;
    CButton m_radioEk;
    CButton m_radioEs;
    CButton m_radioEd;
    afx_msg void OnUpdateEnergyPlot();
    BOOL m_bFreeBC;
    BOOL m_bWipeEk;
    BOOL m_bHardWipe;
    CEdit m_sumEnergyFluct;
};
