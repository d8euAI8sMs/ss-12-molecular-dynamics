
// molecular_dynamicsDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>

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
};
