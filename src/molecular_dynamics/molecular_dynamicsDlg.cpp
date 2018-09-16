
// molecular_dynamicsDlg.cpp : implementation file
//

#include "stdafx.h"
#include "molecular_dynamics.h"
#include "molecular_dynamicsDlg.h"
#include "afxdialogex.h"

#include <omp.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMolecularDynamicsDlg dialog



CMolecularDynamicsDlg::CMolecularDynamicsDlg(CWnd* pParent /*=NULL*/)
    : CSimulationDialog(CMolecularDynamicsDlg::IDD, pParent)
    , m_data(model::make_model_data())
    , m_withVacancy(m_data.system_data.skip_0)
    , m_bFreeBC(m_data.params->freebc)
    , m_bWipeEk(m_data.params->wipeke)
    , m_bHardWipe(m_data.params->hardwipe)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMolecularDynamicsDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_STRUCT, m_gridCtrl);
    DDX_Control(pDX, IDC_ENERGY, m_energyCtrl);
    DDX_Control(pDX, IDC_EDIT4, m_initialEnergy);
    DDX_Control(pDX, IDC_EDIT5, m_vacancyEnergy);
    DDX_Text(pDX, IDC_EDIT1, m_data.params->dt);
    DDX_Text(pDX, IDC_EDIT2, m_data.params->eps);
    DDX_Text(pDX, IDC_EDIT3, m_data.params->n);
    DDX_Check(pDX, IDC_CHECK1, m_withVacancy);
    DDX_Control(pDX, IDC_RADIO1, m_radioEp);
    DDX_Control(pDX, IDC_RADIO2, m_radioEk);
    DDX_Control(pDX, IDC_RADIO3, m_radioEs);
    DDX_Control(pDX, IDC_RADIO4, m_radioEd);
    DDX_Check(pDX, IDC_CHECK2, m_bFreeBC);
    DDX_Check(pDX, IDC_CHECK3, m_bWipeEk);
    DDX_Check(pDX, IDC_CHECK4, m_bHardWipe);
}

BEGIN_MESSAGE_MAP(CMolecularDynamicsDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CMolecularDynamicsDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON2, &CMolecularDynamicsDlg::OnBnClickedButton2)
    ON_BN_CLICKED(IDC_RADIO1, &CMolecularDynamicsDlg::OnUpdateEnergyPlot)
    ON_BN_CLICKED(IDC_RADIO2, &CMolecularDynamicsDlg::OnUpdateEnergyPlot)
    ON_BN_CLICKED(IDC_RADIO3, &CMolecularDynamicsDlg::OnUpdateEnergyPlot)
    ON_BN_CLICKED(IDC_RADIO4, &CMolecularDynamicsDlg::OnUpdateEnergyPlot)
END_MESSAGE_MAP()


// CMolecularDynamicsDlg message handlers

BOOL CMolecularDynamicsDlg::OnInitDialog()
{
    CSimulationDialog::OnInitDialog();

    // Set the icon for this dialog.  The framework does this automatically
    //  when the application's main window is not a dialog
    SetIcon(m_hIcon, TRUE);            // Set big icon
    SetIcon(m_hIcon, FALSE);        // Set small icon

    m_gridCtrl.m_data = &m_data;
    m_gridCtrl.RedrawWindow();

    m_radioEp.SetCheck(BST_UNCHECKED);
    m_radioEk.SetCheck(BST_UNCHECKED);
    m_radioEs.SetCheck(BST_UNCHECKED);
    m_radioEd.SetCheck(BST_CHECKED);

    m_data.penergy_data.plot->visible = false;
    m_data.kenergy_data.plot->visible = false;
    m_data.senergy_data.plot->visible = false;
    m_data.denergy_data.plot->visible = true;

    m_energyCtrl.triple_buffered = true;
    m_energyCtrl.plot_layer.with(
        model::make_root_drawable(m_data.config, {{
            m_data.penergy_data.plot,
            m_data.kenergy_data.plot,
            m_data.senergy_data.plot,
            m_data.denergy_data.plot
        }})
    );

    // TODO: Add extra initialization here

    return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CMolecularDynamicsDlg::OnPaint()
{
    if (IsIconic())
    {
        CPaintDC dc(this); // device context for painting

        SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

        // Center icon in client rectangle
        int cxIcon = GetSystemMetrics(SM_CXICON);
        int cyIcon = GetSystemMetrics(SM_CYICON);
        CRect rect;
        GetClientRect(&rect);
        int x = (rect.Width() - cxIcon + 1) / 2;
        int y = (rect.Height() - cyIcon + 1) / 2;

        // Draw the icon
        dc.DrawIcon(x, y, m_hIcon);
    }
    else
    {
        CSimulationDialog::OnPaint();
    }
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CMolecularDynamicsDlg::OnQueryDragIcon()
{
    return static_cast<HCURSOR>(m_hIcon);
}

void CMolecularDynamicsDlg::OnSimulation()
{
    omp_set_num_threads(omp_get_num_procs() / 2);

    m_data.senergy_data.data->clear();
    m_data.kenergy_data.data->clear();
    m_data.penergy_data.data->clear();
    m_data.denergy_data.data->clear();
    m_data.config.autoworld->clear();

    CString fmt;

    auto & s = m_data.system_data;
    s.skip_0 = false;
    s.init(*m_data.params);
    if (m_withVacancy) s.skip_0 = true;
    m_data.params->freebc = m_bFreeBC == TRUE;
    m_data.params->wipeke = m_bWipeEk == TRUE;
    m_data.params->hardwipe = m_bHardWipe == TRUE;

    m_gridCtrl.RedrawWindow();
    
    double e0 = s.penergy();
    if (m_withVacancy) e0 *= (s.all.size() - 1) / (double) s.all.size();

    fmt.Format(TEXT("%lf"), e0);
    m_initialEnergy.SetWindowText(fmt);

    double t = 0;

    double e1 = e0;

    while (m_bWorking)
    {
        s.next(*m_data.params);
        
        double e2 = s.penergy();
        fmt.Format(TEXT("%lf"), std::abs(e2 - e0));
        m_vacancyEnergy.SetWindowText(fmt);

        double e3 = s.kenergy();

        m_data.penergy_data.data->push_back({ t, e2 });
        m_data.kenergy_data.data->push_back({ t, e3 });
        m_data.senergy_data.data->push_back({ t, e3 + e2 });
        m_data.denergy_data.data->push_back({ t, std::abs(e0 - e2) });

        OnUpdateEnergyPlot();
        
        if (std::abs(e2 - e1) < m_data.params->eps) break;

        e1 = e2;
        t += m_data.params->dt;
    }

    CSimulationDialog::OnSimulation();
}


void CMolecularDynamicsDlg::OnBnClickedButton1()
{
    UpdateData(TRUE);
    StartSimulationThread();
}

void CMolecularDynamicsDlg::OnBnClickedButton2()
{
    StopSimulationThread();
}


void CMolecularDynamicsDlg::OnUpdateEnergyPlot()
{
    m_data.penergy_data.plot->visible = m_radioEp.GetCheck() == BST_CHECKED;
    m_data.kenergy_data.plot->visible = m_radioEk.GetCheck() == BST_CHECKED;
    m_data.senergy_data.plot->visible = m_radioEs.GetCheck() == BST_CHECKED;
    m_data.denergy_data.plot->visible = m_radioEd.GetCheck() == BST_CHECKED;
    m_data.config.autoworld->clear();
    if (m_data.penergy_data.plot->visible)
        m_data.config.autoworld->setup(*m_data.penergy_data.data);
    if (m_data.kenergy_data.plot->visible)
        m_data.config.autoworld->setup(*m_data.kenergy_data.data);
    if (m_data.senergy_data.plot->visible)
        m_data.config.autoworld->setup(*m_data.senergy_data.data);
    if (m_data.denergy_data.plot->visible)
        m_data.config.autoworld->setup(*m_data.denergy_data.data);
    m_energyCtrl.RedrawBuffer();
    m_energyCtrl.SwapBuffers();
    Invoke([this] () { m_energyCtrl.RedrawWindow(); });
}
