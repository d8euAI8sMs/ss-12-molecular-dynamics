
// molecular_dynamics.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
    #error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"        // main symbols


// CMolecularDynamicsApp:
// See molecular_dynamics.cpp for the implementation of this class
//

class CMolecularDynamicsApp : public CWinApp
{
public:
    CMolecularDynamicsApp();

// Overrides
public:
    virtual BOOL InitInstance();

// Implementation

    DECLARE_MESSAGE_MAP()
};

extern CMolecularDynamicsApp theApp;