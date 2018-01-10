
// visualizer101Dlg.cpp : implementation file
//

#include "stdafx.h"
#include "visualizer101.h"
#include "visualizer101Dlg.h"
#include "afxdialogex.h"
#include "Visualizer.h"
#include <string>

using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// Cvisualizer101Dlg dialog



Cvisualizer101Dlg::Cvisualizer101Dlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(IDD_visualizer101_DIALOG, pParent)
	, channelStr(_T(""))
	, widthStr(_T(""))
	, heightStr(_T(""))
	, multisamplingStr(_T(""))
	, orderVal(0)
	, uincVal(0)
	, xincVal(0)
	, skySplitsVal(0)
	, srevSplitsVal(0)
	, tiVal(0)
	, toVal(0)
	, aaVal(FALSE)
	, peakRotateVal(FALSE)
	, peakFovVal(FALSE)
	, animTessVal(FALSE)
	, volumeVal(0)
	, lightingVal(FALSE)
	, explodeVal(FALSE)
	, wireframeVal(FALSE)
	, simpleNormsVal(FALSE)
	, switchNormsVal(FALSE)
	, texSkyVal(FALSE)
	, colourStr(_T(""))
	, shapreStr(_T(""))
	, deviceStr(_T(""))
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void Cvisualizer101Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_COMBO1, channelBox);
	DDX_CBString(pDX, IDC_COMBO1, channelStr);
	DDV_MaxChars(pDX, channelStr, 5);
	DDX_Control(pDX, IDC_EDIT1, widthBox);
	DDX_Text(pDX, IDC_EDIT1, widthStr);
	DDV_MaxChars(pDX, widthStr, 5);
	DDX_Control(pDX, IDC_EDIT2, heightBox);
	DDX_Text(pDX, IDC_EDIT2, heightStr);
	DDV_MaxChars(pDX, heightStr, 5);
	DDX_Control(pDX, IDC_COMBO2, multisamplingBox);
	DDX_CBString(pDX, IDC_COMBO2, multisamplingStr);
	DDV_MaxChars(pDX, multisamplingStr, 5);
	DDX_Control(pDX, IDC_SLIDER1, orderSlider);
	DDX_Slider(pDX, IDC_SLIDER1, orderVal);
	DDV_MinMaxInt(pDX, orderVal, 1, 10);
	DDX_Control(pDX, IDC_STATIC_TXT, orderLabel);
	DDX_Control(pDX, IDC_STATIC_TXT2, uincLabel);
	DDX_Control(pDX, IDC_STATIC_TXT3, xincLabel);
	DDX_Control(pDX, IDC_SLIDER2, uincSlider);
	DDX_Control(pDX, IDC_SLIDER6, xincSlider);
	DDX_Slider(pDX, IDC_SLIDER2, uincVal);
	DDV_MinMaxInt(pDX, uincVal, 1, 500);
	DDX_Slider(pDX, IDC_SLIDER6, xincVal);
	DDV_MinMaxInt(pDX, xincVal, 1, 100);
	DDX_Control(pDX, IDC_COMBO3, sampleRateBox);
	DDX_Control(pDX, IDC_COMBO4, numChannelsBox);
	DDX_Control(pDX, IDC_STATIC_TXT4, skySplitsLabel);
	DDX_Control(pDX, IDC_STATIC_TXT5, srevSplitsLabel);
	DDX_Control(pDX, IDC_STATIC_TXT6, tiLabel);
	DDX_Control(pDX, IDC_STATIC_TXT7, toLabel);
	DDX_Control(pDX, IDC_SLIDER7, skySplitsSlider);
	DDX_Control(pDX, IDC_SLIDER8, srevSplitsSlider);
	DDX_Control(pDX, IDC_SLIDER9, tiSlider);
	DDX_Control(pDX, IDC_SLIDER10, toSlider);
	DDX_Slider(pDX, IDC_SLIDER7, skySplitsVal);
	DDX_Slider(pDX, IDC_SLIDER8, srevSplitsVal);
	DDX_Slider(pDX, IDC_SLIDER9, tiVal);
	DDX_Slider(pDX, IDC_SLIDER10, toVal);
	DDX_Control(pDX, IDC_CHECK1, aaBox);
	DDX_Check(pDX, IDC_CHECK1, aaVal);
	DDX_Control(pDX, IDC_CHECK4, peakRotateBox);
	DDX_Check(pDX, IDC_CHECK4, peakRotateVal);
	DDX_Control(pDX, IDC_CHECK5, peakFovBox);
	DDX_Check(pDX, IDC_CHECK5, peakFovVal);
	DDX_Control(pDX, IDC_CHECK6, animTessBox);
	DDX_Check(pDX, IDC_CHECK6, animTessVal);
	DDX_Control(pDX, IDC_SLIDER11, volumeSlider);
	DDX_Slider(pDX, IDC_SLIDER11, volumeVal);
	DDX_Control(pDX, IDC_CHECK2, lightingBox);
	DDX_Check(pDX, IDC_CHECK2, lightingVal);
	DDX_Control(pDX, IDC_CHECK3, explodeBox);
	DDX_Check(pDX, IDC_CHECK3, explodeVal);
	DDX_Control(pDX, IDC_CHECK7, wireframeBox);
	DDX_Check(pDX, IDC_CHECK7, wireframeVal);
	DDX_Control(pDX, IDC_CHECK8, simpleNormsBox);
	DDX_Check(pDX, IDC_CHECK8, simpleNormsVal);
	DDX_Control(pDX, IDC_CHECK9, switchNormsBox);
	DDX_Check(pDX, IDC_CHECK9, switchNormsVal);
	DDX_Control(pDX, IDC_CHECK10, texSkyBox);
	DDX_Check(pDX, IDC_CHECK10, texSkyVal);
	DDX_Control(pDX, IDC_COMBO5, colourBox);
	DDX_Control(pDX, IDC_COMBO6, shapeBox);
	DDX_CBString(pDX, IDC_COMBO5, colourStr);
	DDX_CBString(pDX, IDC_COMBO6, shapreStr);
	DDX_Control(pDX, IDC_COMBO7, deviceBox);
	DDX_CBString(pDX, IDC_COMBO7, deviceStr);
}

BEGIN_MESSAGE_MAP(Cvisualizer101Dlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDOK, &Cvisualizer101Dlg::OnBnClickedOk)
	ON_CBN_SELCHANGE(IDC_COMBO1, &Cvisualizer101Dlg::OnCbnSelchangeCombo1)
	ON_EN_CHANGE(IDC_EDIT1, &Cvisualizer101Dlg::OnEnChangeEdit1)
	ON_EN_CHANGE(IDC_EDIT2, &Cvisualizer101Dlg::OnEnChangeEdit2)
	ON_CBN_SELCHANGE(IDC_COMBO2, &Cvisualizer101Dlg::OnCbnSelchangeCombo2)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER1, &Cvisualizer101Dlg::OnNMCustomdrawSlider1)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER2, &Cvisualizer101Dlg::OnNMCustomdrawSlider2)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER6, &Cvisualizer101Dlg::OnNMCustomdrawSlider6)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER7, &Cvisualizer101Dlg::OnNMCustomdrawSlider7)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER8, &Cvisualizer101Dlg::OnNMCustomdrawSlider8)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER9, &Cvisualizer101Dlg::OnNMCustomdrawSlider9)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER10, &Cvisualizer101Dlg::OnNMCustomdrawSlider10)
	ON_BN_CLICKED(IDC_CHECK1, &Cvisualizer101Dlg::OnBnClickedCheck1)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER11, &Cvisualizer101Dlg::OnNMCustomdrawSlider11)
	ON_BN_CLICKED(IDC_CHECK4, &Cvisualizer101Dlg::OnBnClickedCheck4)
	ON_BN_CLICKED(IDC_CHECK5, &Cvisualizer101Dlg::OnBnClickedCheck5)
	ON_BN_CLICKED(IDC_CHECK6, &Cvisualizer101Dlg::OnBnClickedCheck6)
	ON_BN_CLICKED(IDC_CHECK2, &Cvisualizer101Dlg::OnBnClickedCheck2)
	ON_BN_CLICKED(IDC_CHECK3, &Cvisualizer101Dlg::OnBnClickedCheck3)
	ON_BN_CLICKED(IDC_CHECK7, &Cvisualizer101Dlg::OnBnClickedCheck7)
	ON_BN_CLICKED(IDC_CHECK8, &Cvisualizer101Dlg::OnBnClickedCheck8)
	ON_BN_CLICKED(IDC_CHECK9, &Cvisualizer101Dlg::OnBnClickedCheck9)
	ON_BN_CLICKED(IDC_CHECK10, &Cvisualizer101Dlg::OnBnClickedCheck10)
	ON_CBN_SELCHANGE(IDC_COMBO6, &Cvisualizer101Dlg::OnCbnSelchangeCombo6)
	ON_CBN_SELCHANGE(IDC_COMBO5, &Cvisualizer101Dlg::OnCbnSelchangeCombo5)
	ON_CBN_SELCHANGE(IDC_COMBO4, &Cvisualizer101Dlg::OnCbnSelchangeCombo4)
	ON_CBN_SELCHANGE(IDC_COMBO3, &Cvisualizer101Dlg::OnCbnSelchangeCombo3)
	ON_CBN_SELCHANGE(IDC_COMBO7, &Cvisualizer101Dlg::OnCbnSelchangeCombo7)
END_MESSAGE_MAP()


// Cvisualizer101Dlg message handlers

BOOL Cvisualizer101Dlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// audio section
	sampleRateBox.SetCurSel(0);
	numChannelsBox.SetCurSel(1);
	channelBox.SetCurSel(3);
	channelStr = L"6";
	deviceBox.SetCurSel(0);
	deviceStr = L"0";
	volumeSlider.SetRange(75, 600, TRUE);
	volumeSlider.SetPos(325);
	volumeVal = 325;

	// window section
	widthBox.SetWindowText(L"1920");
	heightBox.SetWindowText(L"1080");
	multisamplingBox.SetCurSel(4);
	multisamplingStr = L"8";
	aaBox.SetCheck(1);
	aaVal = TRUE;

	// bspline section
	orderSlider.SetRange(1, 9, TRUE);
	orderSlider.SetPos(3);
	orderVal = 3;
	uincSlider.SetRange(1, 200, TRUE);
	uincSlider.SetPos(5);
	uincVal = 5;
	xincSlider.SetRange(1, 20, TRUE);
	xincSlider.SetPos(10);
	xincVal = 10;

	// resolution/subdivision section
	skySplitsSlider.SetRange(1, 50, TRUE);
	skySplitsSlider.SetPos(12);
	skySplitsVal = 12;
	srevSplitsSlider.SetRange(1, 50, TRUE);
	srevSplitsSlider.SetPos(24);
	srevSplitsVal = 24;
	tiSlider.SetRange(1, 12, TRUE);
	tiSlider.SetPos(1);
	tiVal = 1;
	toSlider.SetRange(1, 12, TRUE);
	toSlider.SetPos(1);
	toVal = 1;

	// animations section
	peakRotateBox.SetCheck(1);
	peakRotateVal = TRUE;
	peakFovBox.SetCheck(1);
	peakFovVal = TRUE;
	animTessBox.SetCheck(0);
	animTessVal = FALSE;

	// display section
	lightingBox.SetCheck(1);
	lightingVal = TRUE;
	explodeBox.SetCheck(0);
	explodeVal = FALSE;
	wireframeBox.SetCheck(0);
	wireframeVal = FALSE;
	simpleNormsBox.SetCheck(1);
	simpleNormsVal = TRUE;
	switchNormsBox.SetCheck(1);
	switchNormsVal = TRUE;
	texSkyBox.SetCheck(1);
	texSkyVal = TRUE;
	colourBox.SetCurSel(1);
	colourStr = L"Rainbow";
	shapeBox.SetCurSel(2);
	shapreStr = L"Flat";

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void Cvisualizer101Dlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void Cvisualizer101Dlg::OnPaint()
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
		CDialogEx::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR Cvisualizer101Dlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void Cvisualizer101Dlg::OnBnClickedOk()
{
	// set window values
	int width = _wtoi(widthStr);
	int height = _wtoi(heightStr);
	int samples = _wtoi(multisamplingStr);
	SetGLWindow(width, height, samples);

	// set audio input values
	int chan = _wtoi(channelStr);
	int dev = _wtoi(deviceStr);
	SetAudioInput(44100, 2, chan, dev);

	// set bspline values
	float uinc = float(uincVal) / 1000;
	SetBSpline(orderVal, uinc);
	
	// set display values
	bool draw = false;
	float xinc = float(xincVal) / 100;
	float yfactor = float(volumeVal) / 100;
	int colours = 0, shape = 0;
	if (colourStr == "Rainbow") colours = 1;
	else if (colourStr == "Anim1") colours = 2;
	else if (colourStr == "Anim2") colours = 3;
	else if (colourStr == "Anim3") colours = 4;
	if (shapreStr == "Pointy") shape = 1;
	else if (shapreStr == "Pointy Big") shape = 2;
	else if (shapreStr == "Multiple") shape = 3;
	else if (shapreStr == "Moving") shape = 4;
	else if (shapreStr == "Jewels") shape = 5;
	else if (shapreStr == "Ellipsoid") shape = 6;
	else if (shapreStr == "Spheroid") shape = 7;
	else if (shapreStr == "CrazySphere") shape = 8;
	else if (shapreStr == "Bell") shape = 9;
	else if (shapreStr == "DRAW") {
		shape = -1;
		draw = true;
	}
	else if (shapreStr == "POINTS") {
		shape = -2;
		draw = true;
	}

	SetDisplayValues(yfactor, xinc, srevSplitsVal, skySplitsVal, 
						tiVal, toVal, colours, shape);
	
	// set flags
	SetFlags(draw, animTessVal, aaVal, explodeVal, lightingVal, 
				peakFovVal, peakRotateVal, !simpleNormsVal, 
				switchNormsVal, texSkyVal, !texSkyVal, wireframeVal);

	if (draw) {
		if (shape == -2) {
			TCHAR e[1000] = _T("Click to add or move points, right click to delete them, and press enter when finished");
			MessageBox(e, _T("Control point draw mode"), NULL);
		}
		else {
			TCHAR e[1000] = _T("Click to draw and press enter when finished");
			MessageBox(e, _T("Freehand draw mode"), NULL);
		}
	}

	// clear data arrays
	ClearPoints();
	ClearVertices(0);
	ClearIndices(0);

	// run
	RunVisualizer();
	CDialogEx::OnOK();
}

void Cvisualizer101Dlg::OnBnClickedButton1()
{
	CWnd *label = GetDlgItem(IDC_STATIC);
	label->SetWindowText(channelStr);
}

void Cvisualizer101Dlg::OnCbnSelchangeCombo1()
{
	channelBox.GetLBText(channelBox.GetCurSel(), channelStr);
}

void Cvisualizer101Dlg::OnEnChangeEdit1()
{
	widthBox.GetWindowText(widthStr);
}

void Cvisualizer101Dlg::OnEnChangeEdit2()
{
	heightBox.GetWindowText(heightStr);
}

void Cvisualizer101Dlg::OnCbnSelchangeCombo2()
{
	multisamplingBox.GetLBText(multisamplingBox.GetCurSel(), multisamplingStr);
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	orderVal = orderSlider.GetPos();
	string orderStr = to_string(orderVal);
	wstring orderWStr (orderStr.begin(), orderStr.end());
	orderLabel.SetWindowText(orderWStr.c_str());
	*pResult = 0;
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider2(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	uincVal = uincSlider.GetPos();
	string uincStr = to_string(uincVal);
	wstring uincWStr(uincStr.begin(), uincStr.end());
	uincWStr.append(L"/1000");
	uincLabel.SetWindowText(uincWStr.c_str());
	*pResult = 0;
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider6(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	xincVal = xincSlider.GetPos();
	string xincStr = to_string(xincVal);
	wstring xincWStr(xincStr.begin(), xincStr.end());
	xincWStr.append(L"/100");
	xincLabel.SetWindowText(xincWStr.c_str());
	*pResult = 0;
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider7(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	skySplitsVal = skySplitsSlider.GetPos();
	string skySplitsStr = to_string(skySplitsVal);
	wstring skySplitsWStr(skySplitsStr.begin(), skySplitsStr.end());
	skySplitsLabel.SetWindowText(skySplitsWStr.c_str());
	*pResult = 0;
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider8(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	srevSplitsVal = srevSplitsSlider.GetPos();
	string srevSplitsStr = to_string(srevSplitsVal);
	wstring srevSplitsWStr(srevSplitsStr.begin(), srevSplitsStr.end());
	srevSplitsLabel.SetWindowText(srevSplitsWStr.c_str());
	*pResult = 0;
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider9(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	tiVal = tiSlider.GetPos();
	string tiStr = to_string(tiVal);
	wstring tiWStr(tiStr.begin(), tiStr.end());
	tiLabel.SetWindowText(tiWStr.c_str());
	*pResult = 0;
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider10(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	toVal = toSlider.GetPos();
	string toStr = to_string(toVal);
	wstring toWStr(toStr.begin(), toStr.end());
	toLabel.SetWindowText(toWStr.c_str());
	*pResult = 0;
}

void Cvisualizer101Dlg::OnBnClickedCheck1()
{
	aaVal = aaBox.GetCheck(); 
}

void Cvisualizer101Dlg::OnNMCustomdrawSlider11(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	volumeVal = volumeSlider.GetPos();
	*pResult = 0;
}

void Cvisualizer101Dlg::OnBnClickedCheck4()
{
	peakRotateVal = peakRotateBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck5()
{
	peakFovVal = peakFovBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck6()
{
	animTessVal = animTessBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck2()
{
	lightingVal = lightingBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck3()
{
	explodeVal = explodeBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck7()
{
	wireframeVal = wireframeBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck8()
{
	simpleNormsVal = simpleNormsBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck9()
{
	switchNormsVal = switchNormsBox.GetCheck();
}

void Cvisualizer101Dlg::OnBnClickedCheck10()
{
	texSkyVal = texSkyBox.GetCheck();
}

void Cvisualizer101Dlg::OnCbnSelchangeCombo6()
{
	shapeBox.GetLBText(shapeBox.GetCurSel(), shapreStr);
}

void Cvisualizer101Dlg::OnCbnSelchangeCombo5()
{
	colourBox.GetLBText(colourBox.GetCurSel(), colourStr);
}

void Cvisualizer101Dlg::OnCbnSelchangeCombo4()
{
	// TODO: Add your control notification handler code here
}

void Cvisualizer101Dlg::OnCbnSelchangeCombo3()
{
	// TODO: Add your control notification handler code here
}

void Cvisualizer101Dlg::OnCbnSelchangeCombo7()
{
	deviceBox.GetLBText(deviceBox.GetCurSel(), deviceStr);
}
