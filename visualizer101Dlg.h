
// visualizer101Dlg.h : header file
//

#pragma once
#include "afxwin.h"
#include "afxcmn.h"


// Cvisualizer101Dlg dialog
class Cvisualizer101Dlg : public CDialogEx
{
// Construction
public:
	Cvisualizer101Dlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_visualizer101_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
//	int channel;
	afx_msg void OnBnClickedButton1();
	afx_msg void OnCbnSelchangeCombo1();
	CComboBox channelBox;
	CString channelStr;
	CEdit widthBox;
	CString widthStr;
	CEdit heightBox;
	CString heightStr;
	CComboBox multisamplingBox;
	CString multisamplingStr;
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnEnChangeEdit2();
	afx_msg void OnCbnSelchangeCombo2();
	CSliderCtrl orderSlider;
	int orderVal;
	afx_msg void OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult);
	CStatic orderLabel;
	CStatic uincLabel;
	CStatic xincLabel;
	CSliderCtrl uincSlider;
	CSliderCtrl xincSlider;
	int uincVal;
	afx_msg void OnNMCustomdrawSlider2(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSlider6(NMHDR *pNMHDR, LRESULT *pResult);
	int xincVal;
	CComboBox sampleRateBox;
	CComboBox numChannelsBox;
	CStatic skySplitsLabel;
	CStatic srevSplitsLabel;
	CStatic tiLabel;
	CStatic toLabel;
	CSliderCtrl skySplitsSlider;
	CSliderCtrl srevSplitsSlider;
	CSliderCtrl tiSlider;
	CSliderCtrl toSlider;
	int skySplitsVal;
	int srevSplitsVal;
	int tiVal;
	int toVal;
	afx_msg void OnNMCustomdrawSlider7(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSlider8(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSlider9(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSlider10(NMHDR *pNMHDR, LRESULT *pResult);
	CButton aaBox;
	BOOL aaVal;
	afx_msg void OnBnClickedCheck1();
	CButton peakRotateBox;
	BOOL peakRotateVal;
	CButton peakFovBox;
	BOOL peakFovVal;
	CButton animTessBox;
	BOOL animTessVal;
	CSliderCtrl volumeSlider;
	int volumeVal;
	afx_msg void OnNMCustomdrawSlider11(NMHDR *pNMHDR, LRESULT *pResult);
	CButton lightingBox;
	BOOL lightingVal;
	CButton explodeBox;
	BOOL explodeVal;
	CButton wireframeBox;
	BOOL wireframeVal;
	CButton simpleNormsBox;
	BOOL simpleNormsVal;
	CButton switchNormsBox;
	BOOL switchNormsVal;
	CButton texSkyBox;
	BOOL texSkyVal;
	afx_msg void OnBnClickedCheck4();
	afx_msg void OnBnClickedCheck5();
	afx_msg void OnBnClickedCheck6();
	afx_msg void OnBnClickedCheck2();
	afx_msg void OnBnClickedCheck3();
	afx_msg void OnBnClickedCheck7();
	afx_msg void OnBnClickedCheck8();
	afx_msg void OnBnClickedCheck9();
	afx_msg void OnBnClickedCheck10();
	CComboBox colourBox;
	CComboBox shapeBox;
	CString colourStr;
	CString shapreStr;
	afx_msg void OnCbnSelchangeCombo6();
	afx_msg void OnCbnSelchangeCombo5();
	afx_msg void OnCbnSelchangeCombo4();
	afx_msg void OnCbnSelchangeCombo3();
	CComboBox deviceBox;
	CString deviceStr;
	afx_msg void OnCbnSelchangeCombo7();
};
