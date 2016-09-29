#include "/data/Work/MyCodeJetMass/utils/CommonTools.C"
#include "/data/macros/LoadALICEFigures.C"
#include "/data/Work/MyCodeJetMass/classes/PlotUtilEmcalJetMass.h"
#include <TLine.h>
#include <TParameter.h>

void CompareProjections(Int_t nh2, TH2D** h2d, const Int_t nbins, Double_t binlims[], TString legt[], Bool_t bscale, const char* namepj = "hMass", Double_t ex = 0.001);

void DefineRangeUnfolding(TH2D* hMeasured, Int_t bkgType = 2, Int_t triggerType = 1, Double_t binWPt = 10., Double_t binWMa = 2., Double_t minCounts = 10, Double_t nwidth = 1);

void DrawSquare(TVirtualPad *pad, Double_t xl, Double_t yd, Double_t xr, Double_t yu, Int_t color);

//main

void DefineRangeUnfolding(TString inputData = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/AnalysisResultsMB.root", Int_t bkgType = 2, Int_t triggerType = 1, Double_t binWPt = 10., Double_t binWMa = 2., Double_t minCounts = 10, Double_t nwidth = 1){
	
	// bkgType:  0 = Deriv; 1 = Const sub; 2 = Raw;
	// triggerType: 1 = MB, 2 = J1
	
	
	TString trkStr = "PicoTracks_pT0150"; //"tracks_pT0150" // "MCParticlesSelected_pT0000"
	TString cltStr = ""; //"CaloClustersCorr"
	TString jetTStr = "Charged"; //Full
	Double_t R = 0.4;
	//gROOT->LoadMacro("$gitJetMass/classes/PlotUtilEmcalJetMass.cxx+");
	PlotUtilEmcalJetMass *util = new PlotUtilEmcalJetMass();
	util->SetTrackString(trkStr);
	util->SetClusterString(cltStr);
	util->SetJTypeString(jetTStr);
	util->SetInputFileName(inputData);
	util->SetJetRadius(R);
	util->SetJetType("Charged");
	util->SetCentBin(0);
	util->LoadFile();
	
	//Deriv
	if(triggerType==1)      util->SetTag("_TCDeriv");
	else if(triggerType==2) util->SetTag("_TCDerivJ1");
	else                util->SetTag("Raw");
	util->SetConstTag("");
	util->LoadList();
	//Const
	if(triggerType==1)      util->SetTag("");
	else if(triggerType==2) util->SetTag("ConstJ1");
	else                util->SetTag("Raw");
	util->SetConstTag("ConstSub_TC");
	util->LoadList();
	
	//Raw
	if(triggerType==1)      util->SetTag("_TCRaw");
	if(triggerType==2)      util->SetTag("_TCRawJ1");
	util->SetConstTag("");
	
	util->LoadList();
	
	TH2D *hMeasured = util->GetJetMassVsPt(PlotUtilEmcalJetMass::kJetTagged, bkgType);
	
	DefineRangeUnfolding(hMeasured, bkgType, triggerType, binWPt, binWMa, minCounts, nwidth);
	
	return;
}

//___________________________________________________________________________

void DefineRangeUnfoldingMCClosure(TString inputData = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", Int_t bkgType = 2, Int_t triggerType = 1, Double_t binWPt = 10., Double_t binWMa = 2., Double_t minCounts = 10, Double_t nwidth = 1){
	
	TString hname = "fhResponseFinal";
	
	TFile *fin = new TFile(inputData);
	if(!fin->IsOpen()){
		Printf("File %s not found", inputData.Data());
		return;
	}
	
	THnSparseF *hsph = (THnSparseF*)fin->Get(hname.Data());
	if(!hsph) {
		Printf("%s not found", hname.Data());
		fin->ls();
		return;
	}
	Int_t ndim = hsph->GetNdimensions();
	
	Int_t axMdet = 0, axPtdet = 2;
	
	TH2D* hReco = (TH2D*)hsph->Projection(axMdet, axPtdet);
	
	DefineRangeUnfolding(hReco, bkgType, triggerType, binWPt, binWMa, minCounts, nwidth);
	
	return;
}

//___________________________________________________________________________


//implementations

void DrawSquare(TVirtualPad *pad, Double_t xl, Double_t yd, Double_t xr, Double_t yu, Int_t color){
	
	//x1, y1, x2, y2
	TLine *lineVr = new TLine(xr, yd, xr, yu);
	lineVr->SetLineColor(color);
	lineVr->SetLineWidth(2);
	
	TLine *lineVl = new TLine(xl, yd, xl, yu);
	lineVl->SetLineColor(color);
	lineVl->SetLineWidth(2);
	
	TLine *lineHu = new TLine(xl, yu, xr, yu);
	lineHu->SetLineColor(color);       
	lineHu->SetLineWidth(2);
	
	TLine *lineHd = new TLine(xl, yd, xr, yd);
	lineHd->SetLineColor(color);
	lineHd->SetLineWidth(2);
	
	pad->cd();
	lineVr->Draw();
	lineVl->Draw();
	lineHu->Draw();
	lineHd->Draw();

}

//___________________________________________________________________________

void DefineRangeUnfolding(TH2D* hMeasured, Int_t bkgType, Int_t triggerType, Double_t binWPt, Double_t binWMa, Double_t minCounts, Double_t nwidth){
	
	TString bkgS = "Deriv";
	if(bkgType == 1) bkgS = "Const";
	if(bkgType == 2) bkgS = "NoSub";
	
	TString trgS = "MB";
	if(triggerType == 2) trgS = "J1";
	
	
	Double_t binWOrigPt = hMeasured->GetXaxis()->GetBinWidth(3);
	Double_t binWOrigMa = hMeasured->GetYaxis()->GetBinWidth(3);
	Int_t rebPt = (Int_t) (binWPt / binWOrigPt);
	Int_t rebMa = (Int_t) (binWMa / binWOrigMa);
	Printf("Rebin factor x axis (%s) = %d, Rebin factor y axis (%s) = %d",  hMeasured->GetXaxis()->GetTitle(), rebPt, hMeasured->GetYaxis()->GetTitle(), rebMa);
	
	Int_t rebinF[2] = {rebPt, rebMa};
	TH2D *hMeasReb = (TH2D*)hMeasured->Rebin2D(rebPt, rebMa, Form("%sRb", hMeasured->GetName()));
	
	Double_t maxMa = 30, maxPt = hMeasReb->GetXaxis()->GetBinLowEdge(hMeasReb->GetXaxis()->GetNbins());
	
	TCanvas *cMassPtNumb = new TCanvas(Form("cMassW%.0fPtW%.0fNumbBkg%sTr%s", binWMa, binWPt, bkgS.Data(), trgS.Data()), "Mass vs Pt content", 500, 500);
	cMassPtNumb->cd();
	hMeasReb->GetYaxis()->SetRangeUser(hMeasReb->GetYaxis()->GetBinLowEdge(1), 30.);
	hMeasReb->Draw("text");	
	
	//variable bin width
	// - MB
	const Int_t nbinsXmb = 16;
	Double_t binlimsXmb[nbinsXmb+1] = {0., 5, 10., 15., 20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 100., 120., 140.};
	const Int_t nbinsYmb = 9;
	Double_t binlimsYmb[nbinsYmb+1] = {0, 1., 2., 4., 6., 8., 10., 12., 14., 20.};
	
	// -EJE
	const Int_t nbinsXje = 12;
	Double_t binlimsXje[nbinsXje+1] = {55., 60, 65., 70., 75., 80., 85., 90., 95., 100., 110., 120., 140.};
	const Int_t nbinsYje = 10;
	Double_t binlimsYje[nbinsYje+1] = {0, 1., 2., 4., 6.,  8., 10., 12., 14., 16., 20};
	
	TH2D* hMeasuredRbVarW = 0x0;
	if(triggerType == 1) hMeasuredRbVarW = RebinVariableWidth2D(nbinsXmb, binlimsXmb, nbinsYmb, binlimsYmb, hMeasured);
	
	if(triggerType == 2) hMeasuredRbVarW = RebinVariableWidth2D(nbinsXje, binlimsXje, nbinsYje, binlimsYje, hMeasured);
	
	TCanvas *cMassRebVarW = new TCanvas(Form("cMassRebVarWNumbBkg%sTr%s", bkgS.Data(), trgS.Data()), Form("Mass vs pT with variable bin width"), 500, 500);
	cMassRebVarW->cd();
	hMeasuredRbVarW->Draw("text");
	Double_t ptminT = 10, mminT = 0, ptmaxT = 140, mmaxT = 40;
	//DrawSquare(cMassRebVarW->cd(), ptminT, mminT, ptmaxT, mmaxT, kBlack);
	
	//standard cuts 
	// - MB sample
	Double_t ptupStd = 80., ptdoStd = 20., maupStd = 12, madoStd = 0., ptupStdVw = 80., ptdoStdVw = 20., maupStdVw = 12, madoStdVw = 0.;
	// - EJE sample
	if(triggerType == 2) {
		ptupStd = 120.; 
		ptdoStd = 70.;
		maupStd = 14;
		madoStd = 0.;
		
		ptupStdVw = 120.;
		ptdoStdVw = 65.;
		maupStdVw = 16.;
		madoStdVw = 0.;
	}
	if(bkgType == 0) madoStd = -10; // for derivative
	
	//x1, y1, x2, y2
	TLine *linePtUpStd = new TLine(ptupStd, madoStd, ptupStd, maupStd);
	linePtUpStd->SetLineColor(kMagenta+2);
	linePtUpStd->SetLineWidth(2);
	
	TLine *linePtDoStd = new TLine(ptdoStd, madoStd, ptdoStd, maupStd);
	linePtDoStd->SetLineColor(kMagenta+2);
	linePtDoStd->SetLineWidth(2);
	
	TLine *lineMaUpStd = new TLine(ptdoStd, maupStd, ptupStd, maupStd);
	lineMaUpStd->SetLineColor(kMagenta+2);       
	lineMaUpStd->SetLineWidth(2);
	
	TLine *lineMaDoStd = new TLine(ptdoStd, madoStd, ptupStd, madoStd);
	lineMaDoStd->SetLineColor(kMagenta+2);
	lineMaDoStd->SetLineWidth(2);
	
	cMassPtNumb->cd();
	linePtUpStd->Draw();
	linePtDoStd->Draw();
	lineMaUpStd->Draw();
	lineMaDoStd->Draw();
	
	const Int_t npTbinsStdRange = (Int_t)((ptupStd - ptdoStd)/binWPt);
	const Int_t nMabinsStdRange = (Int_t)((maupStd - madoStd)/binWMa);
	
	Double_t limsptStd[npTbinsStdRange+1];
	for(Int_t ipt = 0; ipt<npTbinsStdRange+1; ipt++) limsptStd[ipt] = ptdoStd + ipt*binWPt;
	Double_t limsmaStd[nMabinsStdRange+1];
	for(Int_t ipt = 0; ipt<nMabinsStdRange+1; ipt++) limsmaStd[ipt] = madoStd + ipt*binWMa;
	TH2D *hMeasStdRange = RebinVariableWidth2D(npTbinsStdRange, limsptStd, nMabinsStdRange, limsmaStd, hMeasured);
	hMeasStdRange->SetName("hMeasStdRange");
	//new TH2D("hMeasStdRange", "Mass vs Pt (standard binning); #it{p}_{T} (GeV/#it{c}; M (GeV/#it{c}^2)", npTbinsStdRange, ptdoStd, ptupStd, nMabinsStdRange, madoStd, maupStd);
	
	TLine *linePtUpStdVw = new TLine(ptupStdVw, madoStdVw, ptupStdVw, maupStdVw);
	linePtUpStdVw->SetLineColor(kMagenta+2);
	linePtUpStdVw->SetLineWidth(2);
	
	TLine *linePtDoStdVw = new TLine(ptdoStdVw, madoStdVw, ptdoStdVw, maupStdVw);
	linePtDoStdVw->SetLineColor(kMagenta+2);
	linePtDoStdVw->SetLineWidth(2);
	
	TLine *lineMaUpStdVw = new TLine(ptdoStdVw, maupStdVw, ptupStdVw, maupStdVw);
	lineMaUpStdVw->SetLineColor(kMagenta+2);       
	lineMaUpStdVw->SetLineWidth(2);
	
	TLine *lineMaDoStdVw = new TLine(ptdoStdVw, madoStdVw, ptupStdVw, madoStdVw);
	lineMaDoStdVw->SetLineColor(kMagenta+2);
	lineMaDoStdVw->SetLineWidth(2);
	
	cMassRebVarW->cd();
	linePtUpStdVw->Draw();
	linePtDoStdVw->Draw();
	lineMaUpStdVw->Draw();
	lineMaDoStdVw->Draw();
	
	//variation UP
	// - MB sample
	Double_t ptupSysU = ptupStd+ nwidth*binWPt, ptdoSysU = ptdoStd+ nwidth*binWPt, maupSysU = maupStd+ binWMa, madoSysU = 0.;
	
	Int_t binSyUPtup = 0;
	Int_t binSyUPtdo = 0;
	Int_t binSyUMaup = 0;
	Int_t binSyUMado = 0;
	
	Double_t ptupSysUVw = 0., ptdoSysUVw = 0., maupSysUVw = 0., madoSysUVw = 0.;
	
	if(triggerType == 1) {
		
		// +1 selects the following limit in the array
		binSyUPtup = FindBinInArray(ptupStdVw, binlimsXmb, nbinsXmb+1) + 1;
		binSyUPtdo = FindBinInArray(ptdoStdVw, binlimsXmb, nbinsXmb+1) + 1;
		binSyUMaup = FindBinInArray(maupStdVw, binlimsYmb, nbinsYmb+1) + 1;
		//binSyUMado = FindBinInArray(madoStdVw, binlimsYmb, nbinsYmb+1) + 1;
		binSyUMado = 0;
		
		ptupSysUVw = binlimsXmb[binSyUPtup];
		ptdoSysUVw = binlimsXmb[binSyUPtdo];
		maupSysUVw = binlimsYmb[binSyUMaup];
		
	}
	if(triggerType == 2) {
		binSyUPtup = FindBinInArray(ptupStdVw, binlimsXje, nbinsXje+1) + 1;
		binSyUPtdo = FindBinInArray(ptdoStdVw, binlimsXje, nbinsXje+1) + 1;
		binSyUMaup = FindBinInArray(maupStdVw, binlimsYje, nbinsYje+1) + 1;
		//binSyUMado = FindBinInArray(madoStdVw, binlimsYje, nbinsYje+1) + 1;
		binSyUMado = 0;
		
		ptupSysUVw = binlimsXje[binSyUPtup];
		ptdoSysUVw = binlimsXje[binSyUPtdo];
		maupSysUVw = binlimsYje[binSyUMaup];
		
	} 
	Printf("binSyUPtdo = %d, ptdoStdVw = %f", binSyUPtdo, ptdoStdVw);
	Printf("binSyUMup = %d, ptdoStdVw = %f", binSyUMaup, maupStdVw);
	// - EJE sample
	//if(triggerType == 2) {
	//	 ptupSysU = 120.; 
	//	 ptdoSysU = 80.;
	//	 maupSysU = 16;
	//	 madoSysU = 0.;
	//}
	if(bkgType == 0) madoSysU = madoStd + binWMa; // for derivative
	
	//x1, y1, x2, y2
	TLine *linePtUpSysU = new TLine(ptupSysU, madoSysU, ptupSysU, maupSysU);
	linePtUpSysU->SetLineColor(kOrange+2);
	linePtUpSysU->SetLineWidth(2);
	
	TLine *linePtDoSysU = new TLine(ptdoSysU, madoSysU, ptdoSysU, maupSysU);
	linePtDoSysU->SetLineColor(kOrange+2);
	linePtDoSysU->SetLineWidth(2);
	
	TLine *lineMaUpSysU = new TLine(ptdoSysU, maupSysU, ptupSysU, maupSysU);
	lineMaUpSysU->SetLineColor(kOrange+2);       
	lineMaUpSysU->SetLineWidth(2);
	
	TLine *lineMaDoSysU = new TLine(ptdoSysU, madoSysU, ptupSysU, madoSysU);
	lineMaDoSysU->SetLineColor(kOrange+2);
	lineMaDoSysU->SetLineWidth(2);
	
	cMassPtNumb->cd();
	//linePtUpSysU->Draw();
	//linePtDoSysU->Draw();
	//lineMaUpSysU->Draw();
	//lineMaDoSysU->Draw();
	
	
	TLine *linePtUpSysUVw = new TLine(ptupSysUVw, madoSysUVw, ptupSysUVw, maupSysUVw);
	linePtUpSysUVw->SetLineColor(kOrange+2);
	linePtUpSysUVw->SetLineWidth(2);
	
	TLine *linePtDoSysUVw = new TLine(ptdoSysUVw, madoSysUVw, ptdoSysUVw, maupSysUVw);
	linePtDoSysUVw->SetLineColor(kOrange+2);
	linePtDoSysUVw->SetLineWidth(2);
	
	TLine *lineMaUpSysUVw = new TLine(ptdoSysUVw, maupSysUVw, ptupSysUVw, maupSysUVw);
	lineMaUpSysUVw->SetLineColor(kOrange+2);       
	lineMaUpSysUVw->SetLineWidth(2);
	
	TLine *lineMaDoSysUVw = new TLine(ptdoSysUVw, madoSysUVw, ptupSysUVw, madoSysUVw);
	lineMaDoSysUVw->SetLineColor(kOrange+2);
	lineMaDoSysUVw->SetLineWidth(2);
	
	cMassRebVarW->cd();
	linePtUpSysUVw->Draw();
	linePtDoSysUVw->Draw();
	lineMaUpSysUVw->Draw();
	lineMaDoSysUVw->Draw();
	
	//variation DOWN
	
	Int_t binSyDPtup = 0;
	Int_t binSyDPtdo = 0;
	Int_t binSyDMaup = 0;
	Int_t binSyDMado = 0;
	
	Double_t ptupSysDVw = 0, ptdoSysDVw = 0, maupSysDVw = 0, madoSysDVw = 0.;
	
	if(triggerType == 1){
		// -1 selects the previous limit in the array
		// the mass limit stays at the default (0)
		binSyDPtup = FindBinInArray(ptupStdVw, binlimsXmb, nbinsXmb+1) - 1;
		binSyDPtdo = FindBinInArray(ptdoStdVw, binlimsXmb, nbinsXmb+1) - 1;
		binSyDMaup = FindBinInArray(maupStdVw, binlimsYmb, nbinsYmb+1) - 1;
		//binSyDMado = FindBinInArray(madoStdVw, binlimsYmb, nbinsYmb+1) - 1;
		binSyDMado = 0;
		
		ptupSysDVw = binlimsXmb[binSyDPtup];
		ptdoSysDVw = binlimsXmb[binSyDPtdo];
		maupSysDVw = binlimsYmb[binSyDMaup];
		
	}
	if(triggerType == 2){
		binSyDPtup = FindBinInArray(ptupStdVw, binlimsXje, nbinsXje+1) - 1;
		binSyDPtdo = FindBinInArray(ptdoStdVw, binlimsXje, nbinsXje+1) - 1;
		binSyDMaup = FindBinInArray(maupStdVw, binlimsYje, nbinsYje+1) - 1;
		//binSyDMado = FindBinInArray(madoStdVw, binlimsYje, nbinsYje+1) - 1;
		binSyDMado = 0;
		
		ptupSysDVw = binlimsXje[binSyDPtup];
		ptdoSysDVw = binlimsXje[binSyDPtdo];
		maupSysDVw = binlimsYje[binSyDMaup];
	}
	
	
	// - MB sample
	Double_t ptupSysD = ptupStd - nwidth*binWPt, ptdoSysD = ptdoStd - nwidth*binWPt, maupSysD = maupStd - binWMa, madoSysD = 0.;
	
	
	
	if(bkgType == 0) madoSysD = madoStd - binWMa; // for derivative
	// - EJE sample
	//if(triggerType == 2) {
	//	 ptupSysD = 100.; 
	//	 ptdoSysD = 60.;
	//	 maupSysD = 12;
	//	 madoSysD = 0.;
	//}
	
	//x1, y1, x2, y2
	TLine *linePtUpSysD = new TLine(ptupSysD, madoSysD, ptupSysD, maupSysD);
	linePtUpSysD->SetLineColor(kGreen+4);
	linePtUpSysD->SetLineWidth(2);
	
	TLine *linePtDoSysD = new TLine(ptdoSysD, madoSysD, ptdoSysD, maupSysD);
	linePtDoSysD->SetLineColor(kGreen+4);
	linePtDoSysD->SetLineWidth(2);
	
	TLine *lineMaUpSysD = new TLine(ptdoSysD, maupSysD, ptupSysD, maupSysD);
	lineMaUpSysD->SetLineColor(kGreen+4);       
	lineMaUpSysD->SetLineWidth(2);
	
	TLine *lineMaDoSysD = new TLine(ptdoSysD, madoSysD, ptupSysD, madoSysD);
	lineMaDoSysD->SetLineColor(kGreen+4);
	lineMaDoSysD->SetLineWidth(2);
	
	cMassPtNumb->cd();
	linePtUpSysD->Draw();
	linePtDoSysD->Draw();
	lineMaUpSysD->Draw();
	lineMaDoSysD->Draw();
	
	TLine *linePtUpSysDVw = new TLine(ptupSysDVw, madoSysDVw, ptupSysDVw, maupSysDVw);
	linePtUpSysDVw->SetLineColor(kGreen+4);
	linePtUpSysDVw->SetLineWidth(2);
	
	TLine *linePtDoSysDVw = new TLine(ptdoSysDVw, madoSysDVw, ptdoSysDVw, maupSysDVw);
	linePtDoSysDVw->SetLineColor(kGreen+4);
	linePtDoSysDVw->SetLineWidth(2);
	
	TLine *lineMaUpSysDVw = new TLine(ptdoSysDVw, maupSysDVw, ptupSysDVw, maupSysDVw);
	lineMaUpSysDVw->SetLineColor(kGreen+4);       
	lineMaUpSysDVw->SetLineWidth(2);
	
	TLine *lineMaDoSysDVw = new TLine(ptdoSysDVw, madoSysDVw, ptupSysDVw, madoSysDVw);
	lineMaDoSysDVw->SetLineColor(kGreen+4);
	lineMaDoSysDVw->SetLineWidth(2);
	
	cMassRebVarW->cd();
	linePtUpSysDVw->Draw();
	linePtDoSysDVw->Draw();
	lineMaUpSysDVw->Draw();
	lineMaDoSysDVw->Draw();
	
	Printf("Standard cuts: %.0f & %.0f & %.0f & %.0f", ptdoStd, ptupStd, madoStd, maupStd);
	Printf("Sys variation: %.0f & %.0f & %.0f & %.0f", ptdoSysU, ptupSysU, madoSysU, maupSysU);
	Printf("Sys variation: %.0f & %.0f & %.0f & %.0f", ptdoSysD, ptupSysD, madoSysD, maupSysD);
	SaveCv(cMassPtNumb);
	
	SaveCv(cMassRebVarW);
	
	// prepare the input for unfolding
	
	Int_t binptupStdVw = FindBinInArray(ptupStdVw, binlimsXmb, nbinsXmb+1);
	Int_t binptdoStdVw = FindBinInArray(ptdoStdVw, binlimsXmb, nbinsXmb+1);
	Int_t binmaupStdVw = FindBinInArray(maupStdVw, binlimsYmb, nbinsYmb+1);
	Int_t binmadoStdVw = FindBinInArray(madoStdVw, binlimsYmb, nbinsYmb+1);
	if(triggerType == 2){
		binptupStdVw = FindBinInArray(ptupStdVw, binlimsXje, nbinsXje+1);
		binptdoStdVw = FindBinInArray(ptdoStdVw, binlimsXje, nbinsXje+1);
		binmaupStdVw = FindBinInArray(maupStdVw, binlimsYje, nbinsYje+1);
		binmadoStdVw = FindBinInArray(madoStdVw, binlimsYje, nbinsYje+1);
	}
	const Int_t nbinsM = binmaupStdVw - binmadoStdVw;
	Double_t limsM[nbinsM+1];
	
	const Int_t nbinsPt = binptupStdVw - binptdoStdVw;
	Double_t limsPt[nbinsPt+1];
	
	const Int_t nbinsSyDM = binSyDMaup - binSyDMado;
	Double_t limsSyDM[nbinsSyDM+1];
	
	const Int_t nbinsSyUM = binSyUMaup - binSyUMado;
	Double_t limsSyUM[nbinsSyUM+1];
	
	const Int_t nbinsSyDPt = binSyDPtup - binSyDPtdo;
	Double_t limsSyDPt[nbinsSyDPt+1];
	const Int_t nbinsSyUPt = binSyUPtup - binSyUPtdo;
	Double_t limsSyUPt[nbinsSyUPt+1];
	
	TFile *foutMassvsPt = new TFile(Form("MassVsPtVarWBkg%sTr%s.root", bkgS.Data(), trgS.Data()), "recreate");
	TParameter<Int_t> nbinsMpar("nbinsM"    , nbinsM);
	TParameter<Int_t> nbinsPpar("nbinsPt"   , nbinsPt);
	TParameter<Int_t> nbinsMSyDpar("nbinsMSyD" , nbinsSyDM);
	TParameter<Int_t> nbinsPSyDpar("nbinsPtSyD", nbinsSyDPt);
	TParameter<Int_t> nbinsMSyUpar("nbinsMSyU" , nbinsSyUM);
	TParameter<Int_t> nbinsPSyUpar("nbinsPtSyU", nbinsSyUPt);
	TParameter<Int_t> nbinsMFixpar("nbinsMFix" , nMabinsStdRange);
	TParameter<Int_t> nbinsPFixpar("nbinsPtFix", npTbinsStdRange);
	
	foutMassvsPt->cd();
	nbinsMpar.Write();
	nbinsPpar.Write();
	nbinsMSyDpar.Write();
	nbinsPSyDpar.Write();
	nbinsMSyUpar.Write();
	nbinsPSyUpar.Write();
	nbinsMFixpar.Write();
	nbinsPFixpar.Write();
	
	Int_t id = 0;
	// mass default
	Printf("Mass default");
	for(Int_t ibm = binmadoStdVw; ibm < binmaupStdVw+1; ibm++){
		if(triggerType == 1) limsM[id] = binlimsYmb[ibm];
		if(triggerType == 2) limsM[id] = binlimsYje[ibm];
		Printf("%d, ibm = %d -> lim %f", ibm, id, limsM[id]);
		
		TParameter<Double_t> mlimpar(Form("limM%d", id), limsM[id]);
		foutMassvsPt->cd();
		mlimpar.Write();
		
		id++;
	}
	id = 0;
	Printf("Mass Systematic down");
	// mass systematic down
	for(Int_t ibm = binSyDMado; ibm < binSyDMaup+1; ibm++){
		if(triggerType == 1) limsSyDM[id] = binlimsYmb[ibm];
		if(triggerType == 2) limsSyDM[id] = binlimsYje[ibm];
		Printf("%d, ibm = %d -> lim %f", ibm, id, limsSyDM[id]);
		
		TParameter<Double_t> mlimpar(Form("limMSyD%d", id), limsSyDM[id]);
		foutMassvsPt->cd();
		mlimpar.Write();
		
		id++;
	}
	id = 0;
	Printf("Mass Systematic up");
	Printf("Here binSyUMado = %d", binSyUMado);
	// mass systematic up
	for(Int_t ibm = binSyUMado; ibm < binSyUMaup+1; ibm++){
		if(triggerType == 1) limsSyUM[id] = binlimsYmb[ibm];
		if(triggerType == 2) limsSyUM[id] = binlimsYje[ibm];
		Printf("%d, ibm = %d -> lim %f", ibm, id, limsSyUM[id]);
		
		TParameter<Double_t> mlimpar(Form("limMSyU%d", id), limsSyUM[id]);
		foutMassvsPt->cd();
		mlimpar.Write();
		
		id++;
	}
	id = 0;
	//pt default
	Printf("pT default");
	for(Int_t ibp = binptdoStdVw; ibp < binptupStdVw+1; ibp++){
		if(triggerType == 1) limsPt[id] = binlimsXmb[ibp];
		if(triggerType == 2) limsPt[id] = binlimsXje[ibp];
		Printf("%d, ibm = %d -> lim %f", ibp, id, limsPt[id]);
		TParameter<Double_t> ptlimpar(Form("limPt%d", id), limsPt[id]);
		foutMassvsPt->cd();
		ptlimpar.Write();
		
		id++;
	}
	id = 0;
	Printf("pT Systematic down");
	// pt systematic down
	for(Int_t ibp = binSyDPtdo; ibp < binSyDPtup+1; ibp++){
		if(triggerType == 1) limsSyDPt[id] = binlimsXmb[ibp];
		if(triggerType == 2) limsSyDPt[id] = binlimsXje[ibp];
		Printf("%d, ibm = %d -> lim %f", ibp, id, limsSyDPt[id]);
		
		TParameter<Double_t> ptlimpar(Form("limPtSyD%d", id), limsSyDPt[id]);
		foutMassvsPt->cd();
		ptlimpar.Write();
		
		id++;
	}
	id = 0;
	Printf("pT Systematic up");
	// pt systematic up
	for(Int_t ibp = binSyUPtdo; ibp < binSyUPtup+1; ibp++){
		if(triggerType == 1) limsSyUPt[id] = binlimsXmb[ibp];
		if(triggerType == 2) limsSyUPt[id] = binlimsXje[ibp];
		Printf("%d, ibm = %d -> lim %f", ibp, id, limsSyUPt[id]);
		TParameter<Double_t> ptlimpar(Form("limPtSyU%d", id), limsSyUPt[id]);
		foutMassvsPt->cd();
		ptlimpar.Write();
		
		id++;
	}
	
	//definition of the 2D histograms...
	TH2D* h2MassPt = RebinVariableWidth2D(nbinsPt, limsPt, nbinsM, limsM, hMeasuredRbVarW);
	h2MassPt->SetName("hMassPt");
	//new TH2D("hMassPt", "Mass vs Pt (default); #it{p}_{T} (GeV/#it{c}); #it{M} (GeV/#it{c}^{2})", nbinsPt, limsPt, nbinsM, limsM);
	Printf("Default range X(%d) = %.0f - %.0f; Y(%d) = %.0f - %.0f", nbinsPt, limsPt[0], limsPt[nbinsPt], nbinsM, limsM[0], limsM[nbinsM]);
	h2MassPt->Sumw2();
	TH2D* h2MassPtSyD = RebinVariableWidth2D(nbinsSyDPt, limsSyDPt, nbinsSyDM, limsSyDM, hMeasuredRbVarW);
	h2MassPtSyD->SetName("hMassPtSyD");
	//new TH2D("hMassPtSyD", "Mass vs Pt (Syst var Down); #it{p}_{T} (GeV/#it{c}); #it{M} (GeV/#it{c}^{2})", nbinsSyDPt, limsSyDPt, nbinsSyDM, limsSyDM);
	h2MassPtSyD->Sumw2();
	
	TH2D* h2MassPtSyU = RebinVariableWidth2D(nbinsSyUPt, limsSyUPt, nbinsSyUM, limsSyUM, hMeasuredRbVarW);
	h2MassPtSyU->SetName("hMassPtSyU");
	//new TH2D("hMassPtSyU", "Mass vs Pt (Syst var Up); #it{p}_{T} (GeV/#it{c}); #it{M} (GeV/#it{c}^{2})", nbinsSyUPt, limsSyUPt, nbinsSyUM, limsSyUM);
	h2MassPtSyU->Sumw2();
	
	//... and filling
	
	//for(Int_t ibg = 0; ibg < hMeasuredRbVarW->GetNbinsX() * hMeasuredRbVarW->GetNbinsY(); ibg++){
	//	Int_t bx, by, bz;
	//	hMeasuredRbVarW->GetBinXYZ(ibg+1, bx, by, bz);
	//	Double_t content = hMeasuredRbVarW->GetBinContent(ibg+1);
	//	Double_t error   = hMeasuredRbVarW->GetBinError(ibg+1);
	//	
	//	Double_t x, y, z;
	//	x = hMeasuredRbVarW->GetXaxis()->GetBinCenter(bx);
	//	y = hMeasuredRbVarW->GetYaxis()->GetBinCenter(by);
	//	z = hMeasuredRbVarW->GetZaxis()->GetBinCenter(bz);
	//	
	//	
	//	Int_t bglStd = h2MassPt->FindBin(x, y, z);
	//	//Printf("binglob = %d, x = %.0f, y = %.0f, z = %.0f, fill with %f", ibg, x, y, z, content);
	//	h2MassPt->GetBinXYZ(bglStd, bx, by, bz);
	//	if(bx>0 && by > 0 && bx<nbinsPt+1 && by < nbinsM+1){
	//	Printf("std filled bin %d/%d (x = %.0f, y = %.0f, z = %.0f) with %e", bglStd, (nbinsPt+1)*(nbinsM+1), x, y, z, content);
	//	h2MassPt->SetBinContent(bglStd, content);
	//	h2MassPt->SetBinError(bglStd, error);
	//	}
	//	
	//	Int_t bglSyU = h2MassPtSyU->FindBin(x, y, z);
	//	h2MassPt->GetBinXYZ(bglSyU, bx, by, bz);
	//	if(bx>0 && by > 0 && bx<nbinsSyUPt+1 && by < nbinsSyUM+1){
	//		//Printf("syU filled bin %d", bglSyU);
	//		h2MassPtSyU->SetBinContent(bglSyU, content);
	//		h2MassPtSyU->SetBinError(bglSyU, error);
	//	}
	//	
	//	Int_t bglSyD = h2MassPtSyD->FindBin(x, y, z);
	//	h2MassPt->GetBinXYZ(bglSyD, bx, by, bz);
	//	if(bx>0 && by > 0 && bx<nbinsSyDPt+1 && by < nbinsSyDM+1){
	//		//Printf("syD filled bin %d", bglSyD);
	//		h2MassPtSyD->SetBinContent(bglSyD, content);
	//		h2MassPtSyD->SetBinError(bglSyD, error);
	//	}
	//	
	//	
	//}
	
	//histogram with fixed binning
	Printf("N bins and Range of the fixed bin width histogram X %d, %.0f-%.0f - Y  %d, %.0f-%.0f", npTbinsStdRange, ptdoStd, ptupStd, nMabinsStdRange, madoStd, maupStd);
	
	TParameter<Double_t> ptminStd = TParameter<Double_t>("limPtMin", ptdoStd);
	TParameter<Double_t> ptmaxStd = TParameter<Double_t>("limPtMax", ptupStd);
	TParameter<Double_t> maminStd = TParameter<Double_t>("limMMin" , madoStd);
	TParameter<Double_t> mamaxStd = TParameter<Double_t>("limMMax" , maupStd);
	foutMassvsPt->cd();
	ptminStd.Write();
	ptmaxStd.Write();
	maminStd.Write();
	mamaxStd.Write();
	
	//for(Int_t ibg = 0; ibg < hMeasReb->GetNbinsX()+1 * hMeasReb->GetNbinsY()+1; ibg++){
	//	Int_t bx, by, bz;
	//	hMeasReb->GetBinXYZ(ibg+1, bx, by, bz);
	//	Double_t content = hMeasReb->GetBinContent(ibg+1);
	//	Double_t error   = hMeasReb->GetBinError(ibg+1);
	//	
	//	Double_t x, y, z;
	//	x = hMeasReb->GetXaxis()->GetBinCenter(bx);
	//	y = hMeasReb->GetYaxis()->GetBinCenter(by);
	//	z = hMeasReb->GetZaxis()->GetBinCenter(bz);
	//	
	//	Int_t bglStd = hMeasStdRange->FindBin(x, y, z);
	//	hMeasStdRange->GetBinXYZ(bglStd, bx, by, bz);
	//	if(bx>0 && by > 0 && bx<npTbinsStdRange+1 && by<nMabinsStdRange+1){
	//		//Printf("std filled bin %d/%d (x = %.0f, y = %.0f, z = %.0f)", bglStd, npTbinsStdRange*nMabinsStdRange, x, y, z);
	//		hMeasStdRange->SetBinContent(bglStd, content);
	//		hMeasStdRange->SetBinError(bglStd, error);
	//	}
	//}
	
	foutMassvsPt->cd();
	h2MassPt->Write();
	h2MassPtSyD->Write();
	h2MassPtSyU->Write();
	hMeasStdRange->Write();
	
	foutMassvsPt->Close();
}

//___________________________________________________________________________

void CompareProjections(Int_t nh2, TH2D** h2d, const Int_t nbins, Double_t binlims[], TString legt[], Bool_t bscale, const char* namepj, Double_t ex){
	
	Int_t nx, ny, dx, dy;
	
	if(nbins>3) CalculatePads(nbins, nx, ny, dx, dy, 2);
	else CalculatePads(nbins, nx, ny, dx, dy);
	
	TCanvas *ccomp = new TCanvas(Form("cNbins%dNh%d%s%s", nbins, nh2, namepj, bscale ? "scaled" : ""), "Comparisons", dx, dy);
	ccomp->Divide(nx, ny);
	
	TLegend *leg = new TLegend(0.1, 0.5, 0.3, 0.8);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	
	for(Int_t ih = 0; ih<nh2; ih++){ //histogram index
		if(!h2d[ih]) continue;
		
		for(Int_t ib = 0; ib<nbins; ib++){ //slice in one variable (X) and project on the other (Y)
			
			
			Int_t range[2] = {h2d[ih]->GetXaxis()->FindBin(binlims[ib]), h2d[ih]->GetXaxis()->FindBin(binlims[ib+1] - ex)};
			
			TH1D* h1d = h2d[ih]->ProjectionY(Form("h%s%d%s_%.0f_%.0f", namepj, ih, bscale ? "scaled" : "", binlims[ib], binlims[ib+1]), range[0], range[1]);
			if(bscale) h1d->SetTitle(Form("%s_scaled", h1d->GetTitle()));
			
			h1d->SetLineColor(colors[ih]);
			h1d->SetMarkerColor(colors[ih]);
			h1d->SetMarkerStyle(20+4*ih);
			
			if(bscale) h1d->Scale(1./h1d->Integral("width"));
			
			if(ib == 0){
				leg->AddEntry(h1d, legt[ih], "PL");
			}
			ccomp->cd(ib+1);
			if(ih == 0) h1d->Draw();
			else  {
				TPaveText *pvlims = new TPaveText(0.3, 0.8, 0.8, 0.9, "NDC");
				pvlims->SetBorderSize(0);
				pvlims->SetFillStyle(0);
				pvlims->AddText(Form("%.0f < %s < %.0f", binlims[ib], h2d[ih]->GetXaxis()->GetTitle(), binlims[ib+1]));
				h1d->Draw("sames");
				pvlims->Draw();
			}
		}
		
	}
	leg->Draw();
	SaveCv(ccomp);
	
}

//___________________________________________________________________________

void DrawRawMassDifferentBinning(Bool_t scaled = kFALSE, TString inputfile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrJ1.root"){
	
	const Int_t nh = 2;
	TString hnames[nh] = {"hMassPt", "hMeasStdRange"};
	TString legtex[nh] = {"Variable", "Fixed"};
	TFile *fin = new TFile(inputfile);
	if(!fin->IsOpen()){
		Printf("File %s not found", inputfile.Data());
		return;
	}
	
	TH2D** h2list = new TH2D*[nh];
	
	for(Int_t ih = 0; ih < nh; ih++){
		h2list[ih] = (TH2D*)fin->Get(hnames[ih]);
		if(!h2list[ih]) Printf("Warning, %s not found", hnames[ih].Data());
	}
	
	//const Int_t nbins = 3;
	//Double_t binlims[nbins+1] = {60, 80, 100, 120};
	const Int_t nbins = 6;
	Double_t binlims[nbins+1] = {60, 70, 80, 90, 100, 110, 120};
	CompareProjections(nh, h2list, nbins, binlims, legtex, scaled);
}
