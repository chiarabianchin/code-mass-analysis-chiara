#ifndef unfoldData_C
#define unfoldData_C
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TParameter.h>
#include <TString.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <TStopwatch.h>

#include "/data/macros/LoadALICEFigures.C"
#include "/data/Work/code-mass-analysis-chiara/classes/PlotUtilEmcalJetMass.h"
#include "/data/Work/code-mass-analysis-chiara/utils/CommonTools.C"

#include "RooUnfoldBayes.h"

//#include "RooUnfoldTestHarness2D.h"


//const Int_t nPtBins = 5;
//Double_t ptmin[nPtBins] = {20.,40.,60.,80.,100.};
//Double_t ptmax[nPtBins] = {40.,60.,80.,100.,120.};
const Int_t nPtBins = 3;//4;
Double_t ptmin[nPtBins] = {60.,80.,100.}; //40.,
Double_t ptmax[nPtBins] = {80.,100.,120.};//60.,

void Rebin2D(TH2 *hOrig, TH2 *hNew);
Int_t CombineTriggersCorrected(const Int_t ninputs, TString files[], Int_t iteration[], TString legs[]);
void ReturnWiderRangeHistogram(TH2D* h, TH2D* g, TH2D*& hnew, TH2D*& gnew);

void unfold(TString str = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/unfolding/response.root", TString strDat = "../AnalysisResults.root", Int_t iterDef = 2, Int_t iterMin =1, Int_t iterMax = 3, Int_t iList = 0, Int_t varyPrior = 0, Double_t R = 0.4, Int_t colType = 1, TString suff = "", Int_t colorSeries = 1, Int_t isVarW = 0) {
	
	TStopwatch watch;
	watch.Start();
	
	//colType
	// 0 = Pb-Pb
	// 1 = p-Pb MB
	// 2 = p-Pb TR J1
	// 3 = MBbelow 80, EJE above from MassOutput.root
	// 4 = MBbelow 80, EJE+MB above from MassOutput.root
	// if(colType == 3 || colType == 4) iList is used for the background subtraction method: 0 = Const sub; 1 = Raw; 2 = Deriv.
	// Note: the unfolding of the full spectrum gives different results than the separate unfolding of the MB and EJE, due to different hypothesis: - smaller stat unc in the high pt part of the spectrum can bias the unfolding,  - the pt cut at reco level will be different for the separate and full spectrum unfolding, not clear how to make it compatible
	// for colType == 1 and 2 implicitly iList is also pointing to the two bkg sub, but 0 = Deriv, 1 = ConstSub, else Raw
	//colorSeries = 1 (Blue tone), 2 (Green tone)
	// isVarW > 0. Response with varible bin width, to be treated differently
	
	Int_t *colorSeriesArray = 0x0;
	if(colorSeries == 1) colorSeriesArray = colorSeriesBlue;
	if(colorSeries == 2) colorSeriesArray = colorSeriesGreen;
	
	TH2D *hMeasured = 0x0;
	PlotUtilEmcalJetMass *util = 0x0;
	TParameter<Double_t> *par = 0x0; // understand how this is used! Not implemented for MassOutput. How to deal with EJE triggers?
	TString cvsuff = TString::Format("%d_Bkg%d_%s", colType, iList, suff.Data());
	
	TCanvas *cMass = new TCanvas(TString::Format("cMassInput%s", cvsuff.Data()), "Input mass", 800, 400);
	cMass->Divide(2,1);
	
	// read from the merged MB+EJE output. This is not to be used, it gives wrong results 
	if(colType == 3 || colType == 4) {
		Printf("** read from the merged MB+EJE output. This is not to be used, it gives wrong results");
		Printf("Input data file type 3, %s", strDat.Data());
		TFile *fin = new TFile(strDat);
		if(!fin->IsOpen()){
			Printf("... not found");
			return;
		}
		TString hname = TString::Format("h2MPtTagged_BkgSub%d_TrgCmb%d", iList, colType-2); //see MBEJETrigComb in plotMassCompareTriggers (plotJetMasspPb.C)
		
		hMeasured = (TH2D*)fin->Get(hname);
		if(!hMeasured){
			Printf("Histograms %s not found, check inputs", hname.Data());
			fin->ls();
		}
		
	} else if(isVarW > 0) { //-------- Read data from the cooked output in DefineRangeUnfolding
		Printf("** Read data from the cooked output in DefineRangeUnfolding");
		TFile *fin = new TFile(strDat);
		if(!fin->IsOpen()){
			Printf("%s not found", strDat.Data());
			return;
		}
		Int_t npossib = 4;
		TString hMeasNames[npossib] = {"hMassPt", "hMassPtSyD", "hMassPtSyU", "hMeasStdRange"};
		if(isVarW>npossib) {
			Printf("isVarW > npossib %d", isVarW);
			return;
		}
		hMeasured = (TH2D*)fin->Get(hMeasNames[isVarW-1]);
		if(!hMeasured){
			Printf("Mass vs Pt histo %s not found in %s", hMeasNames[isVarW-1].Data(), strDat.Data());
			return;
		}
		Printf("Read histogram %s", hMeasured->GetName() );
	} else {
		//read data directly from the task Jet mass output
		Printf("** read data directly from the task Jet mass output");
		TString trkStr = "PicoTracks_pT0150"; //"tracks_pT0150" // "MCParticlesSelected_pT0000"
		TString cltStr = ""; //"CaloClustersCorr"
		TString jetTStr = "Charged"; //Full
		//gROOT->LoadMacro("$gitJetMass/classes/PlotUtilEmcalJetMass.cxx+");
		util = new PlotUtilEmcalJetMass();
		util->SetTrackString(trkStr);
		util->SetClusterString(cltStr);
		util->SetJTypeString(jetTStr);
		util->SetInputFileName(strDat);
		util->SetJetRadius(R);
		util->SetJetType("Charged");
		util->SetCentBin(0);
		util->LoadFile();
		//Deriv
		if(colType==1)      util->SetTag("_TCDeriv");
		else if(colType==2) util->SetTag("_TCDerivJ1");
		else                util->SetTag("Raw");
		util->SetConstTag("");
		util->LoadList();
		//Const
		if(colType==1)      util->SetTag("");
		else if(colType==2) util->SetTag("ConstJ1");
		else                util->SetTag("Raw");
		util->SetConstTag("ConstSub_TC");
		util->LoadList();
		
		//Raw
		if(colType==1)      util->SetTag("_TCRaw");
		if(colType==2)      util->SetTag("_TCRawJ1");
		//if(colType==1)      util->SetTag("_TC");
		//if(colType==2)      util->SetTag("_TCJ1");
		util->SetConstTag("");
		
		util->LoadList();
		
		hMeasured = util->GetJetMassVsPt(PlotUtilEmcalJetMass::kJetTagged,iList);
		//  TH2D *hMeasured = util->GetJetMassVsPt(PlotUtilEmcalJetMass::kJetAll,iList);
		
		
		Double_t nEvt = 1.;
		if(colType==0) nEvt = util->GetNEvents(0);
		else if(colType==1) nEvt = util->GetNEventsAll(0);
		par = new TParameter<Double_t>("nEvt",nEvt);
	}
	if(!hMeasured) {
		Printf("Problem! measured histogram not found!");
		return;
	}
	cMass->cd(1);
	hMeasured->Draw("colz");
	
	TFile *f = new TFile(str.Data());
	if(!f->IsOpen()){
		Printf("Problem! File with response %s not found, cannot unfold", str.Data());
		return;
	}
	RooUnfoldResponse *resp;
	resp = (RooUnfoldResponse*)f->Get("resp");
	TH2F *hSmear = (TH2F*)f->Get("fh2RespDimM");
	TH2F *hTrue = (TH2F*)f->Get("fh2RespDimT");
	TH2F *hPrior = (TH2F*)f->Get("fh2RespDimP");
	
	if(!hPrior) {
		hPrior = (TH2F*)f->Get("fh2RespDimT"); //use truth as prior
		Printf("!!!! Using fh2RespDimT as prior... do you really want that?");
	}
	hSmear->Sumw2();
	hTrue->Sumw2();
	hPrior->Sumw2();
	
	if(hSmear) Printf("Found hSmear (fh2RespDimM)");
	TH2 *hMeas = 0x0;
	if(isVarW == 0) {
		hMeas = dynamic_cast<TH2*>(hSmear->Clone("hMeas"));
		hMeas->Reset();
		Rebin2D(hMeasured,hMeas);
	}
	else {
		hMeas = dynamic_cast<TH2*>(hMeasured->Clone("hMeas"));
	}
	hMeas->SetTitle("Measured; #it{p}_{T}; #it{M}");
	cMass->cd(2);
	hMeas->Draw("colz");
	Printf("GetN bins hMeasured M = %d, Pt = %d", hMeasured->GetNbinsY(), hMeasured->GetNbinsX());
	Printf("GetN bins hMeas M = %d, Pt = %d", hMeas->GetNbinsY(), hMeas->GetNbinsX());
	Printf("GetN bins hSmear M = %d, Pt = %d", hSmear->GetNbinsY(), hSmear->GetNbinsX());
	
	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;//kNoError;//
	if(iterDef-iterMin < 0){
		Printf("%%%%%%%%%%%% Warning! requesting iteration %d when minimum is %d! Changing def to minimum", iterDef, iterMin);
		iterDef = iterMin;
	}
	
	// iterMin and iterMax are included (see for loop), so need one more space, this is the +1
	//Int_t nIterTmp = iterMax - iterMin +1;
	Int_t nIterTmp = iterMax - iterMin;
	const Int_t nIter = nIterTmp;
	TH2D *hReco[nIter];
	TH2D *hFolded[nIter];
	TMatrixD covmat[nIter];
	RooUnfoldBayes unfold[nIter];
	//for(Int_t iter = iterMin; iter<=iterMax; iter++) {
	for(Int_t iter = iterMin; iter<iterMax; iter++) {
		// RooUnfoldBayes    unfold (resp, hMeas, iter);
		unfold[iter-iterMin] = RooUnfoldBayes(resp, hMeas, iter);
		
		hReco[iter-iterMin] = (TH2D*)unfold[iter-iterMin].Hreco(errorTreatment);
		hReco[iter-iterMin]->SetName(TString::Format("hReco_Iter%d",iter));
		hReco[iter-iterMin]->SetTitle(TString::Format("Result Iter%d; #it{p}_{T};#it{M}",iter));
		// unfold[iter-iterMin].Print();
		hFolded[iter-iterMin] = (TH2D*)resp->ApplyToTruth(hReco[iter-iterMin]);
		hFolded[iter-iterMin]->SetName(TString::Format("hFolded_Iter%d",iter));
		hFolded[iter-iterMin]->SetTitle(TString::Format("Folded Iter%d; #it{p}_{T};#it{M}",iter));
		Printf("Before Ereco");
		covmat[iter-iterMin] = unfold[iter-iterMin].Ereco(RooUnfold::kCovariance);
		Printf("After Ereco");
	}
	TH2D *hPriorFolded = (TH2D*)resp->ApplyToTruth(hPrior);
	hPriorFolded->SetName("hPriorFolded");
	
	TCanvas *c2 = new TCanvas(TString::Format("cRespSummary%s",cvsuff.Data()),"Response and Result summary",800,750);
	c2->Divide(2,2);
	TLegend *legSummary = new TLegend(0.1, 0.55, 0.3, 0.8);
	legSummary->SetFillStyle(0);
	//legSummary->AddEntry();
	c2->cd(1);
	hReco[iterDef-iterMin]->Draw("colz");
	c2->cd(2);
	hMeas->Draw("colz");
	
	//hMeas = raw input data distribution with the same bigging as the response
	//hReco = Result of the unfolding
	TH1D *hRecoP[2];
	hRecoP[0] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin]->ProjectionX("hRecoP_x")); //pt
	hRecoP[1] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin]->ProjectionY("hRecoP_y")); //M
	
	
	TH1D *hFoldedP[2];
	hFoldedP[0] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin]->ProjectionX("hFoldedP_x"));
	hFoldedP[1] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin]->ProjectionY("hFoldedP_y"));
	
	TH1D *hSP[2];
	hSP[0] = dynamic_cast<TH1D*>(hMeas->ProjectionX("hSP_x"));
	hSP[1] = dynamic_cast<TH1D*>(hMeas->ProjectionY("hSP_y"));
	
	TH1D *hTP[2];
	hTP[0] = dynamic_cast<TH1D*>(hTrue->ProjectionX("hTP_x"));
	hTP[1] = dynamic_cast<TH1D*>(hTrue->ProjectionY("hTP_y"));
	
	//make-up
	for(Int_t ii = 0; ii<2; ii++){
		hRecoP[ii]->SetMarkerStyle(20);
		hRecoP[ii]->SetMarkerColor(4);
		hRecoP[ii]->SetLineColor(4);
		hRecoP[ii]->SetLineWidth(3);
		
		hFoldedP[ii]->SetMarkerStyle(24);
		hFoldedP[ii]->SetMarkerColor(kMagenta+2);
		hFoldedP[ii]->SetLineColor(kMagenta+2);
		
		hSP[ii]->SetLineColor(2);
		hSP[ii]->SetMarkerColor(2);
		
		hTP[ii]->SetLineColor(1);
		hTP[ii]->SetMarkerColor(1);
		
	}
	
	TLegend *leg2 = new TLegend(0.1, 0.55, 0.3, 0.8);
	//plotting x-axis projection
	c2->cd(3);
	gPad->SetLogy();
	hSP[0]->Draw();
	hTP[0]->Draw("sames");
	hRecoP[0]->Draw("sames");
	hFoldedP[0]->Draw("Psames");
	
	leg2->AddEntry(hTP[0], "True", "PL");
	leg2->AddEntry(hSP[0], "Raw", "PL");
	leg2->AddEntry(hRecoP[0], "Result ", "P");
	leg2->AddEntry(hFoldedP[0], "Folded ", "P");
	
	//plotting y-axis projection
	c2->cd(4);
	gPad->SetLogy();
	
	hSP[1]->Draw();
	hTP[1]->Draw("sames");
	
	hRecoP[1]->Draw("sames");
	
	hFoldedP[1]->Draw("sames");
	
	c2->cd(1);
	leg2->Draw();
	
	
	
	TH1D *hMUnf[nPtBins][nIter]; // projection of hReco on mass in each pT bin and iteration
	TH1D *hPtUnf[nIter];
	TH1D *hMFol[nPtBins][nIter]; // projection of hFolded on mass in each pT bin and iteration
	TH1D *hMTru[nPtBins];        // projection of hTrue on mass in each pT bin
	TH1D *hMSme[nPtBins];        // projection of hMeas on mass in each pT bin
	TH1D *hMPri[nPtBins];        // projection of hPrior on mass in each pT bin
	TH1D *hMUnfRatiotoTru[nPtBins][nIter];
	TH1D *hMFolRatiotoMeas[nPtBins][nIter];
	
	//uncertainties
	TH1D *hMUncUnf[nPtBins][nIter]; // projection of hReco on mass in each pT bin and iteration
	TH1D *hMRelUncUnf[nPtBins][nIter];
	TH1D *hMeasUnc[nPtBins];
	TH1D *hMeasRelUnc[nPtBins];
	
	TLegend *legIter = new TLegend(0.6, 0.3, 0.9, 0.8, "Iter");
	legIter->SetFillStyle(0);
	legIter->SetBorderSize(0);
	
	for(Int_t i = 0; i<nPtBins; i++) {
		Int_t min = hTrue->GetXaxis()->FindBin(ptmin[i]+0.00001);
		Int_t max = hTrue->GetXaxis()->FindBin(ptmax[i]-0.00001);
		hMTru[i] = dynamic_cast<TH1D*>(hTrue->ProjectionY(TString::Format("hMTru_%d",i),min,max));
		
		min = hMeas->GetXaxis()->FindBin(ptmin[i]+0.00001);
		max = hMeas->GetXaxis()->FindBin(ptmax[i]-0.00001);
		hMSme[i] = dynamic_cast<TH1D*>(hMeas->ProjectionY(TString::Format("hMSme_%d",i),min,max));
		
		min = hPrior->GetXaxis()->FindBin(ptmin[i]+0.00001);
		max = hPrior->GetXaxis()->FindBin(ptmax[i]-0.00001);
		hMPri[i] = dynamic_cast<TH1D*>(hPrior->ProjectionY(TString::Format("hMPri_%d",i),min,max));
		
		//make up
		hMTru[i]->SetLineColor(hTP[0]->GetLineColor());
		
		hMSme[i]->SetLineColor(hSP[0]->GetLineColor());
		hMSme[i]->SetMarkerColor(hSP[0]->GetLineColor());
		
		hMPri[i]->SetLineColor(hTP[0]->GetLineColor());
		hMPri[i]->SetMarkerColor(hTP[0]->GetLineColor());
		hMPri[i]->SetMarkerStyle(3);
		hMPri[i]->SetLineWidth(2);
		
		hMeasUnc[i] = (TH1D*)hMSme[i]->Clone(TString::Format("hMeasUnc_%d", i));
		hMeasUnc[i]->Reset();
		hMeasUnc[i]->SetMarkerStyle(28);
		hMeasUnc[i]->SetMarkerColor(kBlue+2);
		hMeasUnc[i]->SetLineColor(kBlue+2);
		
		hMeasRelUnc[i] = (TH1D*)hMSme[i]->Clone(TString::Format("hMeasRelUnc_%d", i));
		hMeasRelUnc[i]->Reset();
		hMeasRelUnc[i]->SetMarkerStyle(34);
		hMeasRelUnc[i]->SetMarkerColor(kBlue+2);
		hMeasRelUnc[i]->SetLineColor(kBlue+2);
		
		for(Int_t ib = 0; ib<hMeasUnc[i]->GetNbinsX(); ib++){
			Double_t err = hMSme[i]->GetBinError(ib+1), counts = hMSme[i]->GetBinContent(ib+1);
			hMeasUnc[i]->SetBinContent(ib+1, err);
			if(TMath::Abs(counts) > 1e-5) hMeasRelUnc[i]->SetBinContent(ib+1,err/counts);
		}
		
		//for(Int_t iter = iterMin; iter<=iterMax; iter++) {
		for(Int_t iter = iterMin; iter<iterMax; iter++) {
			
			min = hReco[iter-iterMin]->GetXaxis()->FindBin(ptmin[i]+0.00001);
			max = hReco[iter-iterMin]->GetXaxis()->FindBin(ptmax[i]-0.00001);
			hMUnf[i][iter-iterMin] = dynamic_cast<TH1D*>(hReco[iter-iterMin]->ProjectionY(TString::Format("hMUnf_%d_Iter%d",i,iter),min,max));
			Printf("Marker size %f", hMUnf[i][iter-iterMin]->GetMarkerSize());
			
			hMUncUnf[i][iter-iterMin] = (TH1D*)hMUnf[i][iter-iterMin]->Clone(TString::Format("hMUncUnf_%d_Iter%d",i,iter));
			hMUncUnf[i][iter-iterMin]->SetMarkerStyle(28);
			hMUncUnf[i][iter-iterMin]->SetMarkerColor(kOrange+7);
			hMUncUnf[i][iter-iterMin]->Reset();
			hMRelUncUnf[i][iter-iterMin] = (TH1D*)hMUnf[i][iter-iterMin]->Clone(TString::Format("hMRelUncUnf_%d_Iter%d",i,iter));
			hMRelUncUnf[i][iter-iterMin]->SetMarkerStyle(34);
			hMRelUncUnf[i][iter-iterMin]->SetMarkerColor(kOrange+7);
			hMRelUncUnf[i][iter-iterMin]->Reset();
			
			for(Int_t ib = 0; ib<hMUnf[i][iter-iterMin]->GetNbinsX(); ib++){
				Double_t err = hMUnf[i][iter-iterMin]->GetBinError(ib+1), counts = hMUnf[i][iter-iterMin]->GetBinContent(ib+1);
				hMUncUnf[i][iter-iterMin]->SetBinContent(ib+1, err);
				if(TMath::Abs(counts) > 1e-5) hMRelUncUnf[i][iter-iterMin]->SetBinContent(ib+1, err/counts);
			}
			
			min = hFolded[iter-iterMin]->GetXaxis()->FindBin(ptmin[i]+0.00001);
			max = hFolded[iter-iterMin]->GetXaxis()->FindBin(ptmax[i]-0.00001);
			hMFol[i][iter-iterMin] = dynamic_cast<TH1D*>(hFolded[iter-iterMin]->ProjectionY(TString::Format("hMFol_%d_Iter%d", i, iter),min,max));
			
			
			//make up
			hMUnf[i][iter-iterMin]->SetMarkerColor(colorSeriesArray[iter-iterMin]);//hRecoP[0]->GetMarkerColor());
			hMUnf[i][iter-iterMin]->SetMarkerStyle(hRecoP[0]->GetMarkerStyle());
			hMUnf[i][iter-iterMin]->SetLineColor(colorSeriesArray[iter-iterMin]);//hRecoP[0]->GetLineColor());
			hMUnf[i][iter-iterMin]->SetLineWidth(hRecoP[0]->GetLineWidth());
			
			hMFol[i][iter-iterMin]->SetMarkerColor(colorSeriesArray[iter-iterMin]);//hFoldedP[0]->GetMarkerColor());
			hMFol[i][iter-iterMin]->SetMarkerStyle(hFoldedP[0]->GetMarkerStyle());
			hMFol[i][iter-iterMin]->SetLineColor(colorSeriesArray[iter-iterMin]);//hFoldedP[0]->GetLineColor());
			hMFol[i][iter-iterMin]->SetLineWidth(hFoldedP[0]->GetLineWidth());
			
			if(i == 0) {
				legIter->AddEntry(hMUnf[i][iter-iterMin], TString::Format("%d", iter), "PL");
			}
			
			
			hReco[iter-iterMin]->GetXaxis()->SetRange(0, -1);
			hReco[iter-iterMin]->GetYaxis()->SetRange(0, -1);
			hFolded[iter-iterMin]->GetXaxis()->SetRange(0, -1);
			hFolded[iter-iterMin]->GetXaxis()->SetRange(0, -1);
		}
		
	}
	
	TFile *fout = new TFile(TString::Format("UnfoldedDistributionsPrior%d%s.root",varyPrior, suff.Data()), "RECREATE");
	hPrior->Write("fh2Prior");
	hTrue->Write("fh2True");
	hMeas->Write("fh2Smear");
	
	if(util) par->Write();
	resp->Write();
	hPriorFolded->Write();
	//for(Int_t iter = iterMin; iter<=iterMax; iter++) {
	for(Int_t iter = iterMin; iter<iterMax; iter++) {
		hReco[iter-iterMin]->Write();
		hFolded[iter-iterMin]->Write();
		covmat[iter-iterMin].Write(TString::Format("covmat%d",iter));
		unfold[iter-iterMin].Write(TString::Format("unfold%d",iter));
		
		// adding this, is it correct?
		hPtUnf[iter-iterMin] = dynamic_cast<TH1D*>(hReco[iter-iterMin]->ProjectionX(TString::Format("hPtUnf_Iter%d",iter)));
		//Double_t Njets = hPtUnf[i][iter-iterMin]->Integral();
		//Double_t bWdt = hMUnf[i][iter-iterMin]->GetBinWidth(1);
		
		//hMUnf[i][iter-iterMin]->Scale(1./Njets, "width");
		hPtUnf[iter-iterMin]->Write();
	}
	
	for(Int_t i = 0; i<nPtBins; i++) {
		hMPri[i]->Write();
		hMSme[i]->Write();
		hMTru[i]->Write();
		hMeasUnc[i]->Write();
		hMeasRelUnc[i]->Write();
		//for(Int_t iter = iterMin; iter<=iterMax; iter++) {
		for(Int_t iter = iterMin; iter<iterMax; iter++) {
			if(hMUnf[i][iter-iterMin]) hMUnf[i][iter-iterMin]->Write();
			if(hMUncUnf[i][iter-iterMin]) hMUncUnf[i][iter-iterMin]->Write();
			if(hMRelUncUnf[i][iter-iterMin]) hMRelUncUnf[i][iter-iterMin]->Write();
			if(hMFol[i][iter-iterMin]) hMFol[i][iter-iterMin]->Write();
		}
	}
	fout->Write();
	fout->Close();
	
	// break;
	
	TCanvas *c3 = new TCanvas(TString::Format("c3%s",cvsuff.Data()),"c3",800,750);
	c3->Divide(3,2);
	TLegend *leg3 = new TLegend(0.2, 0.4, 0.8, 0.8);
	TCanvas *cIter = new TCanvas(TString::Format("cIter%s",cvsuff.Data()), "Iterations", 800, 750);
	cIter->Divide(3, 2);
	TCanvas *cItertoTrue = new TCanvas(TString::Format("cItertoTrue%s",cvsuff.Data()), "Iterations Ratio to True", 800, 750);
	cItertoTrue->Divide(3, 2);
	TCanvas *cIterFoltoMeas = new TCanvas(TString::Format("cIterFoltoMeas%s",cvsuff.Data()), "Iterations Ratio Folded to Measured", 800, 750);
	cIterFoltoMeas->Divide(3, 2);
	TLegend *legUnc = new TLegend(0.2, 0.4, 0.8, 0.8);
	TLegend *legRUnc = new TLegend(0.2, 0.4, 0.8, 0.8);
	
	TCanvas *cUnc = new TCanvas(TString::Format("cUnc%s",cvsuff.Data()), "Uncertainties", 800, 750);
	cUnc->Divide(3, 2);
	
	TCanvas *cRelUnc = new TCanvas(TString::Format("cRelUnc%s",cvsuff.Data()), "Relative Uncertainties", 800, 750);
	cRelUnc->Divide(3, 2);
	
	for(Int_t i = 0; i<nPtBins; i++) {
		TPaveText *pvpT = new TPaveText(0.2, 0.8, 0.8, 0.9, "NDC");
		pvpT->SetBorderSize(0);
		pvpT->SetFillStyle(0);
		pvpT->AddText(TString::Format("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", ptmin[i], ptmax[i]));
		
		c3->cd(i+1);
		gPad->SetLogy();
		
		hMSme[i]->Draw();
		hMTru[i]->Draw("sames");
		hMPri[i]->Draw("Psames");
		hMUnf[i][iterDef-iterMin]->Draw("Psames");
		hMFol[i][iterDef-iterMin]->Draw("Psames");
		pvpT->Draw();
		if(i==0){
			leg3->AddEntry(hMTru[i], "True", "LP");
			leg3->AddEntry(hMUnf[i][iterDef-iterMin], "Unfolded", "LP");
			leg3->AddEntry(hMFol[i][iterDef-iterMin], "Folded", "L");
			leg3->AddEntry(hMSme[i], "Raw", "LP");
			leg3->AddEntry(hMPri[i], "Priors", "LP");
			
		}
		
		if(i == 0 ){
			legUnc->AddEntry(hMeasUnc[i], "#epsilon Meas", "P");
			legUnc->AddEntry(hMUncUnf[i][iterDef-iterMin], "#epsilon Unfol", "P");
			legRUnc->AddEntry(hMeasRelUnc[i], "#epsilon_{Rel} Meas", "P");
			legRUnc->AddEntry(hMRelUncUnf[i][iterDef-iterMin], "#epsilon_{Rel} Unfol", "P");
		}
		cUnc->cd(i+1);
		hMeasUnc[i]->Draw("P");
		hMUncUnf[i][iterDef-iterMin]->Draw("Psames");
		pvpT->Draw();
		
		cRelUnc->cd(i+1);
		hMeasRelUnc[i]->Draw("P");
		hMRelUncUnf[i][iterDef-iterMin]->Draw("Psames");
		
		pvpT->Draw();
		
		Bool_t firstdraw[3] = {kTRUE, kTRUE, kTRUE};
		for(Int_t j = iterMin; j < iterMax; j++){
			cIter->cd(i+1);
			//hMUnf[i][j]->SetLineColor(colors[j]);
			if(hMUnf[i][j-iterMin]->Integral() > 0){
				if(firstdraw[0]) {
					hMUnf[i][j-iterMin]->GetYaxis()->SetRangeUser(0., 1.5 * hMUnf[i][j-iterMin]->GetMaximum());
					hMUnf[i][j-iterMin]->Draw();
					firstdraw[0] = kFALSE;
				}
				else hMUnf[i][j-iterMin]->Draw("sames");
				legIter->Draw();
				
				
				hMUnfRatiotoTru[i][j-iterMin] = (TH1D*)hMUnf[i][j-iterMin]->Clone(TString::Format("hMUnf_It%dRatiotoTru_Pt%d", j, i));
				hMUnfRatiotoTru[i][j-iterMin]->SetTitle(TString::Format("Unfolded/True; #it{M} (GeV); Unfolded_{Iter}/True"));
				hMUnfRatiotoTru[i][j-iterMin] ->Divide(hMTru[i]);
				
				cItertoTrue->cd(i+1);
				if(firstdraw[1]) {
					hMUnfRatiotoTru[i][j-iterMin]->GetYaxis()->SetRangeUser(0., 2);
					hMUnfRatiotoTru[i][j-iterMin]->Draw();
					firstdraw[1] = kFALSE;
				} else hMUnfRatiotoTru[i][j-iterMin]->Draw("sames");
				legIter->Draw();
				
				//hMFol[i][j]->SetLineColor(colors[j]);
				hMFolRatiotoMeas[i][j-iterMin] = (TH1D*)hMFol[i][j-iterMin]->Clone(TString::Format("hMFol_It%dRatiotoMeas_Pt%d", j, i));
				hMFolRatiotoMeas[i][j-iterMin]->SetTitle(TString::Format("Folded/Measured; #it{M} (GeV); Folded_{Iter}/Measured"));
				hMFolRatiotoMeas[i][j-iterMin]->Divide(hMSme[i]);
				cIterFoltoMeas->cd(i+1);
				if(firstdraw[2]) {
					hMFolRatiotoMeas[i][j-iterMin]->GetYaxis()->SetRangeUser(0., 2);
					hMFolRatiotoMeas[i][j-iterMin]->Draw();
					firstdraw[2] = kFALSE;
				}
				else hMFolRatiotoMeas[i][j-iterMin]->Draw("sames");
				legIter->Draw();
			}
		}
	}
	c3->cd(nPtBins+1);
	leg3->Draw();
	cUnc->cd(nPtBins+1);
	legUnc->Draw();
	cRelUnc->cd(nPtBins+1);
	legRUnc->Draw();
	
	SaveCv(c2);
	SaveCv(c3);
	SaveCv(cUnc);
	SaveCv(cRelUnc);
	SaveCv(cIter);
	SaveCv(cItertoTrue);
	SaveCv(cIterFoltoMeas);
	
	
	delete util;
	
	watch.Stop();
	watch.Print();
}

void Rebin2D(TH2 *hOrig, TH2 *hNew) {
	
	for(Int_t i =1; i<=hOrig->GetXaxis()->GetNbins(); i++) {
		Double_t x = hOrig->GetXaxis()->GetBinCenter(i);
		for(Int_t j =1; j<=hOrig->GetYaxis()->GetNbins(); j++) {
			Double_t y = hOrig->GetYaxis()->GetBinCenter(j);
			Double_t err2 = hOrig->GetBinError(i,j)*hOrig->GetBinError(i,j);
			Int_t k = hNew->GetXaxis()->FindBin(x);
			Int_t l = hNew->GetYaxis()->FindBin(y);
			if(hNew->GetBinContent(k,l)>0.) err2 += hNew->GetBinError(k,l)*hNew->GetBinError(k,l);
			hNew->SetBinContent(k,l,hNew->GetBinContent(k,l)+hOrig->GetBinContent(i,j));
			if(err2>0.) hNew->SetBinError(k,l,TMath::Sqrt(err2));
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------
void RunCombineTriggerCorrected(){
	const Int_t ninputs = 2;
	TString files[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgDeriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/EJE/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgDeriv.root"};
	Int_t iterations[ninputs] = {3,3};
	TString leg[ninputs] = {"MB", "EJE"};
	
	CombineTriggersCorrected(ninputs, files, iterations, leg);
	
}
//-----------------------------------------------------------------------------------------------------------------
Int_t CombineTriggersCorrected(const Int_t ninputs, TString files[], Int_t iteration[], TString legs[]){
	
	/// combine the two (ninputs) triggers. Some manual input is needed concerning the bin regions
	/// Wrong -> consider the 2D distributions, sum them up and project
	/// ATTENTION!!! Need to cut the histogram and do not use overlapping regions, but compare them (to do!!!)
	/// returns the bin number where the shift between MB and EJE is performed (it's the first bin of the EJE). The bin array used to calculate it is Double_t ptmin[nPtBins] = {40.,60.,80.,100.};
	
	TH1::AddDirectory(kFALSE);
	
	// 1D distribution copied
	TString basename1d = "hMUnf_", iters = "_Iter";
	TString basenamept = "hPtUnf"; //_%d_Iter%d",i,iter)
	TH1D *hMUnf[ninputs][nPtBins];
	TH1D *hPtUnf[ninputs];
	// 2D distribution to be elaborated
	TString basename = "hReco_Iter";
	TH2D *hUnfol2D[ninputs];
	TH2D *hUnfol2DRB[ninputs];
	TH2D *hUnfol2DSum = 0x0;
	TCanvas *c2D = new TCanvas("c2D", "2D histograms read from file", 800, 600);
	c2D->Divide(ninputs, 1);
	
	TCanvas *c2DN = new TCanvas("c2DN", "2D histograms new range", 800, 600);
	c2DN->Divide(ninputs, 1);
	
	TCanvas *c2DSum = new TCanvas("c2DSum", "2D histograms sum", 600, 600);
	
	TCanvas *cpT = new TCanvas("cpT", "Pt distributions", 600, 600);
	
	TLegend *legpt = new TLegend(0.5, 0.6, 0.8, 0.8, "#it{p}_{T} distribution");
	legpt->SetFillStyle(0);
	legpt->SetBorderSize(0);
	
	Double_t thresh = 79;
	
	for(Int_t ifi = 0 ; ifi< ninputs; ifi++){
		TFile *fin = new TFile(files[ifi]);
		if(!fin->IsOpen()){
			Printf("File %s not found, continue or exit ", files[ifi].Data());
			continue;
		}
		hUnfol2D[ifi] = 0x0;
		hUnfol2D[ifi] = (TH2D*)fin->Get(TString::Format("%s%d", basename.Data(), iteration[ifi]));
		if(!hUnfol2D[ifi]){
			Printf("The 2D histogram %s%d is not there, check!", basename.Data(), iteration[ifi]);
			fin->ls();
			return -1;
		}
		hUnfol2D[ifi]->SetName(TString::Format("%s_%d", hUnfol2D[ifi]->GetName(), ifi));
		
		c2D->cd(1+ifi);
		hUnfol2D[ifi]->Draw("colz");
		
		hPtUnf[ifi] = (TH1D*)fin->Get(TString::Format("%s%s%d", basenamept.Data(), iters.Data(), iteration[ifi]));
		
		
		if(!hPtUnf[ifi]){
			Printf("Histogram %s%s%d not found", basenamept.Data(), iters.Data(), iteration[ifi]);
			continue;
		}
		hPtUnf[ifi]->SetLineColor(colors[ifi]);
		hPtUnf[ifi]->SetMarkerColor(colors[ifi]);
		hPtUnf[ifi]->SetMarkerStyle(24);
		
		legpt->AddEntry(hPtUnf[ifi], legs[ifi], "Pl");
		
		for(Int_t ipt = 0; ipt < nPtBins; ipt++){
			hMUnf[ifi][ipt] = (TH1D*)fin->Get(TString::Format("%s%d%s%d", basename1d.Data(), ipt, iters.Data(), iteration[ifi]));
			hMUnf[ifi][ipt]->SetLineColor(colors[ifi]);
			hMUnf[ifi][ipt]->SetMarkerColor(colors[ifi]);
			// decide wether it's better to change name
			//hMUnf[ifi][ipt]->SetName(TString::Format("%s%s%dPt%d", basename1d.Data(), iters.Data(), iteration[ifi], ipt));
			
			Int_t binrange[2] = {hPtUnf[ifi]->GetXaxis()->FindBin(ptmin[ipt]), hPtUnf[ifi]->GetXaxis()->FindBin(ptmax[ipt] -0.0001)};
			
			Double_t njets = hPtUnf[ifi]->Integral(binrange[0], binrange[1]);
			Printf("Integral from file %d, pT %d = %f", ifi, ipt, njets);
			if(njets > 0) hMUnf[ifi][ipt]->Scale(1./njets, "width");
			
		}
	}
	
	ReturnWiderRangeHistogram(hUnfol2D[0], hUnfol2D[1], hUnfol2DRB[0], hUnfol2DRB[1]);
	c2DN->cd(1);
	hUnfol2DRB[0]->Draw("colz");
	c2DN->cd(2);
	hUnfol2DRB[1]->Draw("colz");
	
	Int_t nminbin = hUnfol2DRB[0]->GetXaxis()->FindBin(thresh);
	Int_t nmaxbinspt = hUnfol2DRB[0]->GetXaxis()->GetNbins();
	Int_t nmaxbinsm = hUnfol2DRB[0]->GetYaxis()->GetNbins();
	// sum histogram starts from the min bias ...
	hUnfol2DSum = (TH2D*)hUnfol2DRB[0]->Clone("hUnfold2DSum");
	
	Printf("Bin threshold %d + 1 = %d", nminbin, nminbin+1);
	for(Int_t ib = nminbin; ib<=nmaxbinspt; ib++){ // .... above threshold the triggered sample substitute the min bias entries
		for(Int_t ibm = 1; ibm<=nmaxbinsm; ibm++){
			hUnfol2DSum->SetBinContent(ib,ibm, hUnfol2DRB[1]->GetBinContent(ib,ibm));
			hUnfol2DSum->SetBinError(ib,ibm, hUnfol2DRB[1]->GetBinError(ib,ibm));
		}
	}
	c2DSum->cd();
	hUnfol2DSum->Draw("colz");
	
	TCanvas *cMass = new TCanvas("cMass", "Mass Projections", 800, 800);
	cMass->Divide(3,2);
	
	TFile *fout = new TFile(TString::Format("MassUnfSum%s%s.root", legs[0].Data(), legs[1].Data()), "recreate");
	fout->cd();
	hUnfol2DRB[0]->Write();
	hUnfol2DRB[1]->Write();
	hUnfol2DSum->Write();
	
	Int_t ptbinSwitch = FindBinInArray(thresh, ptmin, nPtBins) + 1; //3; // check this !!!!! 
	if(ptbinSwitch<0) {
		Printf("ERROR, out of bound");
		return -1;
	} else Printf("Bin threshold %d", ptbinSwitch);
	
	TH1D* hPtProj = hUnfol2DSum->ProjectionX(TString::Format("hUnfPtpj_Itr%d", iteration[0]));
	hPtProj->SetMarkerStyle(20);
	hPtProj->SetMarkerColor(colors[ninputs]);
	legpt->AddEntry(hPtProj, "2Dglued", "Pl");
	
	cpT->cd();
	gPad->SetLogy();
	hPtProj->Draw("P");
	legpt->Draw();
	//Projections and save to file
	for(Int_t ipt = 0; ipt < nPtBins; ipt++){
		
		Int_t ptbinrange[2] = {hUnfol2DSum->GetXaxis()->FindBin(ptmin[ipt]), hUnfol2DSum->GetXaxis()->FindBin(ptmax[ipt] - 0.001)}; //pt axis
		
		TPaveText *pvt = new TPaveText(0.3, 0.8, 0.8, 0.9, "NDC");
		pvt->SetFillStyle(0);
		pvt->SetBorderSize(0);
		pvt->AddText(TString::Format("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", ptmin[ipt], ptmax[ipt]));
		
		TH1D *hMptbin = hUnfol2DSum->ProjectionY(TString::Format("hUnfMpj_Itr%d_ptb%d", iteration[0], ipt), ptbinrange[0], ptbinrange[1]);
		hUnfol2DSum->SetTitle("Mass projection of the 2D unfolded sum");
		hMptbin->SetLineWidth(2);
		hMptbin->SetLineColor(colors[ninputs]);
		hMptbin->SetMarkerColor(colors[ninputs]);
		hMptbin->SetMarkerStyle(21);
		
		
		Double_t njets = hPtProj->Integral(ptbinrange[0], ptbinrange[1]);
		
		Printf("Integral from 2D, pT %d = %f", ipt, njets);
		
		if(njets > 0) hMptbin->Scale(1./njets, "width");
		
		cMass->cd(ipt+1);
		gPad->SetLogy();
		hMptbin->Draw();
		
		pvt->Draw();
		fout->cd();
		hMptbin ->Write();
		if(!hMUnf[0][ipt] || !hMUnf[1][ipt]){
			Printf("Can't save originals, not found");
			continue;	 
		} else {
			if(ipt<ptbinSwitch){
				if(!legs[0].Contains("MB")) Printf("Check what you're saving @@@@@@@@");
				hMUnf[0][ipt] -> Write();
				cMass->cd(ipt+1);
				hMUnf[0][ipt] ->Draw("sames");
				cpT->cd();
				
				hPtUnf[0]->Draw("sames");
			} else {
				hMUnf[1][ipt] -> Write();
				cMass->cd(ipt+1);
				hMUnf[1][ipt] ->Draw("sames");	
				cpT->cd();
				
				hPtUnf[1]->Draw("sames");
			}
		}
		
	}
	
	SaveCv(c2D);
	SaveCv(c2DN);
	SaveCv(c2DSum);
	SaveCv(cMass);
	SaveCv(cpT);
	
	fout->Close();
	
	return ptbinSwitch;
}

//-----------------------------------------------------------------------------------------------------------------

void ReturnWiderRangeHistogram(TH2D* h, TH2D* g, TH2D*& hnew, TH2D*& gnew){
	
	// define new histograms with range equal to the wider of the two. Bin widht must be the same in the two input histograms
	
	if(!h || !g) {
		Printf("Input null, return");
		return;
	}
	// number of bins -> can be different
	Int_t nbinsh[2] = {h->GetNbinsX(), h->GetNbinsY()};
	Int_t nbinsg[2] = {g->GetNbinsX(), g->GetNbinsY()};
	
	//ranges, can be different, take the wider range
	Double_t rangeLh[2] = {h->GetXaxis()->GetBinLowEdge(1), h->GetYaxis()->GetBinLowEdge(1)}; 
	Double_t rangeHh[2] = {h->GetXaxis()->GetBinLowEdge(nbinsh[0]+1), h->GetYaxis()->GetBinLowEdge(nbinsh[1]+1)};
	Double_t rangeLg[2] = {g->GetXaxis()->GetBinLowEdge(1), g->GetYaxis()->GetBinLowEdge(1)}; 
	Double_t rangeHg[2] = {g->GetXaxis()->GetBinLowEdge(nbinsg[0]+1), g->GetYaxis()->GetBinLowEdge(nbinsg[1]+1)};
	Double_t higherR[2] = {rangeHh[0], rangeHh[1]}; //set to one of the two
	Double_t lowerR[2] = {rangeLh[0], rangeLh[1]};
	
	// bin width: they're hopefully the same, otherwise it's more messy... not implemented
	Double_t binWh[2] = {h->GetXaxis()->GetBinWidth(2), h->GetYaxis()->GetBinWidth(2)};
	Double_t binWg[2] = {g->GetXaxis()->GetBinWidth(2), g->GetYaxis()->GetBinWidth(2)};
	
	Printf("Original histograms: %s: Nbins = %d, %d, Ranges (%.2f, %.2f), (%.2f, %.2f)", h->GetName(), nbinsh[0], nbinsh[1],  rangeLh[0], rangeHh[0], rangeLh[1], rangeHh[1]);
	Printf("\t%s: Nbins = %d, %d, Ranges (%.2f, %.2f), (%.2f, %.2f)", g->GetName(), nbinsg[0], nbinsg[1],  rangeLg[0], rangeHg[0], rangeLg[1], rangeHg[1]);
	
	// perform checks here
	for(Int_t i = 0; i < 2; i++){
		if(TMath::Abs(binWh[i] - binWg[i] < 1e-5)) Printf("(%d) Same bin width = %.3f, can continue! :-) ", i, binWh[i]);
		else Printf(":-(");
		
		//if(higherR[i]<rangeHh[i]) higherR[i] = rangeHh[i]; //not needed because set to it
		if(higherR[i]<rangeHg[i]) higherR[i] = rangeHg[i];
		
		//if(lowerR[i]>rangeLh[i]) lowerR[i] = rangeLh[i]; //not needed because set to it
		if(lowerR[i]>rangeLg[i]) lowerR[i] = rangeLg[i];
	}
	
	Int_t nbinsX =  (higherR[0] - lowerR[0])/binWh[0];
	Int_t nbinsY =  (higherR[1] - lowerR[1])/binWh[1];
	Printf("Output histograms: Nbins = %d, %d, Ranges (%.2f, %.2f), (%.2f, %.2f)", nbinsX, nbinsY,  lowerR[0], higherR[0],  lowerR[1], higherR[1]);
	hnew = new TH2D(TString::Format("%s_n", h->GetName()), ";#it{p}_{T} (GeV/#it{c}); #it{M} (GeV)",  nbinsX, lowerR[0], higherR[0],  nbinsY,  lowerR[1], higherR[1]);
	
	gnew = new TH2D(TString::Format("%s_n", g->GetName()), ";#it{p}_{T} (GeV/#it{c}); #it{M} (GeV)", nbinsX, lowerR[0], higherR[0], nbinsY,  lowerR[1], higherR[1]);
	
	for(Int_t i = 0; i < nbinsX; i++){
		for(Int_t j = 0; j < nbinsY; j++){
			
			// filling histogram h
			Double_t bincontenth = 0., binerrorh = 0;
			
			Int_t binh = h->FindBin(hnew->GetXaxis()->GetBinCenter(i+1), 
				hnew->GetYaxis()->GetBinCenter(j+1));
			//Printf("i = %d, j = %d, bin center %f, %f -> old gl bin %d", i+1, j+1, hnew->GetXaxis()->GetBinCenter(i+1), hnew->GetYaxis()->GetBinCenter(j+1), binh);
			Int_t bxh, byh, z;
			h->GetBinXYZ(binh, bxh, byh, z);
			// Printf("H old: Gl Bin  %d (%d, %d); new: Gl Bin  %d (%d, %d)", binh, bxh, byh, hnew->GetBin(i+1, j+1), i+1, j+1);
			if(bxh > 0 && bxh < nbinsh[0] && byh > 0 && byh < nbinsh[1]){
				bincontenth = h->GetBinContent(binh);
				binerrorh = h->GetBinError(binh);
			}
			//Printf("Filling hnew bin %d, %d with %e +- %e", i+1, j+1, bincontenth, binerrorh);
			hnew->SetBinContent(i+1, j+1, bincontenth);
			hnew->SetBinError(i+1, j+1, binerrorh);
			
			// filling histogram g
			Double_t bincontentg = 0., binerrorg = 0;
			Int_t bing = g->FindBin(gnew->GetXaxis()->GetBinCenter(i+1), 
				gnew->GetYaxis()->GetBinCenter(j+1));
			//Printf("i = %d, j = %d, bin center %f, %f -> old gl bin %d", i+1, j+1, gnew->GetXaxis()->GetBinCenter(i+1), gnew->GetYaxis()->GetBinCenter(j+1), bing);
			Int_t bxg, byg;
			g->GetBinXYZ(bing, bxg, byg, z);
			
			//Printf("G old: Gl Bin  %d (%d, %d); new: Gl Bin  %d (%d, %d)", bing, bxg, byg, gnew->GetBin(i+1, j+1), i+1, j+1);
			if(bxg > 0 && bxg < nbinsg[0] && byg > 0 && byg < nbinsg[1]){
				bincontentg = g->GetBinContent(bing);
				binerrorg = g->GetBinError(bing);
			}
			//Printf("Filling gnew bin %d, %d with %e +- %e", i+1, j+1, bincontentg, binerrorg);
			gnew->SetBinContent(i+1, j+1, bincontentg);
			gnew->SetBinError(i+1, j+1, binerrorg);
			
		}
	}
	return;
} 


#endif

