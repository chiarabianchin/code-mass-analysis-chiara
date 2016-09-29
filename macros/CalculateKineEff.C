#ifndef CalculateKineEff_C
#define CalculateKineEff_C
#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>
#include <THnSparse.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TList.h>

THnSparseF* ReadInput(TString file, TString hname);
THnSparseF* ReadInput(TString file, TString listname, TString hname);
THnSparseF* ReadResp(TString inputfile, TString hrespname, TString listname); // this uses the previous two methods and decides which one according to the value of listname
TH1D* ApplyKineEff(TH1D *hMass, TH1D *hKineNum, TH1D *hKineDen);
void CalculateKineEff(Double_t minPtR, Double_t maxPtR, Double_t minMR, Double_t maxMR, TString inputfile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/DetFluc/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root", TString hrespname = "fhResponseFinal", TString listname = "");


//implementation

THnSparseF* ReadResp(TString inputfile, TString hrespname, TString listname){
	THnSparseF* hRespW = 0x0;
	if(listname.IsNull()) hRespW = ReadInput(inputfile, hrespname);
	else hRespW = ReadInput(inputfile, listname, hrespname);
	
	if(!hRespW) {
		Printf("Response not found");
		return 0 ;
	}
	return hRespW;
}

//_____________________________________________________________________________

void CalculateKineEff(Double_t minPtR, Double_t maxPtR, Double_t minMR, Double_t maxMR, TString inputfile, TString hrespname, TString listname){
	
	THnSparseF* hRespW = ReadResp(inputfile, hrespname, listname);
	
	Int_t nbinsPtP = 4;
	Double_t binsPtP[nbinsPtP+1] = {40, 60, 80, 100., 120};
	Double_t massRangePar[2] = {0., 40.};
	Printf("Setting max particle level to %.0f - %.0f", massRangePar[0], massRangePar[1]);
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nbinsPtP, nx, ny, dx, dy);
	
	Int_t axMp = 1, axPtp = 3 , axPtd = 2, axMd = 0;
	
	TH2D* hMvsPtParFullRange = hRespW->Projection(axMp, axPtp);
	
	Double_t rangePtReco[2] = {minPtR, maxPtR};
	Double_t rangeMReco[2]  = {minMR, maxMR};
	
	Int_t binPtrange[2] = {hRespW->GetAxis(axPtd)->FindBin(rangePtReco[0]), hRespW->GetAxis(axPtd)->FindBin(rangePtReco[1])};
	hRespW->GetAxis(axPtd)->SetRange(binPtrange[0], binPtrange[1]);
	
	Int_t binMrange[2] = {hRespW->GetAxis(axMd)->FindBin(rangeMReco[0]), hRespW->GetAxis(axMd)->FindBin(rangeMReco[1])};
	hRespW->GetAxis(axMd)->SetRange(binMrange[0], binMrange[1]);
	
	
	TH2D* hMvsPtParRecRange = hRespW->Projection(axMp, axPtp);
	
	hMvsPtParFullRange->GetYaxis()->SetRangeUser(massRangePar[0], massRangePar[1]);
	hMvsPtParRecRange->GetYaxis()->SetRangeUser(massRangePar[0], massRangePar[1]);
	
	TH2D* hRatioMvsPtPar = (TH2D*)hMvsPtParRecRange->Clone("hRatioMvsPtParRecROverFull");
	hRatioMvsPtPar->Divide(hMvsPtParFullRange);
	
	TCanvas *cMvsPtPar = new TCanvas(Form("cMvsPtParRangePt%.0f%0.fM%.0f%0.f", rangePtReco[0], rangePtReco[1], rangeMReco[0], rangeMReco[1]), "Mass vs Pt (particle level)", 990, 500);
	cMvsPtPar->Divide(3,1);
	
	TCanvas *c1dProjPtBins = new TCanvas(Form("c1dProjPtBinsRangePt%.0f%0.fM%.0f%0.f", rangePtReco[0], rangePtReco[1], rangeMReco[0], rangeMReco[1]), "1D Mass projections in pt bins", 900, 800);
	c1dProjPtBins->Divide(nx, ny);
	
	TPaveText *pvFullR = new TPaveText(0.35, 0.8, 0.75, 0.9, "NDC");
	pvFullR->SetFillStyle(0);
	pvFullR->SetBorderSize(0);
	pvFullR->AddText("Full #it{p}_{T,rec}, #it{M} range");
	
	TPaveText *pvRecoR = new TPaveText(0.35, 0.8, 0.75, 0.9, "NDC");
	pvRecoR->SetFillStyle(0);
	pvRecoR->SetBorderSize(0);
	pvRecoR->AddText("Used #it{p}_{T,rec}, #it{M} range");
	
	TPaveText *pvRatio = new TPaveText(0.35, 0.8, 0.75, 0.9, "NDC");
	pvRatio->SetFillStyle(0);
	pvRatio->SetBorderSize(0);
	pvRatio->AddText("Ratio");
	
	cMvsPtPar->cd(1);
	hMvsPtParFullRange->Draw("colz");
	pvFullR->Draw();
	
	cMvsPtPar->cd(2);
	hMvsPtParRecRange->Draw("colz");
	pvRecoR->Draw();
	
	cMvsPtPar->cd(3);
	hRatioMvsPtPar->Draw("colz");
	pvRatio->Draw();
	
	TFile *fout = new TFile(Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", rangePtReco[0], rangePtReco[1], rangeMReco[0], rangeMReco[1]), "recreate"); 
	fout->cd();
	hMvsPtParFullRange->Write();
	hMvsPtParRecRange->Write();
	hRatioMvsPtPar->Write();
	
	
	for(Int_t iptp = 0; iptp<nbinsPtP; iptp++){
		Int_t binr[2] = {hRatioMvsPtPar->GetXaxis()->FindBin(binsPtP[iptp]), hRatioMvsPtPar->GetXaxis()->FindBin(binsPtP[iptp+1]-0.01)};
		
		hRatioMvsPtPar->GetXaxis()->SetRange(binr[0], binr[1]);
		
		
		
		TH1D* hNum = hMvsPtParRecRange->ProjectionY(Form("hNumMPtPar%.0f%0.f", binsPtP[iptp], binsPtP[iptp+1]), binr[0], binr[1]);
		TH1D* hDen = hMvsPtParFullRange->ProjectionY(Form("hDenMPtPar%.0f%0.f", binsPtP[iptp], binsPtP[iptp+1]), binr[0], binr[1]);
		
		TH1D* hRatioM = (TH1D*)hNum->Clone(Form("hRatioMPtPar%.0f%0.f", binsPtP[iptp], binsPtP[iptp+1]));
		
		hRatioM->Divide(hDen);
		
		c1dProjPtBins->cd(iptp+1);
		hRatioM->Draw();
		
		fout->cd();
		hNum->Write();
		hDen->Write();
		hRatioM->Write();
		
	}
	
	hRatioMvsPtPar->GetXaxis()->SetRange(0, -1);
	
	SaveCv(cMvsPtPar);
	SaveCv(c1dProjPtBins);
	fout->Close();
	
}

//_____________________________________________________________________________

// read THnSparse from input file
THnSparseF* ReadInput(TString file, TString hname){
	
	TFile *fin = new TFile(file);
	
	if(!fin->IsOpen()){
		Printf("File %s not found, return 0", file.Data());
		return 0x0;
	}
	
	THnSparseF *hsp = (THnSparseF*)fin->Get(hname);
	return hsp;
}

//_____________________________________________________________________________

// read THnSparse from input file
THnSparseF* ReadInput(TString file, TString listname, TString hname){
	
	TList *list = ReadFile(file, listname);
	if(!list ){
		Printf("List %s nor found", listname.Data());
		return 0x0;
	}
	
	
	THnSparseF *hsp = (THnSparseF*)list->FindObject(hname);
	return hsp;
}

//_____________________________________________________________________________

// main method to correct unfolded mass spectra for kinematic efficiency
// it requires the input file of the unfolded mass and the input file of the kine eff saved in the method CalculateKineEff
void CorrectForKineEff(TString fileMassName = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/updatedRooUnfold/UnfoldedDistributionsPrior0.root", TString fileKineEffName = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/analysis/WeightNormPerBin/DetFluc/KineEff/KineEffPtR20_120MR0_12.root", TString fileoutSuff = "", Double_t minPtR = 20, Double_t maxPtR = 140, Double_t minMR = 0, Double_t maxMR = 12){
	
	TFile *finD = new TFile(fileMassName);
	if(!finD->IsOpen()){
		Printf("File %s not found, return ", fileMassName.Data());
		return ;
	}
	
	TFile *finE = new TFile(fileKineEffName);
	if(!finE->IsOpen()){
		Printf("File %s not found, return ", fileKineEffName.Data());
		return ;
	}
	TString basenameNum = "hNumMPtPar";
	TString basenameDen = "hDenMPtPar";
	TString basenameDat = "hMUnf_", suffDat = "_Iter3";
	
	Int_t nptbins = 4;
	Double_t binsPt[nptbins+1] = {40, 60, 80, 100., 120};
	Int_t nx, ny, dx, dy;
	CalculatePads(nptbins, nx, ny, dx, dy);
	
	TCanvas *cMassKineCorr = new TCanvas("cMassKineCorr", "Kine Corrected mass", 800, 900);
	cMassKineCorr->Divide(nx, ny);
	
	TFile *fout = new TFile(Form("CorrectedUnfoldedMass%s.root", fileoutSuff.Data()), "recreate");
	
	for(Int_t ipt = 0; ipt<nptbins; ipt++){
		TH1D* hMassUnf = (TH1D*)finD->Get(Form("%s%d%s", basenameDat.Data(), ipt, suffDat.Data()));
		if(!hMassUnf) {
			Printf("Mass Unf not found");
			finD->ls();
		}
		
		TH1D *hENum = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameNum.Data(), binsPt[ipt], binsPt[ipt+1]));
		TH1D *hEDen = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameDen.Data(), binsPt[ipt], binsPt[ipt+1]));
		hMassUnf->SetName(Form("hMassUnfCorrPt%.0f%.0f", binsPt[ipt], binsPt[ipt+1]));
		
		
		TH1D* hMassUnfCor = ApplyKineEff(hMassUnf, hENum, hEDen);
		hMassUnfCor->SetMarkerColor(colors[5]);
		hMassUnfCor->SetLineColor(colors[5]);
		
		cMassKineCorr->cd(ipt+1);
		//hMassUnfCor->GetXaxis()->SetRangeUser(minMR, maxMR);
		hMassUnfCor->GetXaxis()->SetRangeUser(0, 20);
		hMassUnfCor->Draw();
		hMassUnf->Draw("sames");
	
		fout->cd();
		hMassUnf->Write();
		hMassUnfCor->Write();
	}
	
	SaveCv(cMassKineCorr);
}

// apply KineEff to the undolded mass (takes care of rebinning and then produced the efficiency)
TH1D* ApplyKineEff(TH1D *hMass, TH1D *hKineNum, TH1D *hKineDen){
	if(!hMass) {
		Printf("No valid mass histo");
		return 0x0;
	}
	
	if(!hKineNum) {
		Printf("No valid num histo");
		return 0x0;
	}
	
	if(!hKineDen) {
		Printf("No valid den histo");
		return 0x0;
	}
	
	TH1 *hMass2 = 0x0;
	TH1 *hNumRb = 0x0;
	TH1 *hDenRb = 0x0;
	
	UniformTH1FForDivide(hMass,hKineNum, hMass2, hNumRb);
	delete hMass2;
	
	UniformTH1FForDivide(hMass,hKineDen, hMass2, hDenRb);
	
	Printf("n bins %d and %d (mass %d)" , hNumRb->GetNbinsX(), hDenRb->GetNbinsX(), hMass2->GetNbinsX());
	TH1D *hKineEff = (TH1D*)hNumRb->Clone(Form("%sOver%s", hKineNum->GetName(), hKineDen->GetName()));
	hKineEff->Divide(hDenRb);
	
	TCanvas *cDebug = new TCanvas(Form("c%s", hKineEff->GetName()), Form("c%s", hKineEff->GetName()));
	//hNumRb->Draw();
	//hDenRb->Draw("sames");
	hKineEff->Draw();
	hMass2->Divide(hKineEff);
	
	//
	//Printf("INITIAL INTEGRAL %f, INTEGRAL FINAL %f", hMass->Integral("width"), hMass2->Integral("width"));
	Printf("Info: rescaling after Kine correction");
	hMass2->Scale(1./hMass2->Integral("width"));
	//Printf("SCALED INTEGRAL %f", hMass2->Integral("width"));
	return (TH1D*)hMass2;
}

void StudyKineEffEvolution( Int_t settings, Double_t ptMin = 60., Double_t massMin = 0, TString inputResp = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", TString hrespname = "fhResponseFinal",  TString listname = ""){
	// settings"
	// 0 = vary pT max
	// 1 = vary mass max
	// 2 = vary pT and mass max
	THnSparseF *hResp = ReadResp(inputResp, hrespname, listname);
	if(!hResp) return;
	
	//For MB
	//Double_t ptMin = 10., massMin = 0;
	//For EJE
	
	
	const Int_t nmassMaxbins = 8;
	Double_t massMax[nmassMaxbins] = {10, 14, 18, 20, 24, 30, 40, 80};
	const Int_t nptMaxbins = 6;
	Double_t ptMax[nptMaxbins] = {100, 110, 120, 130, 140, 150};
	Int_t nbinsPtP = 4;
	Double_t binsPtP[nbinsPtP+1] = {40, 60, 80, 100., 120};
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nbinsPtP, nx, ny, dx, dy);
	
	if(settings == 0){
		TCanvas *cEffpTMaxVar = new TCanvas("cEffpTMaxVar", "Efficiencies pT max variation", 900, 800);
		cEffpTMaxVar->Divide(nx, ny);
		
		// mass full range, scan pT
		TLegend *legPtScan = new TLegend(0.3, 0.2, 0.7, 0.7);
		legPtScan->SetFillStyle(0);
		legPtScan->SetBorderSize(0);
		
		for(Int_t ipmax = 0; ipmax < nptMaxbins; ipmax++){
			TString filename = Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", ptMin, ptMax[ipmax], massMin, massMax[nmassMaxbins-1]);
			
			//TFile *finE = new TFile(filename);
			//if(!finE->IsOpen()){
			//	finE->Close();
			//	
			//}
			CalculateKineEff(ptMin, ptMax[ipmax], massMin, massMax[nmassMaxbins-1], inputResp, hrespname, listname);
			
			TFile *finE = new TFile(filename);
			if(!finE->IsOpen()){
				Printf("No kine eff file");
				return;
			}
			
			TString basenameNum = "hNumMPtPar";
			TString basenameDen = "hDenMPtPar";
			TString basenameEff = "hRatioMPtPar";
			
			for(Int_t iptp = 0; iptp<nbinsPtP; iptp++){
				TH1D *hENum = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameNum.Data(), binsPtP[iptp], binsPtP[iptp+1]));
				TH1D *hEDen = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameDen.Data(), binsPtP[iptp], binsPtP[iptp+1]));
				TH1D *hEffi = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameEff.Data(), binsPtP[iptp], binsPtP[iptp+1]));
				
				if(!hEffi){
					Printf("Eff %s%.0f%0.f not found", basenameEff.Data(), binsPtP[iptp], binsPtP[iptp+1]);
					finE->ls();
					continue;
				}
				hEffi->SetName(Form("%s_MaxPt%.0f", finE->GetName(), ptMax[ipmax]));
				hEffi->SetLineColor(colors[ipmax]);
				hEffi->SetLineWidth(2);
				
				cEffpTMaxVar->cd(iptp+1);
				if(ipmax>0) hEffi->Draw("sames");
				else {
					hEffi->GetYaxis()->SetRangeUser(0., 1.1);
					hEffi->Draw();
				}
				if(iptp == 0) legPtScan->AddEntry(hEffi, Form("PtMax = %.0f", ptMax[ipmax]), "L");
			}
		}
		cEffpTMaxVar->cd(1);
		legPtScan->Draw();
		
		SaveCv(cEffpTMaxVar);
	}
	
	if(settings == 1){
		
		TCanvas *cEffMassMaxVar = new TCanvas("cEffMassMaxVar", "Efficiencies M max variation", 900, 800);
		cEffMassMaxVar->Divide(nx, ny);
		
		// mass full range, scan pT
		TLegend *legMScan = new TLegend(0.3, 0.2, 0.7, 0.7);
		legMScan->SetFillStyle(0);
		legMScan->SetBorderSize(0);
		
		for(Int_t ipmax = 0; ipmax < nmassMaxbins; ipmax++){
			TString filename = Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", ptMin, ptMax[nptMaxbins-1], massMin, massMax[ipmax]);
			
			
			CalculateKineEff(ptMin, ptMax[nptMaxbins-1], massMin, massMax[ipmax], inputResp, hrespname, listname);
			
			TFile *finE = new TFile(filename);
			if(!finE->IsOpen()){
				Printf("No kine eff file");
				return;
			}
			
			TString basenameNum = "hNumMPtPar";
			TString basenameDen = "hDenMPtPar";
			TString basenameEff = "hRatioMPtPar";
			
			for(Int_t iptp = 0; iptp<nbinsPtP; iptp++){
				TH1D *hENum = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameNum.Data(), binsPtP[iptp], binsPtP[iptp+1]));
				TH1D *hEDen = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameDen.Data(), binsPtP[iptp], binsPtP[iptp+1]));
				TH1D *hEffi = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameEff.Data(), binsPtP[iptp], binsPtP[iptp+1]));
				
				if(!hEffi){
					Printf("Eff %s%.0f%0.f not found", basenameEff.Data(), binsPtP[iptp], binsPtP[iptp+1]);
					finE->ls();
					continue;
				}
				hEffi->SetName(Form("%s_MaxM%.0f", finE->GetName(), massMax[ipmax]));
				hEffi->SetLineColor(colors[ipmax]);
				hEffi->SetMarkerColor(colors[ipmax]);
				hEffi->SetMarkerStyle(markers[ipmax]);
				hEffi->SetLineWidth(2);
				
				cEffMassMaxVar->cd(iptp+1);
				if(ipmax>0) hEffi->Draw("sames");
				else {
					hEffi->GetYaxis()->SetRangeUser(0., 1.1);
					hEffi->Draw();
				}
				if(iptp == 0) legMScan->AddEntry(hEffi, Form("MassMax = %.0f", massMax[ipmax]), "LP");
			}
		}
		cEffMassMaxVar->cd(1);
		legMScan->Draw();
		
		SaveCv(cEffMassMaxVar);
	}
	
	if(settings == 2){
		
		TCanvas *cEffpTMMaxVar = new TCanvas("cEffpTMMaxVar", "Efficiencies pT and M max variation", 900, 800);
		cEffpTMMaxVar->Divide(nx, ny);
		
		// mass full range, scan pT
		TLegend *legPtMScan = new TLegend(0.3, 0.2, 0.7, 0.7);
		legPtMScan->SetFillStyle(0);
		legPtMScan->SetBorderSize(0);
		
		for(Int_t ipmax = 0; ipmax < nptMaxbins; ipmax++){
			for(Int_t immax = 0; immax < nmassMaxbins; immax++){
				TString filename = Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", ptMin, ptMax[ipmax], massMin, massMax[immax]);
				
				//TFile *finE = new TFile(filename);
				//if(!finE->IsOpen()){
				//	finE->Close();
				//	
				//}
				CalculateKineEff(ptMin, ptMax[ipmax], massMin, massMax[immax], inputResp, hrespname, listname);
				
				TFile *finE = new TFile(filename);
				if(!finE->IsOpen()){
					Printf("No kine eff file");
					return;
				}
				
				TString basenameNum = "hNumMPtPar";
				TString basenameDen = "hDenMPtPar";
				TString basenameEff = "hRatioMPtPar";
				
				for(Int_t iptp = 0; iptp<nbinsPtP; iptp++){
					TH1D *hENum = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameNum.Data(), binsPtP[iptp], binsPtP[iptp+1]));
					TH1D *hEDen = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameDen.Data(), binsPtP[iptp], binsPtP[iptp+1]));
					TH1D *hEffi = (TH1D*)finE->Get(Form("%s%.0f%0.f", basenameEff.Data(), binsPtP[iptp], binsPtP[iptp+1]));
					
					if(!hEffi){
						Printf("Eff %s%.0f%0.f not found", basenameEff.Data(), binsPtP[iptp], binsPtP[iptp+1]);
						finE->ls();
						continue;
					}
					hEffi->SetName(Form("%s_MaxPt%.0f_MaxM%.0f", finE->GetName(), ptMax[ipmax], massMax[immax]));
					hEffi->SetLineColor(colors[ipmax]);
					hEffi->SetMarkerColor(colors[ipmax]);
					hEffi->SetMarkerStyle(markers[immax]);
					hEffi->SetLineWidth(2);
					
					cEffpTMMaxVar->cd(iptp+1);
					if(ipmax == 0 && immax == 0) {
						hEffi->GetYaxis()->SetRangeUser(0., 1.1);
						hEffi->Draw();
					}
					else {
						hEffi->Draw("sames");
					}
					if(iptp == 0) legPtMScan->AddEntry(hEffi, Form("PtMax = %.0f, MMax = %.0f", ptMax[ipmax],massMax[immax]), "L");
				}
			}
			cEffpTMMaxVar->cd(1);
			legPtMScan->Draw();
			
			SaveCv(cEffpTMMaxVar);
		}
	}
}
#endif
