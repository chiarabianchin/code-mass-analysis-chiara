#ifndef PrepareRhoWeightingFactors_C
#define PrepareRhoWeightingFactors_C

#include <TGrid.h>
#include <THnSparse.h>

#include "/data/Work/MyCodeJetMass/utils/CommonTools.C"
#include "/data/macros/LoadALICEFigures.C"

void PrepareRhoWeightingFactors(TString fPathRhoDistr, TString fListRhoDistr = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv", TString fNameTHnSparseRhoDistr = "fhnRhoVsRhoMVsLeadJetPtVsMassVsCent", TString suff = "Deriv"){

	const Int_t    fnPtDetBinsForRho = 2;
	Double_t       fPtDetBinsForRho[fnPtDetBinsForRho+1] = {0, 40., 120.};
	TH1D           **fhRhoData  = new TH1D*[fnPtDetBinsForRho];
	TH1D           **fhRhoMData = new TH1D*[fnPtDetBinsForRho];
	TH1D           **fhRhoDataR  = new TH1D*[fnPtDetBinsForRho];
	TH1D           **fhRhoMDataR = new TH1D*[fnPtDetBinsForRho];
	
	Int_t          binNorm = 0; // ASSUMPTION: the bin 0-40 is the same as the response, so we use it as normalization factor
	
	
	TCanvas *cRho = new TCanvas(Form("cRho%s", suff.Data()), "Rho distributions", 600, 600);
	TCanvas *cRhoR = new TCanvas(Form("cRhoR%s", suff.Data()), "Rho ratios", 600, 600);
	
	TFile* fout = new TFile(Form("RhoFromDataRatiosNPtDetbins%d%s.root", fnPtDetBinsForRho, suff.Data()), "recreate");
	
	//prepare rho from data to be used as weighting factor
	if(!fPathRhoDistr.IsNull() && !fListRhoDistr.IsNull() && !fNameTHnSparseRhoDistr.IsNull()){
		if(fPathRhoDistr.Contains("alien")) {
			TGrid::Connect("alien://");
		}
		TFile *f = TFile::Open(fPathRhoDistr);
		if(f->IsOpen()){
			TList *listin = (TList*)f->Get(fListRhoDistr);
			
			if(listin){
				
				THnSparseF* hspRho = (THnSparseF*)listin->FindObject(fNameTHnSparseRhoDistr);
				if(hspRho) {
					// project the thnsparse into the rho and rho m axes in the ranges specified
					Int_t axrho = 0, axrhom = 1;
					Int_t axptdet = 2;
					
					for(Int_t ib = 0; ib < fnPtDetBinsForRho; ib++){
						Int_t rangebins[2] = {hspRho->GetAxis(axptdet)->FindBin(fPtDetBinsForRho[ib]), hspRho->GetAxis(axptdet)->FindBin(fPtDetBinsForRho[ib+1]-0.001)};
						
						hspRho->GetAxis(axptdet)->SetRange(rangebins[0], rangebins[1]);
						fhRhoData[ib] = hspRho->Projection(axrho);
						fhRhoData[ib]->SetName(Form("hRhoData_Bin%d", ib));
						fhRhoData[ib]->Scale(1./fhRhoData[ib]->Integral());
						fhRhoMData[ib] = hspRho->Projection(axrhom);
						fhRhoMData[ib]->SetName(Form("hRhoMData_Bin%d", ib));
						
						fhRhoData[ib]->SetLineColor(colors[ib]);
						fhRhoData[ib]->SetLineWidth(2);
						fhRhoData[ib]->SetMarkerColor(colors[ib]);
						cRho->cd();
						if(ib == 0) fhRhoData[ib]->Draw();
						else fhRhoData[ib]->Draw("sames");
						
						fout->cd();
						fhRhoData[ib]->Write();
					}
					Printf("Created %d rho distributions from data", fnPtDetBinsForRho);
					
				} else Printf("THnSparse %s not found", fNameTHnSparseRhoDistr.Data());
			} else Printf("List %s not found, cannot read the rho distribution", fListRhoDistr.Data());
		} else Printf("File %s not found, cannot read the rho distribution", fPathRhoDistr.Data());
			
		for(Int_t ib = 0; ib < fnPtDetBinsForRho; ib++){
			
			fhRhoDataR[ib] = (TH1D*)fhRhoData[ib]->Clone(Form("hRhoDataR_Bin%d", ib));
			fhRhoDataR[ib]->Divide(fhRhoData[binNorm]);
			cRhoR->cd();
			if(ib == 0){ 
				fhRhoDataR[ib]->GetYaxis()->SetRangeUser(1.e-5, 10.);
				fhRhoDataR[ib]->Draw();
			}
			else fhRhoDataR[ib]->Draw("sames");
			fhRhoDataR[ib]->Write();
			
		}
	}
	
	
}

//------------------------------------------------------------------------
void CheckEmbeddingOutput(TString embfile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/RhoWeighted/AnalysisResults.root", TString emblist1 = "JetShapeConst_JetRhosub_AKTChargedR040_PicoTracksPtH", TString emblist2 = "_pT0150_E_scheme_TC", TString inputRho = "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/RhoFromDataRatiosNPtDetbins2Const.root"){
	
	TString hrespname = "fhnMassResponse_0";
	TString hrhobname = "hRhoData_Bin";
	Int_t   nbins     = 2;
	Int_t   npthb     = 10;
	
	TCanvas *cCmp = new TCanvas("cCmp", "Comparison of rho distributions", 600, 600);
	TLegend *leg = new TLegend(0.5, 0.4, 0.9, 0.8);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	for(Int_t i = 0; i < npthb ; i++ ){
		TString emblist = Form("%s%d%s", emblist1.Data(), i, emblist2.Data());
		TList *list = ReadFile(embfile, emblist);
		if(!list) continue;
		
		THnSparseF *hembResp = (THnSparseF*)list->FindObject(hrespname);
		
		if(!hembResp) {
			Printf("Response not found");
			continue;
		}
		Int_t minptbin = hembResp->GetAxis(3)->FindBin(40.); //pt det
		hembResp->GetAxis(3)->SetRange(minptbin, hembResp->GetAxis(3)->GetNbins());
		TH1D *hRhopjEmbResp = hembResp->Projection(6); //project on rho
		hRhopjEmbResp->SetName(Form("hRhopjEmbResp_PtH%d", i));
		hRhopjEmbResp->SetLineWidth(2);
		hRhopjEmbResp->SetLineColor(colors[i+nbins]);
		if(hRhopjEmbResp->Integral() < 1e-4) {
			hRhopjEmbResp = 0x0;
			continue;
		} hRhopjEmbResp->Scale(1./hRhopjEmbResp->Integral());
		leg->AddEntry(hRhopjEmbResp, Form("Emb pThard %d", i), "L");
		cCmp->cd();
		if(i == 0) hRhopjEmbResp->Draw();
		else hRhopjEmbResp->Draw("sames");
	}
	
	TFile *fin = new TFile(inputRho);
	if(!fin->IsOpen()){
		Printf("Data rho not found");
		return;
	}
	
	
	
	for(Int_t i = 0; i<nbins; i++){
		TH1D* hRho = (TH1D*)fin->Get(Form("%s%d", hrhobname.Data(), i));
		hRho->SetLineStyle(2);
		leg->AddEntry(hRho, Form("Data pT bin %d", i), "L");
		hRho->Draw("sames");
	}
	
	cCmp->cd();
	leg->Draw();
	
	SaveCv(cCmp);

	//TFile *finrespW = new TFile(weightedResp);
	//if(!finrespW->IsOpen()){
	//	Printf("Weighted response not found, quietly exit");
	//	return;
	//}

	
}

#endif
