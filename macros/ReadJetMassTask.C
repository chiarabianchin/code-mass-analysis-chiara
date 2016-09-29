#include "/data/macros/LoadALICEFigures.C"
#include "/data/Work/MyCodeJetMass/utils/CommonTools.C"
#include <TH3F.h>
#include <THnSparse.h>

TList* ReadMassTaskOutput(TString pathfile = "AnalysisResults.root", TString listname = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCRaw");

TObject* ReadMassHistogramFromList(TList *list, TString hnamebase = "fh3PtJet1VsMassVsLeadPt", TString selLev = "Tagged", Int_t centralityBin = 0);

TObject* ReadMassTHnSparseFromList(TList *list, TString hnamebase = "fh3PtJet1VsMassVsLeadPt", TString other = "", Int_t centralityBin = -1);

TH1D* ProjectTH3InRangeWithColorAndMarker(TH3F* h3, Int_t axrg1 = 0, Double_t min1 = 0, Double_t max1 = -1, Int_t axrg2 = 1, Double_t min2 = 0, Double_t max2 = -1, Int_t axpj = 2, Int_t color = 0, Int_t marker = 0, TString namepj = "hPj");

//_____________________________________________________________________________

TList* ReadMassTaskOutput(TString pathfile, TString listname){
	
	/// this method can read TH2 or TH3 from the input file
	/// Need to cast it properly in output... check if it works
	
	
	TList *list = ReadFile(pathfile, listname);
	return list;
}

//_____________________________________________________________________________

TObject* ReadMassHistogramFromList(TList *list, TString hnamebase, TString selLev, Int_t centralityBin){
	TH3F *h3read = 0x0;
	TH2F *h2read = 0x0;
	TString hname = Form("%s%s_%d", hnamebase.Data(), selLev.Data(), centralityBin);
	
	TObject* obj = list->FindObject(hname);
	TClass *cl = obj->IsA();
	TString clname = cl->GetName();
	
	if(clname == "TH2F"){
		h2read = (TH2F*)list->FindObject(hname);
		Printf("Reading TH2F %s", hname.Data());
		return h2read;
	}
	
	if(clname == "TH3F"){
		h3read = (TH3F*)list->FindObject(hname);
		Printf("Reading TH3F %s", hname.Data());
		return h3read;
	}
	Printf("%s is a %s, not supported", hname.Data(), clname.Data());
	return 0x0;
}

//_____________________________________________________________________________

TObject* ReadMassTHnSparseFromList(TList *list, TString hnamebase, TString other, Int_t centralityBin){
	
	THnSparseF *hsparseread = 0x0;
	TString hname = Form("%s%s", hnamebase.Data(), other.Data());
	if(centralityBin > -1) 	{
		hname+="_";
		hname+=centralityBin;
	}
	
	TObject* obj = list->FindObject(hname);
	TClass *cl = obj->IsA();
	TString clname = cl->GetName();
	
	if(clname == "ThnSparseF"){
		hsparseread = (THnSparseF*)list->FindObject(hname);
		Printf("Reading ThnSparseF %s", hname.Data());
		return hsparseread;
	}
	
	
	Printf("%s is a %s, not supported", hname.Data(), clname.Data());
	return 0x0;
}

//_____________________________________________________________________________

TH1D* ProjectTH3InRangeWithColorAndMarker(TH3F* h3, Int_t axrg1, Double_t min1, Double_t max1, Int_t axrg2, Double_t min2, Double_t max2, Int_t axpj, Int_t color, Int_t marker, TString namepj){
	
	/// axrg1 and axrg2 have to be in order 0 = x, 1 = y , 2 = z (skipping the axis that it the projection axis axpj)
	
	Int_t binrange1[2] = {0, -1};
	Int_t binrange2[2] = {0, -1};
	TH1D* hprojection = 0x0;
	Printf("Range axis %d (%.2f, %.2f) and %d (%.2f, %.2f)", axrg1, min1, max1, axrg2, min2, max2);
	if(axpj == 0) {
		if(min1 != 0)  binrange1[0] = h3->GetYaxis()->FindBin(min1);
		if(max1 != -1) binrange1[1] = h3->GetYaxis()->FindBin(max1);
		if(min2 != 0)  binrange2[0] = h3->GetZaxis()->FindBin(min2);
		if(max2 != -1) binrange2[1] = h3->GetZaxis()->FindBin(max2);
		hprojection = h3->ProjectionX(Form("%s_Pj%d", namepj.Data(), axpj), binrange1[0], binrange1[1], binrange2[0], binrange2[1]);
		Printf("(%d, %d) and (%d, %d)", binrange1[0], binrange1[1], binrange2[0], binrange2[1]);
	}
	if(axpj == 1) {
		if(min1 != 0)  binrange1[0] = h3->GetXaxis()->FindBin(min1);
		if(max1 != -1) binrange1[1] = h3->GetXaxis()->FindBin(max1);
		if(min2 != 0)  binrange2[0] = h3->GetZaxis()->FindBin(min2);
		if(max2 != -1) binrange2[1] = h3->GetZaxis()->FindBin(max2);
		hprojection = h3->ProjectionY(Form("%s_Pj%d", namepj.Data(), axpj), binrange1[0], binrange1[1], binrange2[0], binrange2[1]);
	}
	if(axpj == 2) {
		if(min1 != 0)  binrange1[0] = h3->GetXaxis()->FindBin(min1);
		if(max1 != -1) binrange1[1] = h3->GetXaxis()->FindBin(max1);
		if(min2 != 0)  binrange2[0] = h3->GetYaxis()->FindBin(min2);
		if(max2 != -1) binrange2[1] = h3->GetYaxis()->FindBin(max2);
		hprojection = h3->ProjectionZ(Form("%s_Pj%d", namepj.Data(), axpj), binrange1[0], binrange1[1], binrange2[0], binrange2[1]);
	}
	
	hprojection->SetLineColor(colors[color]);
	hprojection->SetLineWidth(2);
	hprojection->SetMarkerColor(colors[color]);
	Int_t markerstyle = color + 20;
	if(marker > 0) markerstyle = marker;
	hprojection->SetMarkerStyle(markerstyle);
	return hprojection;
}

//_____________________________________________________________________________

void OnlyOneTrigger(TString pathfile = "AnalysisResults.root", TString listname = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCRaw"){
	
	const Int_t nptbins = 4;
	Double_t ptbinlims[nptbins+1] = {40., 60., 80., 100., 120.};
	Int_t nx, ny, dx, dy;
	CalculatePads(nptbins, nx, ny, dx, dy);
	
	TCanvas *cMass = new TCanvas("cMass", "Mass", 800, 600);
	cMass->Divide(nx, ny);
	
	TList *list = ReadMassTaskOutput(pathfile, listname);
	if(!list) return;
	
	TFile *fout = new TFile("MassOutput.root", "recreate");
	TLegend *legCent = new TLegend(0.6, 0.4, 0.8, 0.7);
	legCent->SetFillStyle(0);
	legCent->SetBorderSize(0);
	
	for(Int_t ic = 0; ic < 4; ic++){
		TH3F * h3MassPtJPtTr = (TH3F*)ReadMassHistogramFromList(list, "fh3PtJet1VsMassVsLeadPt", "Tagged", ic);
		if(!h3MassPtJPtTr) {
			Printf("Histogram 3D not found");
		}
		fout->cd();
		h3MassPtJPtTr->Write();
		
		TH2D *h2DMassPt = (TH2D*)h3MassPtJPtTr->Project3D("yx");
		h2DMassPt->SetName(Form("hMassPtJ_cent%d", ic));
		fout->cd();
		h2DMassPt->Write();
		TCanvas *c2D = new TCanvas(Form("c2D%d", ic), Form("2D distib centrality bin %d", ic), 600, 500);
		h2DMassPt->Draw();
		for(Int_t ipt = 0; ipt < nptbins; ipt++){
			// mass (x), ptj (y), ptleadtr(z)
			TH1D * hPjM = ProjectTH3InRangeWithColorAndMarker(h3MassPtJPtTr, 0, ptbinlims[ipt], ptbinlims[ipt+1]-0.001, 2, 0, -1, 1, ic, 0, Form("hMassC%dPtB%d", ic, ipt));
			if(!hPjM) {
				Printf("Something went wrong with projection, continue");
				continue;
			}
			
			if(ipt == 0) legCent->AddEntry(hPjM, Form("CentBin %d", ic), "L");
			cMass->cd(ipt+1);
			gPad->SetLogy();
			if(ic == 0) {
				hPjM->Draw();
				TPaveText* pv = new TPaveText(0.35, 0.8, 0.75, 0.9, "NDC");
				pv->SetFillStyle(0);
				pv->SetBorderSize(0);
				pv->AddText(Form("%.0f < #it{p}_{T} (GeV/#it{c}) < %.0f", ptbinlims[ipt], ptbinlims[ipt+1]));
				pv->Draw();
			}
			else hPjM->Draw("sames");
			
			fout->cd();
			
			hPjM->Write();
		}
	}
	cMass->cd(1);
	legCent->Draw();
	SaveCv(cMass);
}


