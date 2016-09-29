#ifndef DrawComparisonsForPaper_C
#define DrawComparisonsForPaper_C
#include <TLegend.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLine.h>
#include <TMath.h>

//global
Double_t maxM = 26.;
//binning
const Int_t nptbins = 4;
Double_t ptlims[nptbins+1] = {40., 60., 80., 100., 120};
Double_t maxRangeMass[nptbins] = {20., 20., 26., 26.};


// Preferred colors and markers
const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed-9, kGreen-10, kBlue - 8};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

//definitions
TList* GetPbPbResults();
TList* GetPythiaPbPb();
void SaveCv(TCanvas* c, TString suffix = "", Int_t format = 2); 
Int_t UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type = "TH1F", Bool_t onlycheck = kFALSE);
// the first histogram contains the numerator and will be filled with the ratio and calculated error. The second histogram contains the denominator
void CalculateSysRatio(TH1D* hRatioPbPbOpPbSys, TH1D* hRDivpPbC);

//implementations
TList* GetPbPbResults(){
	TString pathMartaPbPb = "/data/Work/jets/JetMass/PbPbResults/TranformedIntoTH1.root";
	
	TFile *fResPbPbM = new TFile(pathMartaPbPb);
	if(!fResPbPbM->IsOpen()){
		Printf("File %s not found", pathMartaPbPb.Data());
		return 0;
	}

	TString nameMSt = "hUnfMass_PtBin";
	TString nameMSy = "hUnfMassSyst_PtBin";
	const Int_t nhM = 4;
	
	TList *listPbPb = new TList();
	listPbPb->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fResPbPbM->Get(Form("%s%d", nameMSt.Data(), ih+1));
		TH1D* hMassSy = (TH1D*)fResPbPbM->Get(Form("%s%d", nameMSy.Data(),ih+1));
		
		// set this nicely at some point
		hMassSt->SetMarkerColor(kBlue+2);
		hMassSt->SetLineColor(kBlue+2);
		hMassSy->SetFillColor(kBlue-10);
		hMassSy->SetMarkerColor(kBlue+2);
		hMassSy->GetXaxis()->SetRangeUser(0., maxM);
		
		listPbPb->Add(hMassSt);
		listPbPb->Add(hMassSy);
	}
	
	return listPbPb;
}

//________________________________________________________________________

TList* GetPythiaPbPb(){
	
	TString pathMartaPythia276 = "/data/Work/jets/JetMass/PbPbResults/JetMassPerugia2011_ecms2760.root";
	
	TFile *fMPythia276 = new TFile(pathMartaPythia276);
	if(!fMPythia276->IsOpen()){
		Printf("File %s not found", pathMartaPythia276.Data());
		return 0;
	}

	TString nameMSt = "hM_";
	const Int_t nhM = 4;
	Int_t offset = 1;
	
	TList *listPythia = new TList();
	listPythia->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fMPythia276->Get(Form("%s%d", nameMSt.Data(), ih+offset));
		
		// set this nicely at some point
		hMassSt->SetMarkerStyle(25);
		hMassSt->SetMarkerColor(kRed+2);
		hMassSt->SetLineColor(kRed+2);
		
		listPythia->Add(hMassSt);
	}
	
	return listPythia;
}


//________________________________________________________________________
void CalculateSysRatio(TH1D* hRatioPbPbOpPbSysC, TH1D* hRDivpPbC){
	
	// the first histogram contains the PbPb and will be filled with the ratio
	
	for(Int_t ib = 0; ib<hRatioPbPbOpPbSysC->GetNbinsX(); ib++) {
		Double_t errPbPb = hRatioPbPbOpPbSysC->GetBinError(ib+1), valPbPb =  hRatioPbPbOpPbSysC->GetBinContent(ib+1);
		Double_t errpPb = hRDivpPbC->GetBinError(ib+1), valpPb =  hRDivpPbC->GetBinContent(ib+1);
		//|| TMath::Abs(errPbPb/valPbPb) > 1e3 || TMath::Abs(errpPb/valpPb) > 1e3
		if(valpPb<0 || valPbPb<0 || TMath::Abs(valpPb) < 1e-6 || TMath::Abs(valPbPb) < 1e-6 ) {
			hRatioPbPbOpPbSysC->SetBinContent(ib+1, 0);
			hRatioPbPbOpPbSysC->SetBinError(ib+1, 0);
		} 
		//Printf("bin %d rel err 1 = %f , rel err2 = %f ", ib+1, errPbPb/valPbPb, errpPb/valpPb);
		Double_t sysC = valPbPb/valpPb * TMath::Sqrt((errPbPb/valPbPb * errPbPb/valPbPb) + (errpPb/valpPb * errpPb/valpPb));
		//Printf("Filling with %f +- %f", valPbPb/valpPb, sysC);
		
		hRatioPbPbOpPbSysC->SetBinContent(ib+1, valPbPb/valpPb);
		hRatioPbPbOpPbSysC->SetBinError(ib+1, sysC);
		//Printf(" =  %f +- %f", hRatioPbPbOpPbSysC->GetBinContent(ib+1), hRatioPbPbOpPbSysC->GetBinError(ib+1));
	}
}

//________________________________________________________________________

// mains
void DrawPbPb(){
	TString suff = "";
	
	//get the PbPb results
	TList* listPbPb = GetPbPbResults();
		
	//get Pythia PbPb
	TList* listPythia276 = GetPythiaPbPb();
	
	if(!listPbPb){ 
		Printf("List PbPb not found");
		return;
	}
	
	
	if(!listPythia276){
		Printf("List Pythia for PbPb energy not found");
		return;
	}
	
	const Int_t nhM = 4;
	
	
	TCanvas *cMassPbPb = new TCanvas(Form("cMassPbPb%s", suff.Data()), "Mass PbPb compared with Pythia", 1100, 800);
	cMassPbPb->Divide(2, 2);
	TLegend *legMassPbPb = new TLegend(0.5, 0.5, 0.8, 0.9);
	legMassPbPb->SetBorderSize(0);
	legMassPbPb->SetFillStyle(0);
	
	TCanvas *cRatioPbPbOPy = new TCanvas(Form("cRatioPbPbOPy%s", suff.Data()), "Ratio PbPb over Pythia 2.76 TeV", 1100, 800);
	cRatioPbPbOPy->Divide(2, 2);
	TLegend *legRatioPbPbOPy = new TLegend(0.4, 0.5, 0.8, 0.9);
	legRatioPbPbOPy->SetBorderSize(0);
	legRatioPbPbOPy->SetFillStyle(0);
	
	Int_t n = 2;
	
		
	for(Int_t ih = 0; ih<nhM ; ih++){
		
		//mass
		TH1D* hMPbPb = (TH1D*)listPbPb->At(n*ih);
		
		TH1D* hMPy276= (TH1D*)listPythia276->At(ih);
		
		hMPy276->Scale(1./hMPy276->Integral("width"));
		
		//syst
		TH1D* hSyPbPb = (TH1D*)listPbPb->At(n*ih+1);
		
		
		//Draw comparison PbPb Pythia 2.76 TeV
		cMassPbPb->cd(ih+1);
		
		hSyPbPb->Draw("E2");
		hMPbPb ->Draw("sames");
		hMPy276->Draw("sames");
		if(ih == 0){
			legMassPbPb->AddEntry(hMPbPb, "Mass PbPb", "LP");
			legMassPbPb->AddEntry(hSyPbPb, "Sys PbPb", "F");
			legMassPbPb->AddEntry(hMPy276, "PYTHIA 2.76 TeV", "LP");
			legMassPbPb->Draw();
		}
		
		
		//line for reference 
		TLine *lineOne = new TLine(0, 1., maxRangeMass[ih], 1.);
		lineOne->SetLineStyle(2);
		lineOne->SetLineWidth(2);
		lineOne->SetLineColor(kGray);
		
		
		//ratio PbPb/PYTHIA
		TH1* hRatioPbPbOPy;
		TH1* hDividePy276;
		UniformTH1FForDivide(hMPbPb, hMPy276, hRatioPbPbOPy, hDividePy276, "TH1D"); hRatioPbPbOPy->SetName(Form("hRatiopPbOPy_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		hRatioPbPbOPy->Divide(hDividePy276);
		hRatioPbPbOPy->GetXaxis()->SetRangeUser(0., maxM);
		hRatioPbPbOPy->GetYaxis()->SetRangeUser(0., 4);
		
		cRatioPbPbOPy->cd(ih+1);
		hRatioPbPbOPy->Draw();
		lineOne->Draw();
		
		// sys ratio PbPb/PYTHIA
		TH1* hRatioPbPbOPy276Sys;
		UniformTH1FForDivide(hSyPbPb, hMPy276, hRatioPbPbOPy276Sys, hDividePy276, "TH1D");
		TH1* hDividePy276NoErr = (TH1*)hDividePy276->Clone(Form("hDividePy276NoErr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		for(Int_t ib = 0; ib<hDividePy276NoErr->GetNbinsX(); ib++) hDividePy276NoErr->SetBinError(ib, 0);
		
		hRatioPbPbOPy276Sys->SetName(Form("hRatioPbPbOPy276Sys_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatioPbPbOPy276Sys->Divide(hDividePy276NoErr);
		hRatioPbPbOPy276Sys->SetFillStyle(1);
		hRatioPbPbOPy276Sys->SetLineWidth(2);
		hRatioPbPbOPy276Sys->SetLineColor(kBlue-10);
		hRatioPbPbOPy276Sys->SetFillColor(kBlue-10);
		
		cRatioPbPbOPy->cd(ih+1);
		
		hRatioPbPbOPy276Sys->Draw("E2sames");
		if(ih == 3){
			legRatioPbPbOPy->AddEntry(hRatioPbPbOPy, "PbPb/PYTHIA2.76TeV", "PL");
			legRatioPbPbOPy->AddEntry(hRatioPbPbOPy276Sys, "(Sys pPb)/PYTHIA2.76TeV", "F");
			legRatioPbPbOPy->Draw();
		}
		
		
	}
	
	
	SaveCv(cMassPbPb);
	
	SaveCv(cRatioPbPbOPy);
}

//_____________________________________________________________________________

Int_t UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type, Bool_t onlycheck){
   
   /// returns -1 if input is not valid, 0 if no change in the histograms is needed, 1 if histograms are redefined
   
   if(!h1 || !h2){
      Printf("Need valid input!");
      return -1;
   }
   Printf("UniformTH1FForDivide:: THINK: if the histograms are ratios this method does not do what you want");
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h2->GetNbinsX();
   Double_t binW1 = h1->GetBinWidth(3), binW2 = h2->GetBinWidth(3);
   Double_t low1 = h1->GetBinLowEdge(1), low2 = h2->GetBinLowEdge(1);
   Double_t up1 = h1->GetBinLowEdge(nBins1+1), up2 = h2->GetBinLowEdge(nBins2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBins2, binW2, low2, up2);
   if(((nBins1 == nBins2) && ((binW1-binW2) < 1e-6) &&  ((low1-low2) < 1e-6) && ((up1-up2) < 1e-6)) || onlycheck) {
      if(onlycheck) Printf("Check performed");
      else Printf("Histograms can be divided");
      h1bis = (TH1*)h1->Clone(Form("%sb", h1->GetName()));
      h2bis = (TH1*)h2->Clone(Form("%sb", h2->GetName()));
      return 0;
   }
   
   Double_t commonlow = TMath::Min(low1, low2), commonup = TMath::Max(up1, up2);
   Double_t commonW = TMath::Max(binW1, binW2); // the wider bin width it the one that will be used for both
   Int_t commonNBins = (commonup - commonlow)/commonW;
   Printf("Common) N bins = %d, Width = %f, Range = %f - %f", commonNBins, commonW, commonlow, commonup);
   Double_t minWidth = TMath::Min(binW1, binW2); // pick the smaller among the two bin widths...
      
   Int_t rebin = (Int_t)(commonW/minWidth); //... and calculate the number of bins to be combined
   Printf("%f/%f = %f or %d", commonW, minWidth, commonW/minWidth, rebin);
   
   if(minWidth == binW1) h1->Rebin(rebin); //if the smaller width was for h1, rebin h1
   else h2->Rebin(rebin); //otherwise rebin h2
   
   
   //TCanvas *ctest = new TCanvas(Form("c%s%s", h1->GetName(), h2->GetName() ), Form("canvas %s%s", h1->GetName(), h2->GetName() ));
   //ctest->cd();
   //h1->Draw();
   //h2->Draw("sames");
   
   //Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsX(), h2->GetName(), h2->GetNbinsX());
   //create the new histograms with correct range
   if(type == "TH1F"){
      h1bis = new TH1F(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      h2bis = new TH1F(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
   } else {
      if(type == "TH1D"){
      	 h1bis = new TH1D(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      	 h2bis = new TH1D(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      } else {
      	 Printf("%s Not defined!", type.Data());
      	 return -1;
      }
   }
   h1bis->Sumw2();
   h1bis->SetLineColor(h1->GetLineColor());
   h1bis->SetLineWidth(h1->GetLineWidth());
   h1bis->SetMarkerColor(h1->GetMarkerColor());
   h1bis->SetMarkerStyle(h1->GetMarkerStyle());
   
   h2bis->Sumw2();
   h2bis->SetLineColor(h2->GetLineColor());
   h2bis->SetLineWidth(h2->GetLineWidth());
   h2bis->SetMarkerColor(h2->GetMarkerColor());
   h2bis->SetMarkerStyle(h2->GetMarkerStyle());
   
   //refill the original histograms using the same binning
   for(Int_t i = 0; i < commonNBins; i++){
      
      //First histogram
      Float_t bincentre = h1bis ->GetBinCenter(i+1);
      Int_t fillingbin = h1->FindBin(bincentre);
      Double_t fillWith = 0, error = 0;
      if(fillingbin>0 && fillingbin <= nBins1) { // exclude over- / under-flow bins
      	 //Printf("Here, debug %d", nBins1);
      	 fillWith = (Double_t)h1->GetBinContent(fillingbin);
      	 error = (Double_t)h1->GetBinError(fillingbin);
      }
      h1bis ->SetBinContent(i+1, fillWith);
      h1bis ->SetBinError(i+1, error);
      if(i%10 == 0) Printf("Filling %s x = %.3f new bin %d with old bin %d content %e [...]", h1->GetName(), bincentre, i+1, fillingbin, fillWith);
      
      //second histogram
      bincentre = h2bis ->GetBinCenter(i+1);
      fillingbin = h2->FindBin(bincentre);
      fillWith = 0; error = 0;
      if(fillingbin>0 && fillingbin <= nBins2) {
      	 fillWith = h2->GetBinContent(fillingbin);
      	 error = h2->GetBinError(fillingbin);
      }
      h2bis ->SetBinContent(i+1, fillWith);
      h2bis ->SetBinError(i+1, error);
       if(i%10 == 0) Printf("Filling %s x = %.3f new bin %d with old bin %d content %e [...]", h2->GetName(), bincentre, i+1, fillingbin, fillWith);
   }
   return 1;
   
}
//_____________________________________________________________________________

void SaveCv(TCanvas* c,TString suffix, Int_t format){
   if(format > 0) c->SaveAs(Form("%s%s.png",c->GetName(),suffix.Data()));
   if(format > 1) c->SaveAs(Form("%s%s.pdf",c->GetName(),suffix.Data()));
   if(format > 2) c->SaveAs(Form("%s%s.eps",c->GetName(),suffix.Data()));
   if(format > 3 || format==0) c->SaveAs(Form("%s%s.root",c->GetName(),suffix.Data()));
}
#endif
