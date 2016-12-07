#ifndef DrawComparisonsForPaper_C
#define DrawComparisonsForPaper_C
#include "/data/macros/LoadALICEFigures.C"
#include "/data/Work/code-mass-analysis-chiara/utils/CommonTools.C"
#include "/data/Work/code-mass-analysis-chiara/unfold/SystematicComparisonsUnf.C"
#include <TLegend.h>
//#include <TGaxis.h>


// final figures with RatiopPbPbPb, 
//global

Double_t maxM = 26.;
//definitions
TList* GetPbPbResults();
TList* GetPbPbResultsGraph();
TList* GetpPbResults(Bool_t kinecorr = kTRUE);
TList* GetpPbMarta();
TList* GetPythiapPb(Int_t input = 0);
TList* GetPythiapPbMarta(Bool_t useline = kFALSE);
TList* GetPythiaPbPb(Bool_t useline = kFALSE);
TList* GetJEWELPbPb(Bool_t useline = kFALSE);
TList* GetJEWELPbPbrecoilOff(Bool_t useline = kFALSE);
TList* GetJEWELpp(Bool_t useline = kFALSE);
TList* GetJEWEL(Int_t system, Bool_t useline = kFALSE);
TH1*   Marco_get_shape_hist(TH3F* h_pt_eta_shape, Float_t ptmin, Float_t ptmax, const char *hname = "hshape_pt", Int_t nrebin = 1);
TList* GetJEWELQPYTHIA(Int_t mod, Bool_t useline = kFALSE);
TList* GetHERWIG(Int_t system, Bool_t useline = kFALSE);


TH1D* CalculateMean(Double_t **xrange, Int_t nbins, Double_t lims[], TH1D** hdistr, TGraphErrors *&grmean, TH1D *&hrms, const char* namehmean);
void DrawMeanComparison(Double_t **xrange, const Int_t ninputs, TString inputDistr[], TString h1names[], const Int_t ninputsSy, TString inputSyst[], TString hsysnames[], Int_t offset[], TString leg[], Int_t mrk[], Int_t clr[], Int_t clrsy[], Int_t fillsy[]);
TH1D* CalculateSysMean(Double_t **xrange, Int_t nbins, Double_t lims[], TH1D** hSysdistr, TGraphErrors *&grSysmean, TH1D *&hSysrms, const char* namehsysmean);
TH1D* GetSystFromFile(TFile *fin);
TH1D* GetMeanFromFile(TString filename);
TH1D* GetMeanFromFile(TFile *fin);
TH1D* GetMeanSystFromFile(TString filename);
TH1D* GetMeanAndSystFromFile(TString filename, TH1D*& hSyst);

TH1D* SystEnergyDep(TString filenamepPb);
TList* TreatSystematicsPbPb(Int_t* reject, TString name);

TH1D* GetpPb(TFile* finpPb, Int_t index);
void GetpPbSystematics(Int_t* reject, TString name, TFile* finpPb, Int_t index,  TH1D*& hNewSyspPb, TH1D*& hFullSyspPb);
void RatioWithReducedSyst(Int_t* rejectPbPb, Int_t* rejectpPb, TString name);


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
		TH1D* hMassSt = (TH1D*)fResPbPbM->Get(TString::Format("%s%d", nameMSt.Data(), ih+1));
		TH1D* hMassSy = (TH1D*)fResPbPbM->Get(TString::Format("%s%d", nameMSy.Data(),ih+1));
		
		// set this nicely at some point
		hMassSt->SetMarkerColor(kBlue+2);
		hMassSt->SetMarkerSize(2);
		hMassSt->SetLineColor(kBlue+2);
		
		hMassSy->SetFillColor(kBlue-10);
		hMassSy->SetMarkerStyle(1);
		hMassSy->SetMarkerSize(2);
		hMassSy->SetMarkerColor(hMassSy->GetFillColor());
		hMassSy->SetLineColor(hMassSy->GetFillColor());
		hMassSy->GetXaxis()->SetRangeUser(0., maxM);
		
		hMassSy->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		//hMassSy->GetYaxis()->SetTitleOffset(1.7);
		hMassSt->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		//hMassSy->GetYaxis()->SetTitleOffset(1.7);
		listPbPb->Add(hMassSt);
		listPbPb->Add(hMassSy);
	}
	
	return listPbPb;
}

TList* GetPbPbResultsGraph(){
	TString pathMartaPbPb = "/data/Work/jets/JetMass/PbPbResults/DistJetMassAllSyst_Area.root";
	
	TFile *fResPbPbM = new TFile(pathMartaPbPb);
	if(!fResPbPbM->IsOpen()){
		Printf("File %s not found", pathMartaPbPb.Data());
		return 0;
	}

	TString nameMSt = "grUnfMass_PtBin";
	TString nameMSy = "grUnfMassSyst_PtBin";
	const Int_t nhM = 4;
	Int_t offset = 1;
	
	TList *listPbPb = new TList();
	listPbPb->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		
		TGraphErrors* gMassSt = (TGraphErrors*)fResPbPbM->Get(TString::Format("%s%d", nameMSt.Data(), ih+offset));
		TGraphErrors* gMassSy = (TGraphErrors*)fResPbPbM->Get(TString::Format("%s%d", nameMSy.Data(),ih+offset));
		
		Int_t npoints = gMassSt->GetN();
		if(gMassSy->GetN() != npoints) Printf("Number of points in syst = %d != %d",gMassSy->GetN(), npoints);
		Double_t binlimits[npoints+1];
		Double_t yvalue[npoints];
		Double_t yerror[npoints];
		Double_t ysyste[npoints];
		Double_t last = 0;
		for(Int_t j = 0; j<npoints; j++){
			Double_t x;
			gMassSt->GetPoint(j, x, yvalue[j]);
			binlimits[j] = x - gMassSt->GetErrorXlow(j);
			last = x + gMassSt->GetErrorXhigh(j);
			yerror[j]    = gMassSt->GetErrorY(j);
			ysyste[j]    = gMassSy->GetErrorY(j);
			
			//Printf("binlims[%d] = %.1f, y = %f+-%f+-%f", j, binlimits[j], yvalue[j], yerror[j], ysyste[j]);
		
		}
		binlimits[npoints] = last;
		
		TH1D* hMassSt = new TH1D(TString::Format("hUnfMass_PtBin%d", ih+offset), TString::Format("hUnfMass_PtBin%d", ih+offset), npoints, binlimits);
		hMassSt->Sumw2();
		TH1D* hMassSy = new TH1D(TString::Format("hUnfMassSyst_PtBin%d", ih+offset), TString::Format("hUnfMassSyst_PtBin%d", ih+offset), npoints, binlimits);
		hMassSy->Sumw2();
		// set this nicely at some point
		hMassSt->SetMarkerColor(kBlue+2);
		hMassSt->SetMarkerStyle(20);
		hMassSt->SetMarkerSize(2);
		hMassSt->SetLineColor(kBlue+2);
		
		hMassSy->SetFillColor(kBlue-10);
		hMassSy->SetMarkerStyle(1);
		hMassSy->SetMarkerSize(2);
		hMassSy->SetMarkerColor(hMassSy->GetFillColor());
		hMassSy->SetLineColor(hMassSy->GetFillColor());
		hMassSy->GetXaxis()->SetRangeUser(0., maxM);
		
		hMassSy->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		//hMassSy->GetYaxis()->SetTitleOffset(1.7);
		hMassSt->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		//hMassSy->GetYaxis()->SetTitleOffset(1.7);
		for(Int_t j = 0; j<npoints; j++){
			hMassSt->SetBinContent(j+1, yvalue[j]);
			hMassSt->SetBinError  (j+1, yerror[j]);
			hMassSy->SetBinContent(j+1, yvalue[j]);
			hMassSy->SetBinError  (j+1, ysyste[j]);
		
		}
		
		listPbPb->Add(hMassSt);
		listPbPb->Add(hMassSy);
	}
	
	return listPbPb;
}

//________________________________________________________________________

TList* GetpPbResults(Bool_t kinecorr){
	
	if(!kinecorr) Printf("DOESN'T WORK PROPERLY YET");
	
	TString pathResults = 
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/AllSyst20160819/MasspPbResults.root";
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160916/MasspPbResults.root";
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreAppr/MasspPbResults.root";
	TString pathSyst = 
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/AllSyst20160819/TotalSystematicUnc.root";
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160916/TotalSystematicUnc.root";
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreApprTotalSystematicUnc.root";
	
	TFile *fin = new TFile(pathResults);
	if(!fin->IsOpen()){
		Printf("File %s not found", pathResults.Data());
		return 0;
	}
	//TFile *fSys = new TFile(pathSyst);
	//if(!fSys->IsOpen()){
	//	Printf("File %s not found", pathSyst.Data());
	//	return 0;
	//}
	
	TString nameMSt = "hUnfEffCorM_Itr3_Pt";
	if(!kinecorr) nameMSt = "hUnfM_Itr3_Pt";
	TString nameMSy = "hSystTot_Pt";
	const Int_t nhM = nptbins;
	Printf("Number of bins found: %d", nhM);
	TList *listpPb = new TList();
	listpPb->SetOwner();
	
	// apply the systematics on the mass corrected for kine eff
	TList *listKineEffSy = 0x0;
	if(kinecorr) listKineEffSy = AddSystematicstoMassFromFile(pathResults, nameMSt, kTRUE, pathSyst, nameMSy, kTRUE);
	
	//listKineEffSy->ls();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		Int_t nbinsM = (Int_t)(maxRangeMassFinal[ih]/2.);
		Printf("Filling histogram %d with %d bins from 0 to %f", ih, nbinsM, maxRangeMassFinal[ih]);
		
		TH1D* hMassSt = new TH1D(TString::Format("%sMax%.0f_%.0f", nameMSt.Data(), ptlims[ih],ptlims[ih+1]), TString::Format("Mass pPb Pt %.0f - %.0f; #it{M} (GeV/c^2)", ptlims[ih],ptlims[ih+1]), nbinsM, 0, maxRangeMassFinal[ih]);
		hMassSt->Sumw2();
		hMassSt->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		
		TH1D* hMassSy = (TH1D*)hMassSt->Clone(TString::Format("%sMax%.0f_%.0f", nameMSy.Data(), ptlims[ih],ptlims[ih+1]));
		hMassSy->Sumw2();
		
		TH1D* hMtmp = (TH1D*)fin->Get(TString::Format("%s%.0f_%.0f", nameMSt.Data(), ptlims[ih],ptlims[ih+1]));
		if(!hMtmp) {
			Printf("Problem, mass histogram %s%.0f_%.0f not found", nameMSt.Data(), ptlims[ih],ptlims[ih+1]);
			fin->ls();
		}
		TH1D* hSytmp = 0x0;
		if(kinecorr) hSytmp = (TH1D*)listKineEffSy->At(ih*nhM+1);
		else hSytmp = (TH1D*)fin->Get(TString::Format("%s%.0f_%.0fb", nameMSy.Data(), ptlims[ih],ptlims[ih+1]));
		Printf("Low edge %f and %f (they have to be the same)", hMtmp->GetBinLowEdge(0), hSytmp->GetBinLowEdge(0));
		
		for(Int_t ib = 0; ib < nbinsM; ib++){
			Int_t bin = hMtmp->FindBin(hMassSt->GetBinCenter(ib+1));
			//Printf("BIN = %d", bin);
			if(bin == 0 || bin == -1) {
				hMassSt->SetBinContent(ib+1, 0);
				hMassSt->SetBinError(ib+1, 0);
				hMassSy->SetBinContent(ib+1, 0);
				hMassSy->SetBinError(ib+1, 0);
			} else {
				hMassSt->SetBinContent(ib+1, hMtmp->GetBinContent(bin));
				hMassSt->SetBinError(ib+1, hMtmp->GetBinError(bin));
				hMassSy->SetBinContent(ib+1, hSytmp->GetBinContent(bin));
				hMassSy->SetBinError(ib+1, hSytmp->GetBinError(bin));
			}
		}
		
		// set this nicely at some point
		hMassSt->SetMarkerStyle(hMtmp->GetMarkerStyle());
		hMassSt->SetMarkerColor(hMtmp->GetMarkerColor());
		hMassSt->SetLineColor(hMtmp->GetLineColor());
		hMassSt->SetMarkerSize(2);
		//hMassSt->GetYaxis()->SetTitleOffset(1.7);
		
		//hMassSt->GetXaxis()->SetRangeUser(0., maxRangeMassFinal[ih]);
		hMassSy->SetFillStyle(hSytmp->GetFillStyle());
		hMassSy->SetLineWidth(hSytmp->GetLineWidth());
		hMassSy->SetLineColor(hSytmp->GetLineColor());
		hMassSy->SetMarkerSize(2);
		//hMassSy->GetYaxis()->SetTitleOffset(1.7);
		//hMassSy->GetXaxis()->SetRangeUser(0., maxRangeMassFinal[ih]);
		
		listpPb->Add(hMassSt);
		listpPb->Add(hMassSy);
		
	}
	return listpPb;
}

//________________________________________________________________________
TList* GetpPbMarta(){

	TString pathMartapPb = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root";
	
	TFile *fRespPbM = new TFile(pathMartapPb);
	if(!fRespPbM->IsOpen()){
		Printf("File %s not found", pathMartapPb.Data());
		return 0;
	}

	TString nameMSt = "hUnfMass_PtBin";
	TString nameMSy = "hUnfMassSyst_PtBin";
	const Int_t nhM = 4;
	
	TList *listpPb = new TList();
	listpPb->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fRespPbM->Get(TString::Format("%s%d", nameMSt.Data(), ih+1));
		TH1D* hMassSy = (TH1D*)fRespPbM->Get(TString::Format("%s%d", nameMSy.Data(),ih+1));
		
		// set this nicely at some point
		hMassSt->SetLineColor(kOrange+1);
		hMassSt->SetMarkerColor(kOrange+1);
		hMassSt->SetMarkerSize(2);
		
		hMassSy->SetLineColor(kOrange+1);
		hMassSy->SetFillColor(kOrange+1);
		hMassSy->SetFillStyle(3001);
		hMassSy->GetXaxis()->SetRangeUser(0., maxM);
		
		listpPb->Add(hMassSt);
		listpPb->Add(hMassSy);
	}
	
	return listpPb;

}

//________________________________________________________________________

TList* GetPythiapPb(Int_t input){
	
	// can add other input files here and assign a value for "input" (default = 0). Assign also projection axis and type of projection method (from THnSparse of TH3F) accordingly
	
	// PYTHIA 
	TString listnamePythia = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme";
	TString pathPythia = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/AnalysisResults.root";
	TString thnspname = "fhnMassResponse", projnamepythia = "hMPythiaPar";
	
	Int_t axrange = 3, axproj = 1;
	
	Double_t massbinW = 2.;
	TH1D** hMassPythia = 0x0;
	
	if(input == 0) hMassPythia = GetPythiaOrThnSpaseProjections(axrange, axproj, massbinW, nptbins, ptlims, pathPythia, listnamePythia, thnspname, projnamepythia);
	
	if(input == 1) {
		listnamePythia = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCMerged";
		pathPythia = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root";
		axrange = 0;
		axproj = 1;
		thnspname = "";
		projnamepythia = "hMPythiaDet";
		hMassPythia = GetProjectionsTH3(axrange, axproj, massbinW, nptbins, ptlims, pathPythia, listnamePythia, thnspname, projnamepythia);
		
	}
	if(!hMassPythia) Printf("Not possible to find pythia");
	
	TList *listPythia = new TList();
	listPythia->SetOwner();
	
	for(Int_t ih = 0; ih<nptbins; ih++){
		hMassPythia[ih]->SetMarkerStyle(25);
		hMassPythia[ih]->SetMarkerSize(2);
		hMassPythia[ih]->SetMarkerColor(kGreen-3);
		hMassPythia[ih]->SetLineColor(kGreen-3);
		listPythia->Add(hMassPythia[ih]);
	}
	return listPythia;
}

//________________________________________________________________________

TList* GetPythiaPbPb(Bool_t useline){
	
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
		TH1D* hMassSt = (TH1D*)fMPythia276->Get(TString::Format("%s%d", nameMSt.Data(), ih+offset));
		
		// set this nicely at some point
		hMassSt->SetMarkerStyle(25);
		hMassSt->SetMarkerSize(2);
		hMassSt->SetMarkerColor(kBlue-7);
		hMassSt->SetLineColor(hMassSt->GetMarkerColor());
		//hMassSt->SetLineColor(kBlack);
		if(useline) {
			hMassSt->SetLineStyle(1);
			hMassSt->SetFillColor(hMassSt->GetLineColor());
			hMassSt->SetMarkerStyle(1);
			hMassSt->SetFillStyle(0);
			hMassSt->SetLineWidth(3);
		}
		hMassSt->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		
		listPythia->Add(hMassSt);
	}
	
	return listPythia;
}

//________________________________________________________________________

TList* GetPythiapPbMarta(Bool_t useline){
	
	TString pathMartaPythia502 = "/data/Work/jets/JetMass/DetectorCorrections/Marta/JetMassPerugia2011_ecms5020.root";
	
	TFile *fMPythia502 = new TFile(pathMartaPythia502);
	if(!fMPythia502->IsOpen()){
		Printf("File %s not found", pathMartaPythia502.Data());
		return 0;
	}

	TString nameMSt = "hM_";
	const Int_t nhM = 4;
	Int_t offset = 1;
	
	TList *listPythia = new TList();
	listPythia->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fMPythia502->Get(TString::Format("%s%d", nameMSt.Data(), ih+offset));
		
		// set this nicely at some point
		hMassSt->SetMarkerStyle(25);
		hMassSt->SetMarkerColor(kBlack);
		hMassSt->SetLineColor(kBlack);
		if(useline) {
			hMassSt->SetLineStyle(1);
			hMassSt->SetFillColor(hMassSt->GetLineColor());
			hMassSt->SetMarkerStyle(1);
			hMassSt->SetFillStyle(1);
			hMassSt->SetLineWidth(3);
		}
		hMassSt->SetMarkerSize(2);
		listPythia->Add(hMassSt);
	}
	
	return listPythia;
}

//________________________________________________________________________

TList* GetJEWELPbPb(Bool_t useline){
	return GetJEWEL(1, useline);
}

//________________________________________________________________________

TList* GetJEWELPbPbrecoilOff(Bool_t useline){
	return GetJEWEL(2, useline);
}

//________________________________________________________________________

TList* GetJEWELpp(Bool_t useline){
	return GetJEWEL(0, useline);
}

//________________________________________________________________________
TList* GetJEWEL(Int_t system, Bool_t useline){
	//system: 0 = pp, 1 = 0-10% Pb-Pb, 2 = 0-10% recoils off
	
	TString pathJewel = "/data/Work/jets/JetMass/PbPbJEWEL/alice_jetmass_histograms_chargedJetsShift_withoutsys.root";
	if(system == 2) pathJewel =	"/data/Work/jets/JetMass/PbPbJEWEL/alice_jetmass_norecoils_histograms.root";
	TFile *fJewel = new TFile(pathJewel);
	if(!fJewel->IsOpen()){
		Printf("File %s not found", pathJewel.Data());
		return 0;
	}
	Int_t cent = 0;
	if(system > 0) cent = 1;
	TString nameMSt = TString::Format("Mass_TrackShifted_centbin%d_ptbin", cent);
	TString nameRat = "Ratio_centbin1_ptbin";
	const Int_t nhM = 3;
	Int_t offset = 0;
	
	TList *listJewel = new TList();
	listJewel->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fJewel->Get(TString::Format("%s%d", nameMSt.Data(), ih+offset));
		TH1D* hMassRa = (TH1D*)fJewel->Get(TString::Format("%s%d", nameRat.Data(), ih+offset));
		
		// set this nicely at some point
		hMassSt->SetMarkerStyle(28);
		hMassSt->SetLineWidth(3);
		if(system == 0){
			
			hMassSt->SetMarkerColor(kGreen+2);
			
			if(useline) {
				hMassSt->SetLineStyle(5);
				hMassSt->SetMarkerStyle(1);
				hMassSt->SetFillStyle(3004);
			}
		}
		if(system == 1) {
			hMassSt->SetMarkerColor(kOrange+1);
			if(useline) {
				hMassSt->SetLineStyle(5);
				
				hMassSt->SetMarkerStyle(1);
				hMassSt->SetFillStyle(3005);
			}
		}
		if(system == 2) {
			hMassSt->SetMarkerColor(kRed);
			if(useline) {
				hMassSt->SetLineStyle(3);
				hMassSt->SetFillStyle(3004);
				hMassSt->SetMarkerStyle(1);
			}
		}
		hMassSt->SetLineColor(hMassSt->GetMarkerColor());
		hMassSt->SetFillColor(hMassSt->GetLineColor());
		listJewel->Add(hMassSt);
		hMassSt->SetMarkerSize(2);
		
		if(hMassRa) {
			hMassRa->SetMarkerStyle(hMassSt->GetMarkerStyle());
			hMassRa->SetMarkerColor(hMassSt->GetMarkerColor());
			hMassRa->SetLineColor(hMassSt->GetMarkerColor());
			listJewel->Add(hMassRa);
		}
		
		
	}
	
	return listJewel;

}

//________________________________________________________________________

TH1 *Marco_get_shape_hist(TH3F* h_pt_eta_shape, Float_t ptmin, Float_t ptmax, const char *hname, Int_t nrebin) {
	
	Int_t min_eta = -0.9, max_eta = 0.9;
	
	Int_t minptbin = h_pt_eta_shape->GetXaxis()->FindBin(ptmin+0.001);
	Int_t maxptbin = h_pt_eta_shape->GetXaxis()->FindBin(ptmax-0.001);
	h_pt_eta_shape->GetXaxis()->SetRange(minptbin, maxptbin);
	Int_t minetabin = h_pt_eta_shape->GetYaxis()->FindBin(min_eta+0.001);
	Int_t maxetabin = h_pt_eta_shape->GetYaxis()->FindBin(max_eta-0.001);
	h_pt_eta_shape->GetYaxis()->SetRange(minetabin, maxetabin);
	TH1 *hshape = h_pt_eta_shape->Project3D("z");
	hshape->UseCurrentStyle();
	hshape->SetName(hname);
	hshape->SetLineWidth(2);
	hshape->SetMarkerStyle(20);
	hshape->Rebin(nrebin);
	hshape->Scale(1./hshape->Integral(),"width");
	return hshape;
}

//________________________________________________________________________

TList* GetJEWELQPYTHIA(Int_t mod, Bool_t useline){
	
	//mod: 1 = JEWEL recoil On, 2 = JEWEL recoil Off, 3 = QPYTHIA, 4 = Vacuum
	
	TString pathModels = "/data/Work/jets/JetMass/PbPbJEWEL/Marco/JetMassJEWELR040.root";
	if(mod == 3) {
		
		pathModels = "/data/Work/jets/JetMass/PbPbJEWEL/Marco/JetMassQPYTHIAv11.root";
		
		TFile *ftest = new TFile(pathModels);
		if(!ftest->IsOpen()){
			// make the projections and save them in a file to be read below
			Printf("QPYTHIA: make the projections and save them in a file to be read below");
			
			TString pathOrig = 
			"/data/Work/jets/JetMass/PbPbJEWEL/Marco/qpythia_v11_jet_shapes_50_merged.root";
			
			TFile *finorig = new TFile(pathOrig);
			if(!finorig->IsOpen()){
				Printf("%s not found, cannot proceed.", pathOrig.Data());
				return 0x0;
			}
			TH3F  *h3d     = (TH3F*)finorig->Get("hJetPtEtaMass_R4");
			if(!h3d){
				Printf("3D histo from Marco not fouund in %s", pathOrig.Data());
				finorig->ls();
			}
			
			TFile* fout = new TFile(pathModels, "recreate");
			for(Int_t ipt = 0; ipt < nptbins; ipt++){
				TH1 *hproj = Marco_get_shape_hist(h3d, ptlims[ipt], ptlims[ipt+1], TString::Format("hM_QPYTHIA_%d", ipt), 8); // same name as below for mod == 3
				fout->cd();
				hproj->Write();
			}
			fout->Close();
			delete fout;
		} else delete ftest; //remove the test file
		
	}
	
	TFile *fModels = new TFile(pathModels);
	if(!fModels->IsOpen()){
		Printf("File %s not found", pathModels.Data());
		return 0;
	}

	TString nameMSt = "";
	Int_t marker = 28;
	Int_t color  = kOrange+1;
	Int_t linest = 1;
	Int_t fillst = 1;
	if(mod == 1) nameMSt = "hM_JEWEL RecOn_";
	if(mod == 2) {
		nameMSt = "hM_JEWEL RecOff_";
		marker = 28;
		color  = kRed;
		fillst = 3004;
		linest = 3;
	}
	if(mod == 3) {
		nameMSt = "hM_QPYTHIA_";
		marker = 33;
		color  = kRed+2;
		fillst = 3006;
		linest = 7;
	}
	if(mod == 4) {
		nameMSt = "hM_Vacuum_";
		marker = 21;
		color  = kViolet+5;
		fillst = 0;
		linest = 5;
	}
	
	const Int_t nhM = 5;
	Int_t offset = 0;
	
	TList *listModels = new TList();
	listModels->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fModels->Get(TString::Format("%s%d", nameMSt.Data(), ih+offset));
		if(!hMassSt) continue; // this is important because the new QYTHIA file is built to have less bins, only 3
		hMassSt->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		hMassSt->SetMarkerStyle(marker);
		hMassSt->SetMarkerColor(color);
		hMassSt->SetLineColor(hMassSt->GetMarkerColor());
		if(useline) {
			hMassSt->SetLineStyle(linest);
			hMassSt->SetFillColor(hMassSt->GetLineColor());
			hMassSt->SetMarkerStyle(1);
			hMassSt->SetFillStyle(fillst);
			hMassSt->SetLineWidth(3);
		}
		hMassSt->SetMarkerSize(2);
		listModels->Add(hMassSt);
	}
	
	return listModels;

}
//________________________________________________________________________
TList* GetHERWIG(Int_t system, Bool_t useline){

	TString pathMartapPb = "/data/Work/jets/JetMass/pPbHERWIG/ChJetMassGenerators.root";
	TString pathHWPbPb = "";
	
	TString fileIn = pathMartapPb;
	if(system == 1) fileIn = pathHWPbPb;
	
	TFile *fRespPbM = new TFile(fileIn);
	if(!fRespPbM->IsOpen()){
		Printf("File %s not found", fileIn.Data());
		return 0;
	}

	TString nameMSt = "hM_Herwig_Pt";
	const Int_t nhM = 3;
	
	TList *listpPb = new TList();
	listpPb->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fRespPbM->Get(TString::Format("%s%.0fTo%.0f", nameMSt.Data(), ptlims[ih], ptlims[ih+1]));
		
		// set this nicely at some point
		hMassSt->SetLineColor(kBlue+1);
		hMassSt->SetMarkerColor(hMassSt->GetLineColor());
		hMassSt->SetMarkerSize(2);
		hMassSt->SetMarkerStyle(24);
		if(useline) {
			hMassSt->SetLineStyle(2);
			hMassSt->SetLineWidth(3);
			hMassSt->SetFillColor(hMassSt->GetLineColor());
			hMassSt->SetMarkerStyle(1);
			hMassSt->SetFillStyle(1001);//3004
		}
		
		hMassSt->GetXaxis()->SetRangeUser(0., maxM);
		
		listpPb->Add(hMassSt);
	}
	
	return listpPb;

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
void RatiopPbPbPb(Bool_t useline = kTRUE, Bool_t stylezero = kTRUE, Bool_t drawRaghav = kTRUE, Bool_t pPbpaperprop = kFALSE, Bool_t corrkine = kFALSE, Bool_t show1134 = kFALSE){
	TString suff = "";
	if(!corrkine) suff = "NoKineCor";
	//TGaxis::SetMaxDigits(2);
	if(stylezero) {
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		gStyle->SetTextFont(42);
		
		gStyle->SetLabelSize(0.06, "X");
		gStyle->SetTitleOffset(1.12, "X");
		gStyle->SetTitleSize(0.055, "X");
		
		gStyle->SetNdivisions(705, "X");
		
		gStyle->SetLabelSize(0.06, "Y");
		gStyle->SetTitleOffset(1.34, "Y");
		gStyle->SetTitleSize(0.048, "Y");
		
		gStyle->SetNdivisions(1004, "Y");
		
		
		
	}
	
	gStyle->SetPadBottomMargin(.13);
	gStyle->SetPadLeftMargin(.16);
	gStyle->SetPadRightMargin(.06);
	
	TLatex latext;
	latext.SetTextSize(0.05);
	
	
	// define how to draw
	TString drawsimopt = "P", drawsimoptsame = "Psames";
	TString optlegsim  = "LP";
	if(useline) {
		drawsimopt = "LP";
		drawsimoptsame = "LPsames";
		optlegsim = "L";
	}
	//get the PbPb results
	TList* listPbPb = GetPbPbResultsGraph(); //GetPbPbResults();
		
	//get the pPb results
	TList* listpPb = GetpPbResults(corrkine);
	TList* listpPbM = GetpPbMarta();
	
	//get Pythia pPb
	TList* listPythia502 = GetPythiapPb();
	//get Pythia pPb
	TList* listPythia502M = GetPythiapPbMarta(useline);
	
	//get Pythia PbPb
	TList* listPythia276 = GetPythiaPbPb(useline);
	
	//get Jewel PbPb
	TList* listJewelPbPb = GetJEWELPbPb(useline);
	
	TList* listJewelMRecoilOff = //GetJEWELPbPbrecoilOff(); 
	GetJEWELQPYTHIA(2, useline);
	
	//get Jewel pPb
	TList* listJewelpp = GetJEWELpp(useline);
	
	//get Q-PYTHIA
	TList* listqpythia = GetJEWELQPYTHIA(3, useline);
	Printf("List QPYTHIA %p", listqpythia);
	//get HERWIG
	TList* listHERWIGpPb = GetHERWIG(0, useline);
	
	
	if(!listPbPb){ 
		Printf("List PbPb not found");
		return;
	}
	if(!listpPb){ 
		Printf("List pPb not found");
		return;
	}
	
	if(!listPythia502){
		Printf("List Pythia for pPb energy not found");
		return;
	}
	
	if(!listPythia276){
		Printf("List Pythia for PbPb energy not found");
		return;
	}
	
	if(!listJewelPbPb){ 
		Printf("List Jewel PbPb not found");
		return;
	}
	
	if(!listJewelpp){ 
		Printf("List Jewel pp not found");
		return;
	}
	
	if(!listJewelMRecoilOff){ 
		Printf("List Jewel recoil off not found");
		return;
	}
	
	if(!listqpythia){
		Printf("List qpythia not found");
		return;
	}
	
	if(!listHERWIGpPb){
		Printf("List herwig not found");
		return;
	}
	
	Int_t npPb    = listpPb ->GetEntries()/2;
	Int_t nPbPb   = listPbPb->GetEntries()/2;
	Int_t npPbM   = listpPbM->GetEntries()/2;
	Int_t nPy502  = listPythia502->GetEntries();
	Int_t nPy502M = listPythia502M->GetEntries();
	Int_t nPy276  = listPythia276->GetEntries();
	Int_t nJewpp  = listJewelpp ->GetEntries()/2;
	Int_t nJewPbPb= listJewelPbPb->GetEntries();
	Int_t nJewPbPbMRoff = listJewelMRecoilOff->GetEntries();
	Int_t nqpythia= listqpythia->GetEntries();
	Int_t nhwppb  = listHERWIGpPb->GetEntries();
	
	Int_t offsetPbPb   = nPbPb   - npPb;
	Int_t offsetpPbM   = npPbM   - npPb;
	Int_t offsetPy502  = nPy502  - npPb;
	Int_t offsetPy502M = nPy502M - npPb;
	Int_t offsetPy276  = nPy276  - npPb;
	Int_t offsetJewpp  = nJewpp  - npPb;
	Int_t offsetJewPbPb= nJewPbPb- npPb;
	Int_t offsetJewPbPbMRoff = nJewPbPbMRoff- npPb;
	Int_t offsetqpythia= nqpythia- npPb;
	Int_t offsethwppb  = nhwppb- npPb;
	
	Printf("offsets %d %d %d %d %d %d %d %d %d ", offsetPbPb, offsetpPbM , offsetPy502, offsetPy502M, offsetPy276, offsetJewpp, offsetJewPbPb, offsetJewPbPbMRoff, offsetqpythia);
	
	const Int_t nhM = npPb;
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nhM, nx, ny, dx, dy, 1, 400);
	TCanvas *cMass = new TCanvas(TString::Format("cMass%s", suff.Data()), "Mass histograms", dx, dy);
	cMass->Divide(nx, ny);
	TLegend *legMass = new TLegend(0.2, 0.55, 0.9, 0.75);
	legMass->SetBorderSize(0);
	legMass->SetFillStyle(0);
	
	TCanvas *cRatioPbPbOpPb = new TCanvas(TString::Format("cRatioPbPbOpPb%s", suff.Data()), "Ratio PbPb/pPb", dx, dy);
	cRatioPbPbOpPb->Divide(nx, ny);
	TLegend *legRatioPbPbOpPb = new TLegend(0.2, 0.55, 0.9, 0.75);
	legRatioPbPbOpPb->SetBorderSize(0);
	legRatioPbPbOpPb->SetFillStyle(0);
	
	TCanvas *cMasspPb = new TCanvas(TString::Format("cMasspPb%s", suff.Data()), "Mass pPb compared with Pythia", dx, dy);
	cMasspPb->Divide(nx, ny);
	TLegend *legMasspPb = new TLegend(0.2, 0.5, 0.92, 0.7);
	legMasspPb->SetBorderSize(0);
	legMasspPb->SetFillStyle(0);
	
	TCanvas *cRelUnc = new TCanvas(TString::Format("cRelUncSys%s", suff.Data()), "Relative systematic uncertainties mass spectra", dx, dy);
	cRelUnc->Divide(nx, ny);
	TLegend *legRelUnc = new TLegend(0.4, 0.5, 0.8, 0.9);
	legRelUnc->SetBorderSize(0);
	legRelUnc->SetFillStyle(0);
	
	TCanvas *cRatiopPbOPy = new TCanvas(TString::Format("cRatiopPbOPy%s", suff.Data()), "Ratio pPb over Pythia 5.02 TeV", dx, dy);
	cRatiopPbOPy->Divide(nx, ny);
	TLegend *legRatiopPbOPy = new TLegend(0.22, 0.55, 0.9, 0.75);
	legRatiopPbOPy->SetBorderSize(0);
	legRatiopPbOPy->SetFillStyle(0);
	
	TCanvas *cMassPbPb = new TCanvas(TString::Format("cMassPbPb%s", suff.Data()), "Mass PbPb compared with Pythia", dx, dy);
	cMassPbPb->Divide(nx, ny);
	TLegend *legMassPbPb1 = new TLegend(0.2, 0.55, 0.9, 0.75);
	legMassPbPb1->SetBorderSize(0);
	legMassPbPb1->SetFillStyle(0);
	TLegend *legMassPbPb2 = new TLegend(0.2, 0.53, 0.8, 0.68);
	legMassPbPb2->SetBorderSize(0);
	legMassPbPb2->SetFillStyle(0);
	
	TLegend* legdummy = new TLegend(0.2, 0.6, 0.8, 0.75);
	legdummy->SetBorderSize(0);
	legdummy->SetFillStyle(0);
	
	TCanvas *cMassPbPbOnly = new TCanvas(TString::Format("cMassPbPbOnly%s", suff.Data()), "Mass PbPb", dx, dy);
	cMassPbPbOnly->Divide(nx, ny);
	TLegend *legMassPbPbOnly = new TLegend(0.15, 0.55, 0.9, 0.75);
	legMassPbPbOnly->SetBorderSize(0);
	legMassPbPbOnly->SetFillStyle(0);
	
	TCanvas *cMassPbPbPy = new TCanvas(TString::Format("cMassPbPbPy%s", suff.Data()), "Mass PbPb and Pythia", dx, dy);
	cMassPbPbPy->Divide(nx, ny);
	TLegend *legMassPbPbPy = new TLegend(0.15, 0.55, 0.9, 0.75);
	legMassPbPbPy->SetBorderSize(0);
	legMassPbPbPy->SetFillStyle(0);
	
	TCanvas *cMassPbPbModels = new TCanvas(TString::Format("cMassPbPbModels%s", suff.Data()), "Mass PbPb in models", dx, dy);
	cMassPbPbModels->Divide(nx, ny);
	TLegend *legMassPbPbModels = new TLegend(0.15, 0.55, 0.9, 0.75);
	legMassPbPbModels->SetBorderSize(0);
	legMassPbPbModels->SetFillStyle(0);
	
	TCanvas *cRatioPbPbOPy = new TCanvas(TString::Format("cRatioPbPbOPy%s", suff.Data()), "Ratio PbPb over Pythia 2.76 TeV", dx, dy);
	cRatioPbPbOPy->Divide(nx, ny);
	TLegend *legRatioPbPbOPy = new TLegend(0.22, 0.63, 0.9, 0.82);
	legRatioPbPbOPy->SetBorderSize(0);
	legRatioPbPbOPy->SetFillStyle(0);
	
	TCanvas *cRatioDataOPy = new TCanvas(TString::Format("cRatioDataOPy%s", suff.Data()), "Ratio data over Pythia", dx, dy);
	cRatioDataOPy->Divide(nx, ny);
	TLegend *legRatioDataOPy = new TLegend(0.22, 0.58, 0.9, 0.78);
	legRatioDataOPy->SetBorderSize(0);
	legRatioDataOPy->SetFillStyle(0);
	
	Int_t n = 2;
	
	//TPaveText **pvpt = new TPaveText*[nhM];
	
	TFile *fOutput = new TFile("FinalResults.root", "recreate");
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		maxM = maxRangeMassFinal[ih];
		//pvpt[ih] = new TPaveText(0.33, 0.8, 0.83, 0.9, "NDC");
		//pvpt[ih]->SetFillStyle(0);
		//pvpt[ih]->SetBorderSize(0);
		//pvpt[ih]->AddText(TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		
		//mass
		TH1D* hMpPb     = (TH1D*)listpPb ->At(n*ih);
		TH1D* hMPbPb    = (TH1D*)listPbPb->At(n*(ih+offsetPbPb));
		TH1D* hMpPbM    = (TH1D*)listpPbM->At(n*(ih+offsetpPbM));
		                
		TH1D* hMPy502   = (TH1D*)listPythia502 ->At(ih+offsetPy502 );
		TH1D* hMPy502M  = (TH1D*)listPythia502M->At(ih+offsetPy502M);
		TH1D* hMPy276   = (TH1D*)listPythia276 ->At(ih+offsetPy276 );
		TH1D* hMJewpp   = (TH1D*)listJewelpp  ->At((ih+offsetJewpp));
		TH1D* hMJewPbPb = (TH1D*)listJewelPbPb->At((ih+offsetJewPbPb));
		
		TH1D* hMJewPbPbMRoff= (TH1D*)listJewelMRecoilOff ->At(ih+offsetJewPbPbMRoff);
		TH1D* hMQPy     = (TH1D*)listqpythia->At(ih+offsetqpythia);
		TH1D* hMHwpPb   = (TH1D*)listHERWIGpPb->At(ih+offsethwppb);
		
		hMPy502->Scale(1./hMPy502->Integral("width"));
		hMPy502M->Scale(1./hMPy502M->Integral("width"));
		
		hMPy276->Scale(1./hMPy276->Integral("width"));
		
		//syst
		TH1D* hSypPb  = (TH1D*)listpPb ->At(n*ih+1);
		TH1D* hSypPbM = (TH1D*)listpPbM->At(n*(ih+offsetpPbM)+1);
		TH1D* hSyPbPb = (TH1D*)listPbPb->At(n*(ih+offsetPbPb)+1);
		//ratio Jewel
		//TH1D* hPbPboppJewel = (TH1D*)listJewelPbPb->At(n*(ih+offsetJewPbPb)+1);
		
		//transform all models in TGraph (needed for unclear drawing reasons)
		
		TGraph* gMPy502   =      new TGraph(hMPy502);
		TGraph* gMPy502M  =      new TGraph(hMPy502M);
		TGraph* gMPy276   =      new TGraph(hMPy276);
		TGraph* gMJewpp   =      new TGraph(hMJewpp);
		TGraph* gMJewPbPb =      new TGraph(hMJewPbPb);
		TGraph* gMJewPbPbMRoff = new TGraph(hMJewPbPbMRoff);
		TGraph* gMQPy     =      new TGraph(hMQPy);
		TGraph* gMHwpPb   =      new TGraph(hMHwpPb);
		
		
		//Draw comparison pPb, PbPb
		
		cMass->cd(ih+1);
		gPad->SetTicks(1,1);
		Printf("Drawing %s, %s, %s, %s, %s", hSyPbPb->GetName(), hMPbPb->GetName(), hSypPb->GetName(), hSypPbM->GetName(), hMpPbM->GetName());
		hSyPbPb->GetYaxis()->SetRangeUser(0., 0.25);
		
		hSyPbPb->GetXaxis()->SetRangeUser(0., maxM);
		hSyPbPb->Draw("E2");
		hMPbPb ->Draw("sames");
		hSypPb ->Draw("E2sames");
		hMpPb  ->Draw("sames");
		if(pPbpaperprop){
			hSypPbM->Draw("E2sames");
			hMpPbM ->Draw("sames");
		}
		if(ih == 1) {
			latext.DrawLatex(4, 0.18, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
			
		}
		latext.DrawLatex(4, 0.22, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih]->Draw();
		if(ih == nhM-1){
			legMass->AddEntry(hSyPbPb , "Pb-Pb 0-10% #sqrt{#it{s}}_{NN} = 2.76 TeV");
			//legMass->AddEntry(hSyPbPb, "Systematic Pb-Pb" , "F");
			legMass->AddEntry(hSypPb ,  "pPb #sqrt{#it{s}}_{NN} = 5.02 TeV");
			//legMass->AddEntry(hSypPb,  "Systematic p-Pb" , "F");
			if(pPbpaperprop){
				legMass->AddEntry(hMpPbM , "Mass p-Pb paper prop", "PL");
				legMass->AddEntry(hSypPbM,  "Sys p-Pb paper prop" , "F");
			}
			legMass->Draw();
		}
		
		//relative uncertainties
		TH1D* hRelSysUncpPb  = (TH1D*)hMpPb->Clone(TString::Format("hRelSysUncMpPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRelSysUncpPb->SetLineColor(kPink-9);
		hRelSysUncpPb->SetLineWidth(2);
		
		TH1D* hRelSysUncpPbM = (TH1D*)hMpPbM->Clone(TString::Format("hRelSysUncMpPbM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		TH1D* hRelSysUncPbPb = (TH1D*)hMPbPb->Clone(TString::Format("hRelSysUncMPbPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRelSysUncPbPb->SetTitle("Relative systematic Uncertainties; #it{M} (GeV/#it{c}^2); Syst Unc / #it{M}");
		for(Int_t ib = 0; ib<hRelSysUncpPb->GetNbinsX(); ib++){
			hRelSysUncpPb->SetBinContent(ib+1, hSypPb->GetBinError(ib+1)/hMpPb->GetBinContent(ib+1));
			hRelSysUncpPb->SetBinError(ib+1, 0 );
		}
		for(Int_t ib = 0; ib<hRelSysUncpPbM->GetNbinsX(); ib++){
			hRelSysUncpPbM->SetBinContent(ib+1, hSypPbM->GetBinError(ib+1)/hMpPbM->GetBinContent(ib+1));
			hRelSysUncpPbM->SetBinError(ib+1, 0 );
		
		}
		for(Int_t ib = 0; ib<hRelSysUncpPb->GetNbinsX(); ib++){
			hRelSysUncPbPb->SetBinContent(ib+1, hSyPbPb->GetBinError(ib+1)/hMPbPb->GetBinContent(ib+1));
			hRelSysUncPbPb->SetBinError(ib+1, 0 );
		
		}
		
		cRelUnc->cd(ih+1);
		gPad->SetTicks(1,1);
		hRelSysUncPbPb ->GetYaxis()->SetRangeUser(0., 1.);
		hRelSysUncPbPb ->Draw("hist");
		hRelSysUncpPb  ->Draw("histsames");
		hRelSysUncpPbM ->Draw("histsames");
		if(ih == nhM-1){
			legRelUnc->AddEntry(hRelSysUncpPb, "p-Pb", "L");
			legRelUnc->AddEntry(hRelSysUncpPbM, "p-Pb paper prop", "L");
			legRelUnc->AddEntry(hRelSysUncPbPb, "Pb-Pb", "L");
			legRelUnc->Draw();
		}
		
		//Draw comparison pPb Pythia 5.02 TeV
		cMasspPb->cd(ih+1);
		
		gPad->SetTicks(1,1);
		hSypPb ->GetYaxis()->SetRangeUser(0., 0.25);
		hSypPb->SetMarkerStyle(hMpPb->GetMarkerStyle());
		hSypPb->SetMarkerColor(hMpPb->GetMarkerColor());
		hSypPb ->Draw("E2");
		hMpPb  ->Draw("sames");
		if(show1134) hMPy502->Draw(drawsimopt);
		gMPy502M->Draw(drawsimoptsame);
		gMHwpPb->Draw(drawsimoptsame);
		//hMJewpp-> Draw("sames");
		if(ih == nhM-1){
			legMasspPb->AddEntry(hSypPb, "p-Pb #sqrt{s_{NN}} = 5.02 TeV");
			//legMasspPb->AddEntry(hSypPb, "Systematic p-Pb", "F");
			if(show1134) legMasspPb->AddEntry(hMPy502, "PYTHIA Perugia 2011 #sqrt{s_{NN}} = 5.02 TeV", optlegsim);
			legMasspPb->AddEntry(hMPy502M, "PYTHIA Perugia 2011 #sqrt{s_{NN}} = 5.02 TeV", optlegsim); //  #sqrt{s_{NN}} = 5.02 TeV // paper prop
			//legMasspPb->AddEntry(hMJewpp, "JEWEL+PYTHIA pp", "LP");
			legMasspPb->AddEntry(hMHwpPb, "HERWIG EE5C #sqrt{s_{NN}} = 5.02 TeV", optlegsim);
			legMasspPb->Draw();
		}
		if(ih == 1) {
			latext.DrawLatex(4, 0.18, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
		}
		latext.DrawLatex(4, 0.22, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih]->Draw();
		
		
		cMassPbPbOnly->cd(ih+1);
		gPad->SetTicks(1,1);
		hSyPbPb  ->Draw("E2");
		hMPbPb   ->Draw("sames");
		if(ih == 1) {
			latext.DrawLatex(4, 0.18, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
		}
		latext.DrawLatex(4, 0.22, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih]->Draw();
		
		cMassPbPbPy->cd(ih+1);
		
		gPad->SetTicks(1,1);
		hSyPbPb  ->Draw("E2");
		hMPbPb   ->Draw("sames");
		gMPy276  ->Draw(drawsimoptsame);
		if(ih == 1) {
			latext.DrawLatex(4, 0.18, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
		}
		latext.DrawLatex(4, 0.22, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih] ->Draw();
		//Draw comparison PbPb Pythia 2.76 TeV
		cMassPbPb->cd(ih+1);
		
		gPad->SetTicks(1,1);
		hSyPbPb->SetMarkerStyle(hMPbPb->GetMarkerStyle());
		hSyPbPb->SetMarkerColor(hMPbPb->GetMarkerColor());
		hSyPbPb  ->Draw("E2");
		hMPbPb   ->Draw("sames");
		gMPy276  ->Draw(drawsimoptsame);
		//temporary don't draw this because wrong curve from Raghav
		gMJewPbPb->GetXaxis()->SetRangeUser(0., 25.);
		if(drawRaghav) gMJewPbPb->Draw(drawsimoptsame);
		gMJewPbPbMRoff->Draw(drawsimoptsame);
		gMQPy    ->Draw(drawsimoptsame);
		latext.DrawLatex(4, 0.22, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih] ->Draw();
		
		cMassPbPbModels->cd(ih+1);
		
		gPad->SetTicks(1,1);
		gMJewPbPbMRoff  ->GetYaxis()->SetNdivisions(1004, "Y");
		gMJewPbPbMRoff  ->GetYaxis()->SetRangeUser(0, 0.25);
		gMJewPbPbMRoff  ->GetXaxis()->SetRangeUser(0, 25);
		gMJewPbPbMRoff->Draw(TString::Format("A%s", drawsimopt.Data()));
		gMPy276  ->Draw(drawsimoptsame);
		if(drawRaghav) gMJewPbPb->Draw(drawsimoptsame);
		
		gMQPy    ->Draw(drawsimoptsame);
		latext.DrawLatex(4, 0.22, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih] ->Draw();
		
		if(ih == nhM-1){
			legMassPbPbOnly->AddEntry(hSyPbPb,    "Pb-Pb 0-10% #sqrt{s_{NN}} = 2.76 TeV");
			//legMassPbPbOnly->AddEntry(hSyPbPb,   "Systematic Pb-Pb" , "F");
			cMassPbPbOnly->cd(ih+1); legMassPbPbOnly->Draw();
			
			
			legMassPbPbModels->AddEntry(hMPy276,   "PYTHIA  Perugia 2011 #sqrt{s_{NN}} = 2.76 TeV", optlegsim);
			if(drawRaghav) legMassPbPbModels->AddEntry(hMJewPbPb, "JEWEL+PYTHIA Recoil on 0-10% PbPb", optlegsim);
			legMassPbPbModels->AddEntry(hMJewPbPbMRoff, "JEWEL+PYTHIA Recoil off 0-10% PbPb", "LP");
			legMassPbPbModels->AddEntry(hMQPy,      "Q-PYTHIA", optlegsim);
			cMassPbPbModels->cd(ih+1); legMassPbPbModels->Draw();
			
			legMassPbPbPy->AddEntry(hSyPbPb,    "Pb-Pb 0-10% #sqrt{s_{NN}} = 2.76 TeV");
			//legMassPbPbPy->AddEntry(hSyPbPb,   "Systematic Pb-Pb" , "F");
			legMassPbPbPy->AddEntry(hMPy276,   "PYTHIA  Perugia 2011", optlegsim);
			cMassPbPbPy->cd(ih+1); legMassPbPbPy->Draw();
			
			
			if(drawRaghav) legMassPbPb2->AddEntry(hMJewPbPb, "Recoil on", optlegsim);
			legMassPbPb2->AddEntry(hMJewPbPbMRoff, "Recoil off", optlegsim);
			//latext.DrawLatex(6, 0.19, "Recoil on");
			//latext.DrawLatex(6, 0.17, "Recoil off");
			cMassPbPb->cd(ih+1); latext.DrawLatex(2, 0.18, "JEWEL + PYTHIA 0-10% Pb-Pb"); legMassPbPb2->Draw();
		}
		//if(ih == 0){
		//	legdummy->AddEntry(hMJewPbPb, "Recoil on", "P");
		//	cMassPbPb->cd(ih+1); 
		//	
		//	legdummy->Draw();
		//}
		if(ih == nhM-2){
			legMassPbPb1->AddEntry(hSyPbPb,    "Pb-Pb 0-10% #sqrt{s_{NN}} = 2.76 TeV");
			//legMassPbPb->AddEntry(hSyPbPb,   "Systematic Pb-Pb" , "F");
			legMassPbPb1->AddEntry(hMPy276,   "PYTHIA  Perugia 2011", optlegsim); //#sqrt{s_{NN}} = 2.76 TeV
			legMassPbPb1->AddEntry(hMQPy,      "Q-PYTHIA", optlegsim);
			cMassPbPb->cd(ih+1); legMassPbPb1->Draw();
		}
		//ratio PbPb/pPb
		TH1* hRatioPbPbOpPb;
		TH1* hDividepPb;
		UniformTH1FForDivide(hMPbPb, hMpPb, hRatioPbPbOpPb, hDividepPb, "TH1D"); hRatioPbPbOpPb->SetName(TString::Format("hRatioPbPbOpPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatioPbPbOpPb->GetYaxis()->SetTitleOffset(1.15);
		
		hRatioPbPbOpPb->SetMarkerColor(kMagenta+2);
		hRatioPbPbOpPb->SetLineColor(hRatioPbPbOpPb->GetMarkerColor());
		hRatioPbPbOpPb->Divide(hDividepPb);
		hRatioPbPbOpPb->GetXaxis()->SetRangeUser(0., maxRangeMassFinal[ih]);
		
		//line for reference 
		TLine *lineOne = new TLine(0, 1., maxRangeMassFinal[ih], 1.);
		lineOne->SetLineStyle(2);
		lineOne->SetLineWidth(2);
		lineOne->SetLineColor(kGray);
		
		//syst ratio PbPb/pPb
		
		// - given by ROOT
		TH1* hRatioPbPbOpPbSys;
		TH1* hRDivpPb;
		UniformTH1FForDivide(hSyPbPb, hSypPb, hRatioPbPbOpPbSys, hRDivpPb, "TH1D");
		hRatioPbPbOpPbSys->SetName(TString::Format("hRatioPbPbOpPbSys_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatioPbPbOpPbSys->GetYaxis()->SetTitle("#Rgothic_{#sqrt{#it{s}}}");
		hRatioPbPbOpPbSys->GetYaxis()->SetTitleOffset(1.15);
		//hRatioPbPbOpPbSys->GetYaxis()->SetTitleOffset(1.7);
		hRatioPbPbOpPbSys->Divide(hRDivpPb);
		hRatioPbPbOpPbSys->SetFillStyle(1001);
		hRatioPbPbOpPbSys->SetFillColor(kMagenta-10);
		hRatioPbPbOpPbSys->SetMarkerStyle(hRatioPbPbOpPb->GetMarkerStyle());
		hRatioPbPbOpPbSys->SetMarkerColor(hRatioPbPbOpPb->GetMarkerColor());
		hRatioPbPbOpPbSys->SetLineWidth(2);
		hRatioPbPbOpPbSys->SetLineColor(kBlack);
		
				
		// - Calculated as R x sqrt(relerrPbPb^2 + relerrpPb^2)
		// gives the same result but it's more messy...comment
		/*
		TH1* hRatioPbPbOpPbSysC;
		TH1* hRDivpPbC;
		Printf("N bins %d and %d", hSyPbPb->GetNbinsX(), hSypPb->GetNbinsX());
		UniformTH1FForDivide(hSyPbPb, hSypPb, hRatioPbPbOpPbSysC, hRDivpPbC, "TH1D");
		hRatioPbPbOpPbSysC->SetName(TString::Format("hRatioPbPbOpPbSysC_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		Printf("BNow N bins %d and %d", hRatioPbPbOpPbSysC->GetNbinsX(), hRDivpPbC->GetNbinsX());
		hRatioPbPbOpPbSysC->SetFillStyle(0);
		hRatioPbPbOpPbSysC->SetLineStyle(2);
		hRatioPbPbOpPbSysC->SetLineWidth(2);
		hRatioPbPbOpPbSysC->SetLineColor(kGreen+3);
		
		CalculateSysRatio((TH1D*)hRatioPbPbOpPbSysC, (TH1D*)hRDivpPbC);
		for(Int_t ib = 0; ib<hRatioPbPbOpPbSysC->GetNbinsX(); ib++) Printf(" and =  %f +- %f", hRatioPbPbOpPbSysC->GetBinContent(ib+1), hRatioPbPbOpPbSysC->GetBinError(ib+1));
		cRatioPbPbOpPb->cd(ih+1);
		hRatioPbPbOpPbSysC->Draw("E2sames");
		*/
		
		// ratio PbPb/pPbMarta
		
		TH1* hRatioPbPbOpPbM;
		TH1* hDividepPbM;
		UniformTH1FForDivide(hMPbPb, hMpPbM, hRatioPbPbOpPbM, hDividepPbM, "TH1D"); 
		hRatioPbPbOpPbM->SetName(TString::Format("hRatioPbPbOpPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		hRatioPbPbOpPbM->GetYaxis()->SetTitleOffset(1.15);
		hRatioPbPbOpPbM->SetMarkerColor(kBlack);
		hRatioPbPbOpPbM->SetLineColor(kBlack);
		hRatioPbPbOpPbM->Divide(hDividepPbM);
		hRatioPbPbOpPbM->GetXaxis()->SetRangeUser(0., maxRangeMassFinal[ih]);
		
		//syst ratio PbPb/pPbM
		
		// - given by ROOT
		TH1* hRatioPbPbOpPbSysM;
		TH1* hRDivpPbSyM;
		UniformTH1FForDivide(hSyPbPb, hSypPbM, hRatioPbPbOpPbSysM, hRDivpPbSyM, "TH1D");
		hRatioPbPbOpPbSysM->SetName(TString::Format("hRatioPbPbOpPbSysM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRDivpPbSyM->SetName(TString::Format("hRDivpPbSyM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		hRatioPbPbOpPbSysM->GetYaxis()->SetTitleOffset(1.15);
		//hRatioPbPbOpPbSysM->GetYaxis()->SetTitle("Ratio");
		hRatioPbPbOpPbSysM->Divide(hRDivpPbSyM);
		hRatioPbPbOpPbSysM->SetFillStyle(1001);
		hRatioPbPbOpPbSysM->SetFillColor(kGray);
		hRatioPbPbOpPbSysM->SetMarkerStyle(24);
		hRatioPbPbOpPbSysM->SetLineColor(kBlack);
		hRatioPbPbOpPbSysM->SetMarkerColor(kBlack);
		
		// ratio Pythia 276/ Pythia 502
		TH1* hRatio276O502;
		TH1* hDividePyPy502;
		UniformTH1FForDivide(hMPy276, hMPy502, hRatio276O502, hDividePyPy502, "TH1D");
		hRatio276O502->Divide(hDividePyPy502);
		
		// ratio Pythia 276/ Pythia 502 (from Marta)
		TH1* hRatio276O502M;
		TH1* hDividePyPy502M;
		UniformTH1FForDivide(hMPy276, hMPy502M, hRatio276O502M, hDividePyPy502M, "TH1D");
		hRatio276O502M->SetName(TString::Format("hRatio276O502M_Pt%d", ih));
		hRatio276O502M->GetYaxis()->SetTitleOffset(1.15);
		hRatio276O502M->Divide(hDividePyPy502M);
		hRatio276O502M->SetLineColor(kRed);
		hRatio276O502M->SetMarkerColor(kRed);
		hRatio276O502M->SetMarkerStyle(20);
		if(useline){
			hRatio276O502M->SetMarkerStyle(1);
			hRatio276O502M->SetFillStyle(1001);
			hRatio276O502M->SetFillColor(hRatio276O502M->GetLineColor());
		}
		
		//for drawing
		TGraph *gRatio276O502M = new TGraph(hRatio276O502M);
		
		cRatioPbPbOpPb->cd(ih+1);
		
		gPad->SetTicks(1,1);
		
		if(pPbpaperprop){
			//Old result paper draft
			hRatioPbPbOpPbSysM->GetYaxis()->SetRangeUser(0., 4);
			hRatioPbPbOpPbSysM->GetXaxis()->SetRangeUser(0., maxM);
			hRatioPbPbOpPbSysM->Draw("E2");
			hRatioPbPbOpPbM->Draw("sames");
			hRatioPbPbOpPbSys->Draw("E2sames");
		}
		//new result
		if(!pPbpaperprop) {
			hRatioPbPbOpPbSys->GetYaxis()->SetRangeUser(0., 4);
			hRatioPbPbOpPbSys->GetXaxis()->SetRangeUser(0., maxM);
			hRatioPbPbOpPbSys->Draw("E2");
		}
		hRatioPbPbOpPb->Draw("sames");
		// pythia ratio from paper draft
		hRatio276O502M->Draw("E3sames");//drawsimoptsame
		// pythia ratio using 1134 output for pPb
		if(show1134) hRatio276O502->Draw(drawsimoptsame);
		//hPbPboppJewel->Draw("Lsames");
		//reference line
		lineOne->Draw();
		if(ih == 1) {
			latext.DrawLatex(4, 3, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 2.5, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
		}
		latext.DrawLatex(4, 3.6, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih]->Draw();
		
		if(ih == nhM-1){
			legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPbSys, "Data Pb-Pb / p-Pb");// #sqrt{s_{NN}}
			//legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPbSys, "Sys Pb-Pb/p-Pb", "F");
			if(pPbpaperprop){
				legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPbM, "Ratio PbPb/pPb paper prop", "LP");
				legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPbSysM, "Sys PbPb/pPb paper prop", optlegsim);
			}
			if(show1134) legRatioPbPbOpPb->AddEntry(hRatio276O502, "PYTHIA(2.76)/PYTHIA(5.02)Matched", optlegsim);
			
			legRatioPbPbOpPb->AddEntry(hRatio276O502M, "PYTHIA 2.76TeV / 5.02TeV", optlegsim); // paper prop
			//legRatioPbPbOpPb->AddEntry(hPbPboppJewel, "Ratio JEWEL", "LP");
			legRatioPbPbOpPb->Draw();
			
		}
		
		//ratio pPb/PYTHIA
		TH1* hRatiopPbOPyM;
		TH1* hDividePy502M;
		UniformTH1FForDivide(hMpPb, hMPy502M, hRatiopPbOPyM, hDividePy502M, "TH1D"); hRatiopPbOPyM->SetName(TString::Format("hRatiopPbOPyM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOPyM->GetYaxis()->SetTitle("Data/PYTHIA");
		
		hRatiopPbOPyM->GetYaxis()->SetTitleOffset(1.15);
		hRatiopPbOPyM->SetMarkerColor(hMpPb->GetMarkerColor());
		hRatiopPbOPyM->SetLineColor(hMpPb->GetLineColor());
		            
		hRatiopPbOPyM->Divide(hDividePy502M);
		hRatiopPbOPyM->GetXaxis()->SetRangeUser(0., maxRangeMassFinal[ih]);
		hRatiopPbOPyM->GetYaxis()->SetRangeUser(0., 4);
		
		// sys ratio pPb/PYTHIA
		TH1* hRatiopPbOPy502SysM;
		TH1* hDivideSysPy502M;
		UniformTH1FForDivide(hSypPb, hMPy502M, hRatiopPbOPy502SysM, hDivideSysPy502M, "TH1D");
		TH1* hDividePy502MNoErr = (TH1*)hDivideSysPy502M->Clone(TString::Format("hDividePy502MNoErr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		for(Int_t ib = 0; ib<hDividePy502MNoErr->GetNbinsX(); ib++) hDividePy502MNoErr->SetBinError(ib, 0);
		
		hRatiopPbOPy502SysM->SetName(TString::Format("hRatiopPbOPy502SysM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOPy502SysM->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatiopPbOPy502SysM->GetYaxis()->SetTitleOffset(1.15);
		
		//hRatiopPbOPy502SysM->GetYaxis()->SetTitleOffset(1.7);
		hRatiopPbOPy502SysM->Divide(hDividePy502MNoErr);
		hRatiopPbOPy502SysM->SetFillStyle(hSypPb->GetFillStyle());
		hRatiopPbOPy502SysM->SetLineWidth(hSypPb->GetLineWidth());
		hRatiopPbOPy502SysM->SetLineColor(hSypPb->GetLineColor());
		hRatiopPbOPy502SysM->SetMarkerColor(hSypPb->GetLineColor());
		hRatiopPbOPy502SysM->SetMarkerStyle(hSypPb->GetMarkerStyle());
		//ratio pPb paper prop /PYTHIA
		TH1* hRatiopPbOPypappr;
		TH1* hDividePy502Mbis;
		UniformTH1FForDivide(hMpPbM, hMPy502M, hRatiopPbOPypappr, hDividePy502Mbis, "TH1D"); hRatiopPbOPypappr->SetName(TString::Format("hRatiopPbMOPypappr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOPypappr->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatiopPbOPypappr->GetYaxis()->SetTitleOffset(1.15);
		
		//hRatiopPbOPypappr->GetYaxis()->SetTitleOffset(1.7);
		hRatiopPbOPypappr->SetMarkerColor(hMpPbM->GetMarkerColor());
		hRatiopPbOPypappr->SetLineColor(hMpPbM->GetLineColor());
		hRatiopPbOPypappr->SetLineStyle(1);
		hRatiopPbOPypappr->SetMarkerStyle(24);
		hRatiopPbOPypappr->Divide(hDividePy502Mbis);
		hRatiopPbOPypappr->GetXaxis()->SetRangeUser(0., maxRangeMassFinal[ih]);
		hRatiopPbOPypappr->GetYaxis()->SetRangeUser(0., 4);
		
		// sys ratio pPb  paper prop /PYTHIA
		TH1* hRatiopPbOPy502Syspappr;
		TH1* hDivideSysPy502Mbis;
		UniformTH1FForDivide(hSypPbM, hMPy502M, hRatiopPbOPy502Syspappr, hDivideSysPy502Mbis, "TH1D");
		TH1* hDividePy502MbisNoErr = (TH1*)hDivideSysPy502Mbis->Clone(TString::Format("hDividePy502MbisNoErr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		for(Int_t ib = 0; ib<hDividePy502MbisNoErr->GetNbinsX(); ib++) hDividePy502MbisNoErr->SetBinError(ib, 0);
		
		hRatiopPbOPy502Syspappr->SetName(TString::Format("hRatiopPbOPy502Syspappr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOPy502Syspappr->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatiopPbOPy502Syspappr->GetYaxis()->SetTitleOffset(1.15);
		
		
		//hRatiopPbOPy502Syspappr->GetYaxis()->SetTitleOffset(1.7);
		hRatiopPbOPy502Syspappr->Divide(hDividePy502MbisNoErr);
		hRatiopPbOPy502Syspappr->SetFillStyle(hSypPbM->GetFillStyle());
		hRatiopPbOPy502Syspappr->SetLineWidth(hSypPbM->GetLineWidth());
		hRatiopPbOPy502Syspappr->SetLineColor(hSypPbM->GetLineColor());
		hRatiopPbOPy502Syspappr->SetFillColor(hSypPbM->GetLineColor());
		hRatiopPbOPy502Syspappr->SetMarkerStyle(1);
		
		
		// ratio pPb/Herwig
		
		TH1* hRatiopPbOHerwig;
		TH1* hDivideHerwig;
		UniformTH1FForDivide(hMpPb, hMHwpPb, hRatiopPbOHerwig, hDivideHerwig, "TH1D");
		hRatiopPbOHerwig->SetName(TString::Format("hRatiopPbOHerwig_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hDivideHerwig->SetName(TString::Format("hDivideHerwig_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOHerwig->Divide(hDivideHerwig);
		hRatiopPbOHerwig->SetMarkerStyle(hMHwpPb->GetMarkerStyle());
		hRatiopPbOHerwig->SetMarkerColor(hMHwpPb->GetMarkerColor());
		hRatiopPbOHerwig->SetLineColor(hMHwpPb->GetLineColor());
		hRatiopPbOHerwig->SetLineWidth(hMpPb->GetLineWidth());
		// systematics ratio pPb/Herwig
		TH1* hRatiopPbOHerwigSy;
		TH1* hDivideHerwigSy;
		UniformTH1FForDivide(hSypPb, hMHwpPb, hRatiopPbOHerwigSy, hDivideHerwigSy, "TH1D");
		hRatiopPbOHerwigSy->SetName(TString::Format("hRatiopPbOHerwigSy_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hDivideHerwigSy->SetName(TString::Format("hDivideHerwigSy_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOHerwigSy->Divide(hDivideHerwigSy);
		hRatiopPbOHerwigSy->SetFillStyle(hSypPb->GetFillStyle());
		hRatiopPbOHerwigSy->SetLineWidth(hSypPb->GetLineWidth());
		hRatiopPbOHerwigSy->SetLineColor(hMHwpPb->GetLineColor());
		hRatiopPbOHerwigSy->SetFillColor(hMHwpPb->GetLineColor());
		hRatiopPbOHerwigSy->SetMarkerStyle(hMHwpPb->GetMarkerStyle());
		hRatiopPbOHerwigSy->SetMarkerColor(hMHwpPb->GetMarkerColor());
		
		cRatiopPbOPy->cd(ih+1);
		gPad->SetTicks(1,1);
		hRatiopPbOPy502SysM->Draw("E2");
		hRatiopPbOPyM->Draw("sames");
		hRatiopPbOHerwigSy->Draw("E2sames");
		hRatiopPbOHerwig->Draw("sames");
		
		if(pPbpaperprop) {
			hRatiopPbOPy502Syspappr->Draw("E2sames");
			hRatiopPbOPypappr->Draw("sames");
			
		}
		lineOne->Draw();
		if(ih == 1) {
			latext.DrawLatex(4, 3, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 2.5, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
		}
		latext.DrawLatex(4, 3.6, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih]->Draw();
		
		if(ih == nhM-1){
			legRatiopPbOPy->AddEntry(hRatiopPbOPy502SysM, "p-Pb / PYTHIA(5.02TeV)");
			//legRatiopPbOPy->AddEntry(hRatiopPbOPy502SysM, "(Sys p-Pb) / PYTHIA(5.02TeV)", "F");
			legRatiopPbOPy->AddEntry(hRatiopPbOHerwigSy, "p-Pb / HERWIG(5.02TeV)");
			if(pPbpaperprop) {
				legRatiopPbOPy->AddEntry(hRatiopPbOPypappr, "pPb paper pr / PYTHIA(5.02TeV)", "PL");
				legRatiopPbOPy->AddEntry(hRatiopPbOPy502Syspappr, "(Sys pPb paper pr) / PYTHIA(5.02TeV)", optlegsim);
			}
			legRatiopPbOPy->Draw();
		}
		if(show1134) {}
		//ratio PbPb/PYTHIA
		TH1* hRatioPbPbOPy;
		TH1* hDividePy276;
		UniformTH1FForDivide(hMPbPb, hMPy276, hRatioPbPbOPy, hDividePy276, "TH1D"); hRatioPbPbOPy->SetName(TString::Format("hRatiopPbOPy_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		hRatioPbPbOPy->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatioPbPbOPy->GetYaxis()->SetTitleOffset(1.15);
		hRatioPbPbOPy->Divide(hDividePy276);
		hRatioPbPbOPy->GetXaxis()->SetRangeUser(0., maxM);
		hRatioPbPbOPy->GetYaxis()->SetRangeUser(0., 4);
		
		
		// sys ratio PbPb/PYTHIA
		TH1* hRatioPbPbOPy276Sys;
		UniformTH1FForDivide(hSyPbPb, hMPy276, hRatioPbPbOPy276Sys, hDividePy276, "TH1D");
		TH1* hDividePy276NoErr = (TH1*)hDividePy276->Clone(TString::Format("hDividePy276NoErr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		for(Int_t ib = 0; ib<hDividePy276NoErr->GetNbinsX(); ib++) hDividePy276NoErr->SetBinError(ib, 0);
		
		hRatioPbPbOPy276Sys->SetName(TString::Format("hRatioPbPbOPy276Sys_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatioPbPbOPy276Sys->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatioPbPbOPy276Sys->GetYaxis()->SetTitleOffset(1.15);
		
		//hRatioPbPbOPy276Sys->GetYaxis()->SetTitleOffset(1.7);
		hRatioPbPbOPy276Sys->Divide(hDividePy276NoErr);
		hRatioPbPbOPy276Sys->SetFillStyle(1001);
		hRatioPbPbOPy276Sys->SetLineWidth(2);
		hRatioPbPbOPy276Sys->SetLineColor(hSyPbPb->GetLineColor());
		hRatioPbPbOPy276Sys->SetFillColor(hSyPbPb->GetFillColor());
		hRatioPbPbOPy276Sys->SetMarkerColor(hSyPbPb->GetMarkerColor());
		hRatioPbPbOPy276Sys->SetMarkerStyle(hSyPbPb->GetMarkerStyle());
		
		cRatioPbPbOPy->cd(ih+1);
		
		gPad->SetTicks(1,1);
		hRatioPbPbOPy276Sys->GetYaxis()->SetRangeUser(0., 4);
		hRatioPbPbOPy276Sys->Draw("E2");
		hRatioPbPbOPy->Draw("sames");
		lineOne->Draw();
		if(ih == nhM-1){
			legRatioPbPbOPy->AddEntry(hRatioPbPbOPy276Sys, "Pb-Pb / PYTHIA(2.76TeV)");
			//legRatioPbPbOPy->AddEntry(hRatioPbPbOPy276Sys, "Sys Pb-Pb / PYTHIA(2.76TeV)", "F");
			legRatioPbPbOPy->Draw();
		}
		//pvpt[ih]->Draw();
		latext.DrawLatex(4, 3.6, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		if(ih == 1) {
			latext.DrawLatex(4, 3, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 2.5, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
		}
		// data/ PYTHIA overlapped
		cRatioDataOPy->cd(ih+1);
		
		gPad->SetTicks(1,1);
		hRatiopPbOPy502SysM->GetYaxis()->SetRangeUser(0., 4.);
		hRatiopPbOPy502SysM->GetXaxis()->SetRangeUser(0., maxM);
		hRatioPbPbOPy276Sys->GetYaxis()->SetRangeUser(0., 4.);
		hRatioPbPbOPy276Sys->GetXaxis()->SetRangeUser(0., maxM);
		hRatioPbPbOPy276Sys->Draw("samesE2");
		hRatioPbPbOPy->Draw("sames");
		hRatiopPbOPy502SysM->Draw("samesE2");
		hRatiopPbOPyM->Draw("sames");
		latext.DrawLatex(4, 3.6, TString::Format("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		//pvpt[ih]->Draw();
		lineOne->Draw();
		if(ih == 1) {
			latext.DrawLatex(4, 3, "Charged jets, Anti-#it{k}_{T}");
			latext.DrawLatex(4.5, 2.5, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
		}
		if(ih == nhM-1){
			legRatioDataOPy->AddEntry(hRatiopPbOPy502SysM, "p-Pb / PYTHIA(5.02TeV)");
			//legRatioDataOPy->AddEntry(hRatiopPbOPy502SysM, "Sys p-Pb / PYTHIA(5.02TeV)", "F");
			legRatioDataOPy->AddEntry(hRatioPbPbOPy276Sys, "Pb-Pb / PYTHIA(2.76TeV)");
			//legRatioDataOPy->AddEntry(hRatioPbPbOPy276Sys, "Sys Pb-Pb / PYTHIA(2.76TeV)", "F");
			legRatioDataOPy->Draw();
		}
		
		//write into a root file
		fOutput->cd();
		hRatioPbPbOpPbSys     ->Write();
		hRatioPbPbOpPb        ->Write();
		hRatio276O502M        ->Write();
		hSyPbPb               ->Write();
		hMPbPb                ->Write();
		hMPy276               ->Write();
		hSypPb                ->Write();
		hMpPb                 ->Write();
		hMPy502M              ->Write();
	}
	cMass->cd(1);
	DrawLogo(0, 12, 0.18, "", 42, "");
	cMasspPb->cd(1);
	DrawLogo(0, 12, 0.18, "", 42, "");
	cMassPbPb->cd(1);
	DrawLogo(0, 12, 0.18, "", 42, "");
	cMassPbPbOnly->cd(1);
	DrawLogo(0, 12, 0.18, "", 42, "");
	cMassPbPbPy->cd(1);  
	DrawLogo(0, 12, 0.18, "", 42, "");
	cRatioPbPbOpPb->cd(1);
	DrawLogo(0, 4, 0.2, "", 42, "");
	cRatioDataOPy->cd(1);
	DrawLogo(0, 4, 0.2, "", 42, "");
	cRatioPbPbOPy->cd(1);
	DrawLogo(0, 4, 0.2, "", 42, "");
	cRatiopPbOPy->cd(1); 
	DrawLogo(0, 4, 0.2, "", 42, "");
	
	cRatiopPbOPy->cd(1);
	latext.DrawLatex(2.8, 3, "p-Pb #sqrt{#it{s}}_{NN} = 5.02 TeV");
	//pvSystpPb->Draw();
	//cRatiopPbOPy->cd(2);
	latext.DrawLatex(2.8, 2.5, "PYTHIA Perugia 2011");
	//pvSystPyth->Draw();
	
	cRatioDataOPy->cd(1);
	latext.DrawLatex(2.8, 3., "Pb-Pb 0-10% #sqrt{#it{s}}_{NN} = 2.76 TeV");
	latext.DrawLatex(2.8, 2.5, "p-Pb #sqrt{#it{s}}_{NN} = 5.02 TeV");
	//pvSystpPb->Draw();
	latext.DrawLatex(2.8, 2, "PYTHIA Perugia 2011");
	//pvSystPyth->Draw();
	
	cRatioPbPbOPy->cd(1);
	latext.DrawLatex(2.8, 3., "Pb-Pb 0-10% #sqrt{#it{s}}_{NN} = 2.76 TeV");
	latext.DrawLatex(2.8, 2.5, "PYTHIA Perugia 2011");
	
	cRatioPbPbOpPb->cd(1);
	latext.DrawLatex(2.8, 3., "Pb-Pb 0-10% #sqrt{#it{s}}_{NN} = 2.76 TeV");
	latext.DrawLatex(4.1, 2.5, "p-Pb #sqrt{#it{s}}_{NN} = 5.02 TeV");
	latext.DrawLatex(4.1, 2, "PYTHIA Perugia 2011");
	//pvSystPbPb->Draw();
	
	SaveCv(cMass);
	SaveCv(cRelUnc);
	SaveCv(cMasspPb);
	SaveCv(cMassPbPb);
	SaveCv(cMassPbPbOnly);
	SaveCv(cMassPbPbPy);
	SaveCv(cMassPbPbModels);
	SaveCv(cRatioPbPbOpPb);
	SaveCv(cRatioPbPbOPy);
	SaveCv(cRatiopPbOPy);
	SaveCv(cRatioDataOPy);
}

//_____________________________________________________________________________

TH1D* CalculateMean(Double_t **xrange, Int_t nbins, Double_t lims[], TH1D** hdistr, TGraphErrors *&grmean, TH1D *&hrms, const char* namehmean){
	
	Double_t mean[nbins]; 
	Double_t merr[nbins];
	Double_t rms[nbins];
	Double_t rer[nbins];
	Double_t deltalims[nbins];
	
	TH1D* hmean = new TH1D(namehmean, "Mean;#it{p}_{T}(GeV/#it{c}); #LT M_{ch jet}#GT (GeV/#it{c}^{2})", nbins, lims);
	hmean->Sumw2();
	
	hrms = new TH1D(TString::Format("%sStdDev", namehmean), "StdDev;#it{p}_{T}(GeV/#it{c}); Std Dev(GeV/#it{c}^{2})", nbins, lims);
	hrms->Sumw2();
	
	for(Int_t ib = 0; ib < nbins; ib++){
		mean[ib] = 0;
		merr[ib] = 0;
		rms [ib] = 0;
		deltalims[ib] = (lims[ib+1] - lims[ib])/2.;
		if(!hdistr[ib]) continue;
		hdistr[ib]->GetXaxis()->SetRange(hdistr[ib]->GetXaxis()->FindBin(xrange[ib][0]), hdistr[ib]->GetXaxis()->FindBin(xrange[ib][1]));
		mean[ib] = hdistr[ib]->GetMean();
		merr[ib] = hdistr[ib]->GetMeanError();
		rms [ib] = hdistr[ib]->GetRMS();
		rer [ib] = hdistr[ib]->GetRMSError();
		Printf("Result[%d] = %f+-%f", ib, mean[ib], merr[ib]);
		hmean->SetBinContent(ib+1, mean[ib]);
		hmean->SetBinError(  ib+1, merr[ib]);
		
		hrms->SetBinContent(ib+1, rms[ib]);
		hrms->SetBinError(  ib+1, rer[ib]);
	}
	
	grmean = new TGraphErrors(nbins, lims, mean, deltalims, merr);
	
	return hmean;

}

//_____________________________________________________________________________

TH1D* CalculateSysMean(Double_t **xrange, Int_t nbins, Double_t lims[], TH1D** hSysdistr, TGraphErrors *&grSysmean, TH1D *&hSysrms, const char* namehsysmean){
	
	Double_t mean[nbins]; 
	Double_t merr[nbins];
	Double_t rms[nbins];
	Double_t rer[nbins];
	Double_t deltalims[nbins];
	
	TH1D* hSysmean = new TH1D(namehsysmean, "Mean;#it{p}_{T}(GeV/#it{c}); <M_{ch jet}>(GeV/#it{c}^{2})", nbins, lims);
	hSysmean->Sumw2();
	hSysmean->SetFillStyle(3001);
	
	hSysrms = new TH1D(TString::Format("%sStdDev", namehsysmean), "StdDev;#it{p}_{T}(GeV/#it{c}); Std Dev(GeV/#it{c}^{2})", nbins, lims);
	hSysrms->Sumw2();
	
	for(Int_t ib = 0; ib < nbins; ib++){
		mean[ib] = 1;
		merr[ib] = 0;
		rms [ib] = 0;
		deltalims[ib] = (lims[ib+1] - lims[ib])/2.;
		if(!hSysdistr[ib]) continue;
		hSysdistr[ib]->GetXaxis()->SetRange(hSysdistr[ib]->GetXaxis()->FindBin(xrange[ib][0]), hSysdistr[ib]->GetXaxis()->FindBin(xrange[ib][1]));
		//Int_t massbins = hSysdistr[ib]->GetNbinsX();
		//Double_t tot = 0;
		//for(Int_t j = 0; j< massbins; j++){
		//	merr[ib] += hSysdistr[ib]->GetBinError(j+1);
		//	if(hSysdistr[ib]->GetBinError(j+1) > 1e-05) tot++;
		//	Printf("+ %f", hSysdistr[ib]->GetBinError(j+1));
		//}
		//Printf("%f/%d", merr[ib], massbins);
		//merr[ib] /= tot;
		//
		mean[ib] = hSysdistr[ib]->GetMean();
		merr[ib] = hSysdistr[ib]->GetMeanError();
		rms [ib] = hSysdistr[ib]->GetRMS();
		rer [ib] = hSysdistr[ib]->GetRMSError();
		Printf("Result[%d] = %f+-(sys)%f", ib, mean[ib], merr[ib]);
		hSysmean->SetBinContent(ib+1, mean[ib]);
		hSysmean->SetBinError(  ib+1, merr[ib]);
		
		hSysrms->SetBinContent(ib+1, rms [ib]);
		hSysrms->SetBinError(  ib+1, rer[ib]);
	}
	
	grSysmean = new TGraphErrors(nbins, lims, mean, deltalims, merr);
	
	hSysmean->SetMarkerStyle(24);
	hSysmean->SetFillStyle(3001);
	hSysrms->SetMarkerStyle(24);
	hSysrms->SetFillStyle(3001);	
	return hSysmean;

}

//_____________________________________________________________________________

TH1D* GetMeanAndSystFromFile(TString filename, TH1D*& hSyst){
	
	TFile *fin = new TFile(filename);
	if(!fin->IsOpen()){
		Printf("File %s not found", filename.Data());
		return 0;
	}
	
	hSyst = GetSystFromFile(fin);

	return GetMeanFromFile(fin);
}

//_____________________________________________________________________________

TH1D* GetMeanFromFile(TString filename){
	TFile *fin = new TFile(filename);
	if(!fin->IsOpen()){
		Printf("File %s not found", filename.Data());
		return 0;
	}
	return GetMeanFromFile(fin);
}

//_____________________________________________________________________________

TH1D* GetMeanSystFromFile(TString filename){
	TFile *fin = new TFile(filename);
	if(!fin->IsOpen()){
		Printf("File %s not found", filename.Data());
		return 0;
	}
	return GetSystFromFile(fin);
}

//_____________________________________________________________________________

TH1D* GetMeanFromFile(TFile *fin){
	if(!fin->IsOpen()){
		Printf("File not found");
		return 0;
	}
	TString hnameMas = "hMeanCentral";
	TH1D* hMass = (TH1D*)fin->Get(hnameMas);
	
	return hMass;
}

//_____________________________________________________________________________

TH1D* GetSystFromFile(TFile *fin){
	if(!fin->IsOpen()){
		Printf("File not found");
		return 0;
	}
	TString hnameSys = "hMeanSystTot";
	TH1D* hSyst = (TH1D*)fin->Get(hnameSys);
	
	return hSyst;
}
//_____________________________________________________________________________

void DrawMeanComparison(Bool_t stylezero = kTRUE){
	
	if(stylezero) {
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		gStyle->SetTextFont(42);
	}
	
	const Int_t ninputs = 4;
	TString inputDistr[ninputs] = {
		"/data/Work/jets/JetMass/PbPbResults/TranformedIntoTH1.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160916/MasspPbResults.root",
		"/data/Work/jets/JetMass/PbPbResults/JetMassPerugia2011_ecms2760.root",
		"/data/Work/jets/JetMass/PbPbJEWEL/alice_jetmass_histograms.root"
	}; 
	TString h1names[ninputs] = {
		"hUnfMass_PtBin",
		"hUnfM_Itr3_Pt",
		"hM_",
		"Mass_centbin1_ptbin"
	}; 
	
	const Int_t ninputsSy = 2;
	TString inputSyst[ninputsSy] = {
		"/data/Work/jets/JetMass/PbPbResults/TranformedIntoTH1.root",
		//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160911/MasspPbResults.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160916/TotalSystematicUnc.root"
	};
	TString hsysnames[ninputsSy] = {
		"hUnfMassSyst_PtBin",
		//"hSystTot_Pt",
		"hMeanSystTot"
	};
	
	//3 bins in this output for pPb, PbPb has 4 bins and it starts from bin 1, so we need an offset of 2
	Int_t offset[ninputs] = {
		2,
		0
	};
	TString leg[ninputs] = {
		"Pb-Pb 0-10\% #sqrt{#it{s}_{NN}} = 2.76 TeV",
		"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV",
		"PYTHIA #sqrt{#it{s}_{NN}} = 2.76 TeV",
		"JEWEL+PYTHIA, Pb-Pb 0-10\% #sqrt{#it{s}_{NN}} = 2.76 TeV"
	};
	Int_t  mrk[ninputs] = {
		20,
		21,
		25,
		28
	};
	Int_t clr[ninputs] = {
		kBlue+2,
		kOrange+8,
		kRed+2,
		kOrange+1
	};
	Int_t clrsy[ninputsSy] = {
		kBlue-10,
		kOrange+8,
	};
	Int_t fillsy[ninputsSy] = {
		1001,
		0
	};
	
	Double_t **xrange = new Double_t*[nptbins];
	for(Int_t ipt = 0; ipt<nptbins; ipt++){
		xrange[ipt] = new Double_t[2];
		xrange[ipt][0] = 0;
		xrange[ipt][1] = maxRangeMassFinal[ipt];
	}
	
	DrawMeanComparison(xrange, ninputs, inputDistr, h1names, ninputsSy, inputSyst, hsysnames, offset, leg, mrk, clr, clrsy, fillsy);
	
	// do the same for the pPb raw and the embedded?one has background, the other not. how to define the variation we need for the prior?
	//we can use the particle level ratio from the RatiopPb method
}

//_____________________________________________________________________________

void DrawMeanComparison(Double_t **xrange, const Int_t ninputs, TString inputDistr[], TString h1names[], const Int_t ninputsSy, TString inputSyst[], TString hsysnames[], Int_t offset[], TString legtx[], Int_t mrk[], Int_t clr[], Int_t clrsy[], Int_t fillsy[]){
	// need to add the tratment of the systematics
		
	// text for the final plots
	TPaveText *pvgeneral = new TPaveText(0.45, 0.7, 0.65, 0.8, "NDC");
	pvgeneral->SetFillStyle(0);
	pvgeneral->SetBorderSize(0);
	pvgeneral->AddText("Anti-#it{k}_{T}, #it{R} = 0.4");
   	   
	TH1D* hmean[ninputs];
	TGraphErrors *grmean[ninputs];
	TH1D *hrms[ninputs];
	
	TH1D* hSysmean[ninputs];
	TGraphErrors *grSysmean[ninputs];
	TH1D *hSysrms[ninputs];
	
	TCanvas **cInput = new TCanvas*[ninputs];
	for(Int_t ifi = 0; ifi< ninputs; ifi++){
		cInput[ifi] = new TCanvas(TString::Format("c%s", h1names[ifi].Data()), TString::Format("%s", h1names[ifi].Data()), 800, 800);
		cInput[ifi]->Divide(3,1);
	}
	TCanvas* cMean = new TCanvas("cMean", "Mean", 800, 800);
	TCanvas* cRms = new TCanvas("cRms", "Standard Deviation", 800, 800);
	
	TLegend *leg = new TLegend(0.2, 0.6, 0.7, 0.75);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	TLegend *legSy = new TLegend(0.7, 0.68, 0.9, 0.75);
	legSy->SetBorderSize(0);
	legSy->SetFillStyle(0);
	Bool_t drawn = kFALSE;
	
	
	for(Int_t ifi = 0; ifi< ninputsSy; ifi++){
	TFile* finS = 0x0;
		if(!inputSyst[ifi].IsNull()){
			finS = new TFile(inputSyst[ifi]);
			if(!finS->IsOpen()){
				Printf("File %s not found", inputSyst[ifi].Data());
				continue;
			}
		} else {
			Printf("No Systematics");
		}
		
		if(hsysnames[ifi] == "hMeanSystTot"){
			hSysmean[ifi] = (TH1D*)finS->Get(hsysnames[ifi]);
			
			
		} else{
			Bool_t isdifferentnamingscheme = kFALSE;
			TString hnamesuff = "", hnamepref = "";
			if(h1names[ifi].Contains("Iter")){
				isdifferentnamingscheme = kTRUE;
				hnamesuff = h1names[ifi](6,11);
				hnamepref = h1names[ifi](0, 6);
				Printf("Pref %s Suff %s -> %d", hnamepref.Data(), hnamesuff.Data(), isdifferentnamingscheme);
			}
			TH1D **hSysdistr = new TH1D*[nptbins];
			
			for(Int_t ib = 0; ib<nptbins; ib++){
				if(!finS) continue;
				//systematics
				TString hnamethisloop = "", hnamedoublethisloop = "";
				hnamethisloop =  TString::Format("%s%d", hsysnames[ifi].Data(),  ib + offset[ifi]);
				hnamedoublethisloop = TString::Format("%s%.0f_%.0fb", hsysnames[ifi].Data() ,  ptlims[ib + offset[ifi]],  ptlims[ib + 1 + offset[ifi]]);
				
				hSysdistr[ib] = (TH1D*)finS->Get(hnamethisloop);
				if(!hSysdistr[ib]){
					Printf("Histogram %s not found, try %s ", hnamethisloop.Data(), hnamedoublethisloop.Data());
					hSysdistr[ib] = (TH1D*)finS->Get(hnamedoublethisloop);
					if(!hSysdistr[ib]) {
						Printf("Histogram %s and %s not found ", hnamethisloop.Data(), hnamedoublethisloop.Data());
						finS->ls();
						continue;
					}
				}
				
				cInput[ifi]->cd(ib+1);
				hSysdistr[ib]->Draw("E2");
			}	
      	 	
			//call the method that calculates the mean systematics
			
			hSysmean[ifi] = CalculateSysMean(xrange, nptbins, ptlims, hSysdistr, grSysmean[ifi], hSysrms[ifi], TString::Format("hSysMeanMass_%d_", ifi));
			
			hSysmean[ifi]->SetFillColor(clrsy[ifi]);
			hSysmean[ifi]->SetMarkerStyle(1);
			hSysmean[ifi]->SetFillStyle(fillsy[ifi]);
			hSysmean[ifi]->SetLineStyle(1);
			hSysmean[ifi]->SetLineWidth(2);
			hSysmean[ifi]->SetLineColor(clrsy[ifi]);
			
			hSysrms[ifi]->SetFillColor(clrsy[ifi]);
			hSysrms[ifi]->SetMarkerStyle(1);
			hSysrms[ifi]->SetFillStyle(fillsy[ifi]);
			//hSysrms[ifi]->SetFillStyle(3001);
			hSysrms[ifi]->SetLineWidth(2);
			hSysrms[ifi]->SetLineStyle(1);
			hSysrms[ifi]->SetLineColor(clrsy[ifi]);
		
		//draw all
		legSy->AddEntry(hSysmean[ifi], TString::Format("Systematics"), "F");
		
		if(hSysmean[ifi] && !drawn) {
			cMean->cd();
			hSysmean[ifi]->GetYaxis()->SetRangeUser(2, 30);
			hSysmean[ifi]->Draw("E2");
			drawn = kTRUE;
			
			cRms->cd();
			if(hSysrms[ifi]) {
				hSysrms[ifi]->GetYaxis()->SetRangeUser(0, 6);
				hSysrms[ifi]->Draw("E2");
			}
		} else {
			cMean->cd();
			if(hSysmean[ifi]) hSysmean[ifi]->Draw("E2sames");
			cRms->cd();
			if(hSysrms[ifi]) hSysrms[ifi]->Draw("E2sames");
		}
		}
	}
	
	
	for(Int_t ifi = 0; ifi< ninputs; ifi++){
		
		hmean [ifi] = 0x0;
		grmean[ifi] = 0x0;
		hrms  [ifi] = 0x0;
		TFile* fin = new TFile(inputDistr[ifi]);
		if(!fin->IsOpen()){
			Printf("File %s not found", inputDistr[ifi].Data());
			continue;
		}
		
		Bool_t isdifferentnamingscheme = kFALSE;
		TString hnamesuff = "", hnamepref = "";
		if(h1names[ifi].Contains("Iter")){
			isdifferentnamingscheme = kTRUE;
			hnamesuff = h1names[ifi](6,11);
			hnamepref = h1names[ifi](0, 6);
			Printf("Pref %s Suff %s -> %d", hnamepref.Data(), hnamesuff.Data(), isdifferentnamingscheme);
		}
		
		TH1D **h1Mass    = new TH1D*[nptbins];
		
		for(Int_t ib = 0; ib<nptbins; ib++){
			
			//mass
			TString hnamethisloop = "", hnamedoublethisloop = "";
			if(isdifferentnamingscheme) {
				Printf("isdifferentnamingscheme");
				hnamethisloop = TString::Format("%s%d%s", hnamepref.Data(),  ib + offset[ifi], hnamesuff.Data());
				
				// for the comparison of the unfolded new
				
			}
			else {
				hnamethisloop = TString::Format("%s%d", h1names[ifi].Data(),  ib + offset[ifi]);
				hnamedoublethisloop = TString::Format("%s%.0f_%.0f", h1names[ifi].Data() ,  ptlims[ib + offset[ifi]],  ptlims[ib + 1 + offset[ifi]]);
			}
			Printf("Reading %s", hnamethisloop.Data());
			h1Mass[ib] = (TH1D*)fin->Get(hnamethisloop);
			
			if(!h1Mass[ib]){
				Printf("Histogram %s not found, try %s ", hnamethisloop.Data(), hnamedoublethisloop.Data());
				h1Mass[ib] = (TH1D*)fin->Get(hnamedoublethisloop);
				if(!h1Mass[ib]) {
					Printf("Histogram %s and %s not found ", hnamethisloop.Data(), hnamedoublethisloop.Data());
					fin->ls();
					continue;
				}
      	 	}
      	 	cInput[ifi]->cd(ib+1);
      	 	h1Mass[ib]->Draw("sames");
      	 	
      	 	
		}
		
		//call the method that calculates the mean
		
		hmean[ifi] = CalculateMean(xrange, nptbins, ptlims, h1Mass, grmean[ifi], hrms[ifi], TString::Format("hMeanMass_%d_", ifi));
		hmean[ifi]->SetMarkerStyle(mrk[ifi]);
		hmean[ifi]->SetMarkerColor(clr[ifi]);
		hmean[ifi]->SetLineColor(hmean[ifi]->GetMarkerColor());
		
		hrms[ifi]->SetMarkerStyle(mrk[ifi]);
		hrms[ifi]->SetMarkerColor(clr[ifi]);
		hrms[ifi]->SetLineColor(clr[ifi]);
		leg->AddEntry(hmean[ifi], legtx[ifi], "PL");
		
		if(hsysnames[ifi] == "hMeanSystTot"){
			SetMassValueInSystematic(hSysmean[ifi], hmean[ifi]);
			cMean->cd();
			if(!drawn) hSysmean[ifi]->Draw("E2");
			else hSysmean[ifi]->Draw("E2sames");
		}
		//draw all
		
		if(hmean[ifi] && hSysmean[ifi] && !drawn) {
			cMean->cd();
			hmean[ifi]->Draw("sames");
			drawn = kTRUE;
			
			cRms->cd();
			if(hrms[ifi]) {
				hrms[ifi]->Draw("sames");
			}
		} else {
			cMean->cd();
			if(hmean[ifi]) hmean[ifi]->Draw("sames");
			cRms->cd();
			if(hrms[ifi]) hrms[ifi]->Draw("sames");
		}
		
		
	}
	
	cMean->cd();
	leg->Draw();
	legSy->Draw();
	pvgeneral->Draw();
	DrawLogo(1, 0.6, 0.8, 0.9, 0.9, "", 42, ""); // 
	cRms->cd();
	leg->Draw();
	legSy->Draw();
	pvgeneral->Draw();
	
	
	SaveCv(cMean);
	SaveCv(cRms);
}

//_____________________________________________________________________________

void DrawMeanFromSystOutput(Bool_t stylezero = kTRUE){
	
	// this is the method to be used for the paper
	
	if(stylezero) {
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		
	}
	gStyle->SetTextFont(42);
	
	gStyle->SetLabelSize(0.05, "X");
	gStyle->SetTitleOffset(1.1, "X");
	gStyle->SetTitleSize(0.048, "X");
	
	gStyle->SetLabelSize(0.05, "Y");
	gStyle->SetTitleOffset(1.34, "Y");
	gStyle->SetTitleSize(0.048, "Y");
	
	gStyle->SetNdivisions(1004, "Y");
	
	gStyle->SetPadBottomMargin(.13);
	gStyle->SetPadLeftMargin(.16);
	gStyle->SetPadRightMargin(.06);
	
	//gStyle->SetErrorX(.2);
	
	TString filenamepPb =
	// results first version
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160916/MasspPbResults.root";
	
	// pPb mean jet mass
	TH1D* hSystpPb = 0x0;
	TH1D* hMeanjmpPb = GetMeanAndSystFromFile(filenamepPb, hSystpPb);
	hSystpPb->SetMarkerStyle(hMeanjmpPb->GetMarkerStyle());
	hSystpPb->SetMarkerColor(hMeanjmpPb->GetMarkerColor());
	hSystpPb->SetMarkerSize(2);
	
	// sqrts dependence
	TH1D* hsqrtdiff = SystEnergyDep(filenamepPb);
	
	TH1D** hMPbPb = new TH1D*[nptbins];
	TH1D** hMSysPbPb = new TH1D*[nptbins];
	TList* listPbPb = GetPbPbResults();
	Int_t offsetPbPb = 1;
	
	TList* list502 = GetPythiapPbMarta();
	TH1D** hpythia502 = new TH1D*[nptbins];
	
	TList* list276 = GetPythiaPbPb();
	TH1D** hpythia276 = new TH1D*[nptbins];
	
	TList* listJewel = GetJEWELPbPb();
	TH1D** hJewel = new TH1D*[nptbins];
	
	TList* listqpy = GetJEWELQPYTHIA(3);
	TH1D** hQPythia = new TH1D*[nptbins];
	Int_t offsetQPy = 0; // was 2
	
	Double_t **xrange = new Double_t*[nptbins];
	for(Int_t ipt = 0; ipt < nptbins; ipt++){
		xrange[ipt] = new Double_t[2];
		xrange[ipt][0] = 0;
		xrange[ipt][1] = maxRangeMassFinal[ipt];
		
		hMPbPb[ipt] = (TH1D*)listPbPb->At(2*(ipt+offsetPbPb));
		//needed for style
		hMPbPb[ipt]->SetMarkerSize(2);
		hMSysPbPb[ipt] = (TH1D*)listPbPb->At(2*(ipt+offsetPbPb)+1);
		hMSysPbPb[ipt]->SetMarkerSize(2);
		
		hpythia502[ipt] = (TH1D*)list502->At(ipt+offsetPbPb);
		hpythia276[ipt] = (TH1D*)list276->At(ipt+offsetPbPb);
		
		hJewel[ipt] = (TH1D*)listJewel->At(ipt);
		
		hQPythia[ipt] = (TH1D*)listqpy->At(ipt+offsetQPy);
		
	}
	
	TGraphErrors *grmean = 0x0;
	TH1D *hrms = 0x0;
	// PbPb mean jet mass
	TH1D* hMeanPbPb = CalculateMean(xrange, nptbins, ptlims, hMPbPb, grmean, hrms, "hMeanjmassPbPb");
	hMeanPbPb->SetMarkerStyle(hMPbPb[0]->GetMarkerStyle());
	hMeanPbPb->SetMarkerColor(hMPbPb[0]->GetMarkerColor());
	hMeanPbPb->SetLineColor(hMPbPb[0]->GetLineColor());
	hMeanPbPb->SetMarkerSize(2);
	
	TFile *finPbPbSys = new TFile("/data/Work/jets/JetMass/PbPbResults/MeanJetMassAllSyst_Area.root");
	TGraphErrors *grMeanPbPbSys = (TGraphErrors*)finPbPbSys->Get("grUnfMeanSyst");
	grMeanPbPbSys->SetFillColor(hMSysPbPb[0]->GetFillColor());
	grMeanPbPbSys->SetLineColor(hMSysPbPb[0]->GetLineColor());
	grMeanPbPbSys->SetMarkerStyle(hMeanPbPb->GetMarkerStyle());
	grMeanPbPbSys->SetMarkerColor(hMeanPbPb->GetMarkerColor());
	grMeanPbPbSys->SetMarkerSize(2);
	
	// PYTHIA mean jet mass
	TH1D* hMeanPy276 = CalculateMean(xrange, nptbins, ptlims, hpythia276, grmean, hrms, "hMeanjmassPy276");
	hMeanPy276->SetMarkerStyle(hpythia276[0]->GetMarkerStyle());
	hMeanPy276->SetMarkerSize(2);
	hMeanPy276->SetMarkerColor(hpythia276[0]->GetMarkerColor());
	hMeanPy276->SetLineColor  (hpythia276[0]->GetLineColor());
	
	TH1D* hMeanPy502 = CalculateMean(xrange, nptbins, ptlims, hpythia502, grmean, hrms, "hMeanjmassPy502");
	
	// JEWEL mean jet mass
	TH1D* hMeanJewel = CalculateMean(xrange, nptbins, ptlims, hJewel, grmean, hrms, "hMeanjmassJewel");
	hMeanJewel->SetMarkerStyle(hJewel[0]->GetMarkerStyle());
	hMeanJewel->SetMarkerColor(hJewel[0]->GetMarkerColor());
	hMeanJewel->SetLineColor  (hJewel[0]->GetLineColor());
	hMeanJewel->SetMarkerSize(2);
	
	// QPYTHIA mean jet mass
	TH1D* hMeanQPythia = CalculateMean(xrange, nptbins, ptlims, hQPythia, grmean, hrms, "hMeanjmassQPythia");
	hMeanQPythia->SetMarkerStyle(hQPythia[0]->GetMarkerStyle());
	hMeanQPythia->SetMarkerColor(hQPythia[0]->GetMarkerColor());
	hMeanQPythia->SetLineColor  (hQPythia[0]->GetLineColor());
	hMeanQPythia->SetMarkerSize(2);
	
	TCanvas *cMeanjetmass = new TCanvas("cMeanjetmass", "Mean jet mass", 600, 600);
	TCanvas *cMeanjetmassdata = new TCanvas("cMeanjetmassdata", "Mean jet mass", 600, 600);
	
	TLegend *leg = new TLegend(0.45, 0.15, 0.95, 0.4);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry(hMeanPbPb, "Pb-Pb 0-10%");
	leg->AddEntry(hMeanPy276, "PYTHIA Perugia11");
	leg->AddEntry(hMeanJewel, "JEWEL+PYTHIA recoil on");
	leg->AddEntry(hMeanQPythia, "Q-PYTHIA");
	
	TLegend *legda = new TLegend(0.45, 0.15, 0.95, 0.4);
	legda->SetBorderSize(0);
	legda->SetFillStyle(0);
	legda->AddEntry(grMeanPbPbSys, "Pb-Pb 0-10% #sqrt{s} = 2.76 TeV");
	legda->AddEntry(hSystpPb, "p-Pb #sqrt{s} = 5.02 TeV");
	legda->AddEntry(hsqrtdiff, "#sqrt{s} difference");
	
	for(Int_t ipt = 0; ipt < nptbins; ipt++){
		cMeanjetmass->cd();
		gPad->SetTicks(1,1);
		grMeanPbPbSys->GetYaxis()->SetRangeUser(0., 22.);
		grMeanPbPbSys->GetXaxis()->SetRangeUser(ptlims[0], ptlims[nptbins]);
		grMeanPbPbSys->SetTitle("; #it{p}_{T, ch jet} (GeV/#it{c}); #LT #it{M}_{ch jet} #GT (GeV/#it{c}^{2})");
		grMeanPbPbSys->DrawClone("AE2");
		hMeanPy276->Draw("sames");
		hMeanPbPb->Draw("sames");
		hMeanJewel->Draw("sames");
		hMeanQPythia->Draw("sames");
		leg->Draw();
		
		
		cMeanjetmassdata->cd();
		gPad->SetTicks(1,1);
		grMeanPbPbSys->GetYaxis()->SetRangeUser(0., 16.);
		grMeanPbPbSys->Draw("AE2");
		hsqrtdiff->Draw("PE2sames");
		hMeanjmpPb->Draw("sames");
		hSystpPb->Draw("E2sames");
		hMeanPbPb->Draw("sames");
		legda->Draw();
	}
	
	TLatex latext;
	latext.SetTextSize(0.05);
	// text for the final plots
	TPaveText *pvgeneral = new TPaveText(0.2, 0.18, 0.4, 0.3, "NDC");
	pvgeneral->SetFillStyle(0);
	pvgeneral->SetBorderSize(0);
	pvgeneral->AddText("Charged jets");
	pvgeneral->AddText("Anti-#it{k}_{T}, #it{R} = 0.4");
	
	
	cMeanjetmass->cd();
	DrawLogo(0, 65, 2, "", 42, "");
	latext.DrawLatex(65, 20, "Charged jets, Anti-#it{k}_{T}");
	latext.DrawLatex(65, 18, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
	//pvgeneral->Draw();
	cMeanjetmassdata->cd();
	DrawLogo(0, 65, 2, "", 42, "");
	latext.DrawLatex(65, 14, "Charged jets, Anti-#it{k}_{T}");
	latext.DrawLatex(65, 12, "#it{R} = 0.4, |#eta_{jet}| < 0.5");
	//pvgeneral->Draw();
	
	SaveCv(cMeanjetmass);
	SaveCv(cMeanjetmassdata);
}

//_____________________________________________________________________________

TH1D* SystEnergyDep(TString filenamepPb){
	
	// get histograms pythia at the two energies
	TList* list502 = GetPythiapPbMarta();
	TList* list276 = GetPythiaPbPb();
	
	// put the systematic for sqrt{s} on top of these points
	TH1D* hMeanpPb = GetMeanFromFile(filenamepPb);
	if(!hMeanpPb) return 0;
	
	TH1D* hMean502 = new TH1D("hMean502", "Mean jet mass vs #it{p}_{T} #sqrt{s} = 5.02 TeV; #it{p}_{T} (GeV/#it{c}); #LT #it{M} #GT (GeV/#it{c}^{2})", nptbins, ptlims);
	
	TH1D* hMean276 = new TH1D("hMean276", "Mean jet mass vs #it{p}_{T} #sqrt{s} = 2.76 TeV; #it{p}_{T} (GeV/#it{c}); #LT #it{M} #GT (GeV/#it{c}^{2})", nptbins, ptlims);
	
	//histogram syst error
	Int_t div = 3;
	const Int_t fakenptbins = nptbins*div;
	
	
	TH1D* hDiffsqrts = new TH1D("hDiffsqrts", "Difference mean jet mass #sqrt{s} = 5.02 TeV - #sqrt{s} = 2.76 TeV vs #it{p}_{T} ; #it{p}_{T} (GeV/#it{c}); #LT #it{M} #GT (GeV/#it{c}^{2})", fakenptbins, ptlims[0], ptlims[nptbins]);
	hDiffsqrts->SetFillStyle(1001);
	hDiffsqrts->SetFillColor(kGray+2);
	hDiffsqrts->SetLineColor(hDiffsqrts->GetFillColor());
	hDiffsqrts->SetMarkerColor(hDiffsqrts->GetFillColor());
	//calculate the average
	for(Int_t ipt = 0; ipt< nptbins; ipt++){
		Int_t fakebin = hDiffsqrts->FindBin((ptlims[ipt+1] + ptlims[ipt])/2.);
		TH1D* htmpM502 = (TH1D*)list502->At(ipt);
		TH1D* htmpM276 = (TH1D*)list276->At(ipt);
		Double_t jetmdiff = htmpM502->GetMean()-htmpM276->GetMean();
		
		hDiffsqrts->SetBinContent(fakebin, hMeanpPb->GetBinContent(ipt+1)-jetmdiff/2.);
		hDiffsqrts->SetBinError(fakebin, jetmdiff/2.);
	
	}
	
	return hDiffsqrts;
}

//_____________________________________________________________________________
TList* TreatSystematicsPbPb(Int_t* reject, TString name){
	TString inputfile = "/data/Work/jets/JetMass/PbPbResults/DistJetMassPriorVar_Area.root";
	TString results = "/data/Work/jets/JetMass/PbPbResults/DistJetMassAllSyst_Area.root";
	
	TList *lout = new TList();
	lout->SetOwner(kTRUE);
	
	
	TFile* fsyst = new TFile(inputfile);
	if(!fsyst->IsOpen()){
		Printf("File %s not found", inputfile.Data());
		return 0x0;
	}
	fsyst->ls();
	
	TString nRelDif = "grUnfMassRelDiff";
	
	const Int_t nsyst = 4;
	TString sufRelD[nsyst] = {"1_PtBin", "2_PtBin", "3_PtBin", "Iter_PtBin"};
	
	TString lgndtxt[nsyst] = {"Prior", "Background", "TrackEff", "Iterations"};
	
	TFile* fRes = new TFile(results);
	if(!fRes->IsOpen()){
		Printf("File %s not found", results.Data());
		return 0x0;
	}
	
	TString nameSys = "grUnfMassSyst_PtBin";
	TString nameUnf = "grUnfMass_PtBin";
	
	Int_t nbinspt = 4;
	Int_t nx, ny, dx, dy;
	CalculatePads(nbinspt, nx, ny, dx, dy);
	
	TCanvas *cRelDiff = new TCanvas("cRelDiff", "Raltive difference PbPb systematics", dx, dy);
	cRelDiff->Divide(nx, ny);
	
	
	TLegend *legRelD = new TLegend(0.2, 0.6, 0.5, 0.8);
	legRelD->SetBorderSize(0);
	legRelD->SetFillStyle(0);
	
	
	
	TH1D* hNewSum[nbinspt];
	
	
	for(Int_t ipt = 0; ipt<nbinspt; ipt++){
		
		TH1D* hRelDiff[nsyst];
		Printf("Running pT bin %d/%d", ipt, nbinspt);
		TGraphErrors *gUnfMass = (TGraphErrors*)fRes->Get(TString::Format("%s%d", nameUnf.Data(),  ipt+1));
		gUnfMass->SetName(TString::Format("%s%d", nameUnf.Data(),  ipt+1));
		Printf("Info graph name %s, npoints = %d", gUnfMass->GetName(), gUnfMass->GetN());
		if(!gUnfMass) {
			Printf("--------------------%s%d not found", nameUnf.Data(), ipt+1);
			fsyst->ls();
			continue;
		}
		TGraphErrors *gSystAll = (TGraphErrors*)fRes->Get(TString::Format("%s%d", nameSys.Data(),  ipt+1));
		gSystAll->SetName(TString::Format("%s%d", nameSys.Data(),  ipt+1));
		
		for(Int_t isy = 0 ; isy < nsyst; isy++){
			if(reject[isy]) {
				hRelDiff[isy] = 0x0; 
				continue;
			}
			TGraphErrors *gRelDif = (TGraphErrors*)fsyst->Get(TString::Format("%s%s%d", nRelDif.Data(), sufRelD[isy].Data(), ipt+1));
			gRelDif->SetName(TString::Format("%s%s%d", nRelDif.Data(), sufRelD[isy].Data(), ipt+1));
			
			if(!gRelDif) continue;
			
			
			gRelDif->SetLineColor(colors[isy]);
			gRelDif->SetMarkerColor(colors[isy]);
			if(ipt == 0) legRelD->AddEntry(gRelDif, lgndtxt[isy], "lP");
			
			
			cRelDiff->cd(ipt+1);
			if(isy == 0) gRelDif->Draw("PA");
			else gRelDif->Draw("P");
			if(isy == nsyst-1) legRelD->Draw();
			Int_t nbinmass = gRelDif->GetN();
			Double_t mass, reldiff, error;
			gRelDif->GetPoint(0, mass, reldiff);
			gRelDif->GetPoint(1,  error, reldiff);
			//Printf("Point 0, mass = %f, error = %f", mass, error);
			error = (error - mass)*0.5;
			Double_t minh, maxh;
			minh = mass - error;
			gRelDif->GetPoint(nbinmass-1, mass, reldiff);
			gRelDif->GetPoint(nbinmass-2, error, reldiff);
			//Printf("Point %d, mass = %f, error = %f", nbinmass-2, mass, error);
			error = (mass - error)*0.5;
			maxh = mass + error;
			Printf("%d points, min = %f, max = %f", nbinmass, minh, maxh);
			
			hRelDiff[isy] = new TH1D(TString::Format("hRelDiff_sys%d_PtBin%d", isy, ipt+1), TString::Format("hRelDiff_sys%d_PtBin%d", isy, ipt+1), nbinmass, minh, maxh);
			hRelDiff[isy]->SetLineColor(gRelDif->GetLineColor());
			hRelDiff[isy]->SetMarkerColor(gRelDif->GetMarkerColor());
			hRelDiff[isy]->SetMarkerStyle(gRelDif->GetMarkerStyle());
			hRelDiff[isy]->GetYaxis()->SetRangeUser(0, 2);
			hRelDiff[isy]->Sumw2();
			for(Int_t ib = 0; ib<nbinmass; ib++){
				
				gRelDif->GetPoint(ib, mass, error);
				hRelDiff[isy]->SetBinContent(ib+1, 1);
				hRelDiff[isy]->SetBinError(ib+1, error);
			
			}
			cRelDiff->cd(ipt+1);
			hRelDiff[isy]->Draw("sames");
		}
		
		
		hNewSum[ipt] = AddInQuadrature(hRelDiff, nsyst,-1, TString::Format("hNewSum_PtBin%d", ipt+1));
		hNewSum[ipt]->SetMarkerStyle(24);
		hNewSum[ipt]->SetMarkerColor(kMagenta+10);
		hNewSum[ipt]->SetLineColor(kMagenta+10);
		if(ipt == 0) legRelD->AddEntry(hNewSum[ipt], "SumQ", "lP");
		cRelDiff->cd(ipt+1);
		hNewSum[ipt]->Draw("sames");
		
		TGraphErrors* gNewSyst = new TGraphErrors(hNewSum[ipt]);
		gNewSyst->SetName(TString::Format("gSyst%s_PtBin%d", name.Data(), ipt+1));
		lout->Add(gNewSyst);
		lout->Add(hNewSum[ipt]);
		lout->Add(gUnfMass);
		lout->Add(gSystAll);
		
	}
	
	Printf("List contains:");
	lout->ls();
	return lout;
}

//_____________________________________________________________________________

TH1D* GetpPb(TFile* finpPb, Int_t index){
	const char* nameMpPb = "hUnfM_Itr3_Pt";
	
	if(!finpPb->IsOpen()){
		Printf("GetpPb::File is not open");
		return 0x0;
	}
	TString name = TString::Format("%s%.0f_%.0f", nameMpPb, ptlims[index], ptlims[index+1]);
	Printf("Looking for %s", name.Data());
	TH1D* hMasspPb = (TH1D*)finpPb->Get(name)->Clone();
	
	if(!hMasspPb){
		Printf("%s%.0f_%.0f not found", nameMpPb, ptlims[index], ptlims[index+1]);
		finpPb->ls();
		
	}
	//debug
	Printf("Returning %p, named %s%.0f_%.0f", hMasspPb, nameMpPb, ptlims[index], ptlims[index+1]);
	return hMasspPb;
}

//_____________________________________________________________________________

void GetpPbSystematics(Int_t* reject, TString name, TFile* finpPb, Int_t index,  TH1D*& hNewSyspPb, TH1D*& hFullSyspPb){
	
	if(!finpPb->IsOpen()){
		Printf("GetpPbSystematics::File is not open");
		return ;
	}
	const Int_t totsys = 6;
	TString nameSyspPb[totsys] = {"hSyst0BkgSub_Pt", "hRespRangesmoothSysOnErr",  "hSystMxIter_Pt", "hOvlExclusionsmoothSysOnErr", "hTrkEffsmoothSysOnErr", "hSyst1Prior_Pt"};
	Int_t usefloat[totsys] = {1, 0, 2, 0, 0, 1};
	
	TH1D* hSys[totsys];
	
	for(Int_t isy = 0; isy<totsys; isy++){
		if(reject[isy]) {
			hSys[isy] = 0x0;
			continue;
		}
		TString nameh = nameSyspPb[isy];
		if(!usefloat[isy]) nameh+=index;
		if(usefloat[isy]>0) nameh+=TString::Format("%.0f_%.0f", ptlims[index], ptlims[index+1]);
		if(usefloat[isy]==2) nameh+=TString::Format("_Pt%.0f_%.0f", ptlims[index], ptlims[index+1]);
		
		Printf("Reading %s", nameh.Data());
		
		hSys[isy] = (TH1D*)finpPb->Get(nameh);
		if(!hSys[isy]) Printf("Problem, %s not found", nameh.Data());
		//else Printf("Found!!!! file is %d", finpPb->IsOpen());
		hSys[isy]->SetName(TString::Format("hSystSource%d_Pt%.0f_%.0f", isy, ptlims[index], ptlims[index+1]));
		
	}
	
	hNewSyspPb = (TH1D*)AddInQuadrature(hSys, totsys, index, TString::Format("hSyst%s_Pt%.0f_%.0f", name.Data(), ptlims[index], ptlims[index+1]))->Clone();
	
	TString nameTotSyst = TString::Format("hSystTot_Pt%.0f_%.0f", ptlims[index], ptlims[index+1]);
	hFullSyspPb = (TH1D*)finpPb->Get(nameTotSyst)->Clone();
	
	if(!hFullSyspPb) Printf("Total syst %s not found", nameTotSyst.Data());
	
	
	
}
//_____________________________________________________________________________
void RunRatioReducedSyst(){
	
	TString name = "All";

	//for the order see legend in method above ("Prior", "Background", "TrackEff", "Iterations")
	Int_t nsysttot = 4;
	Int_t nsysttotpPb = 6;
	Int_t rejectPbPb[nsysttot] = {0, 0, 0, 0};
	
	//pPb: bkg sub, range, iter, ovl exlu, tr eff, prior
	Int_t rejectpPb[nsysttotpPb] = {0, 0, 0, 0, 0, 0};
	
	RatioWithReducedSyst(rejectPbPb, rejectpPb, name);
	
	name = "NoTrEff";
	rejectPbPb[2] = 1;
	rejectpPb[4] = 1; 
	RatioWithReducedSyst(rejectPbPb, rejectpPb, name);
	
	name = "NoPrior";
	
	rejectPbPb[3] = 1; rejectPbPb[2] = 0;
	rejectpPb[5] = 1; rejectpPb[4] = 0; 
	
	RatioWithReducedSyst(rejectPbPb, rejectpPb, name);
	
	name = "NoPriorNoTrEff";
	
	rejectPbPb[3] = 1; rejectPbPb[2] = 1;
	rejectpPb[5] = 1; rejectpPb[4] = 1; 
	
	RatioWithReducedSyst(rejectPbPb, rejectpPb, name);
}

//_____________________________________________________________________________

void RatioWithReducedSyst(Int_t* rejectPbPb, Int_t* rejectpPb, TString name){
	
	
	
	TList* listSystPbPb = TreatSystematicsPbPb(rejectPbPb, name);
	Int_t offsetPbPb = 2;
	
	TString hnewsysn = TString::Format("gSyst%s_PtBin", name.Data());
	TString nameUnf = "grUnfMass_PtBin";
	TString nameSys = "grUnfMassSyst_PtBin";
	
	Int_t nbinspt = 3;
	Int_t nx, ny, dx, dy;
	CalculatePads(nbinspt, nx, ny, dx, dy);
	
	TCanvas *cMassPbPb = new TCanvas(TString::Format("cMassPbPb_Cmp%s", name.Data()), TString::Format("cMassPbPb_Cmp%s", name.Data()), dx, dy);
	cMassPbPb->Divide(nx, ny);
	
	TCanvas *cMasspPb = new TCanvas(TString::Format("cMasspPb_Cmp%s", name.Data()), TString::Format("cMasspPb_Cmp%s", name.Data()), dx, dy);
	cMasspPb->Divide(nx, ny);
	
	TCanvas *cMassRatio = new TCanvas(TString::Format("ccMassRatio_Cmp%s", name.Data()), TString::Format("ccMassRatio_Cmp%s", name.Data()), dx, dy);
	cMassRatio->Divide(nx, ny);
	
	TH1D *hRatio[nbinspt];
	TH1D *hRatioSys[nbinspt];
	//read the pPb results
	TString filenameSyspPb = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystRatio/TotalSystematicUnc.root";
	TString filenamepPb = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystRatio/MasspPbResults.root";
	TFile *finpPb = new TFile(filenamepPb);
	if(!finpPb->IsOpen()){
		Printf("File %s is not open", filenamepPb.Data());
		return;
	}
	finpPb->ls();
	TFile *finpPbSy = new TFile(filenameSyspPb);
	if(!finpPbSy->IsOpen()){
		Printf("File %s is not open", filenameSyspPb.Data());
		return;
	}
	finpPbSy->ls();
	
	
	//read final results where the ratio is
	TString nameRatioSysDef = "hRatioPbPbOpPbSys_Pt", namePbPbSys = "hUnfMassSyst_PtBin";
	TString filenameFinal = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystRatio/FinalResults.root";
	TFile *fFinal = new TFile(filenameFinal);
	if(!fFinal->IsOpen()){
		Printf("file %s not found", filenameFinal.Data());
		return;
	}
	
	TLegend *legMPbPb = new TLegend(0.5, 0.5, 0.8, 0.8);
	legMPbPb->SetFillStyle(0);
	legMPbPb->SetBorderSize(0);
	TLegend *legMpPb = new TLegend(0.5, 0.5, 0.8, 0.8);
	legMpPb->SetFillStyle(0);
	legMpPb->SetBorderSize(0);
	
	TLegend *legRatio = new TLegend(0.5, 0.5, 0.8, 0.8);
	legRatio->SetFillStyle(0);
	legRatio->SetBorderSize(0);
	
	for(Int_t ipt = 0; ipt < nbinspt; ipt++){
		
		TGraphErrors *gMassPbPb = (TGraphErrors*)listSystPbPb->FindObject(TString::Format("%s%d", nameUnf.Data(), ipt+offsetPbPb));
		if(!gMassPbPb) {
			Printf("%s%d not found", nameUnf.Data(), ipt+offsetPbPb);
			listSystPbPb->ls();
			continue;
		}
		gMassPbPb->SetName(TString::Format("%s%d", nameUnf.Data(), ipt+offsetPbPb));
		
		TGraphErrors *gMassSysO = (TGraphErrors*)listSystPbPb->FindObject(TString::Format("%s%d", nameSys.Data(), ipt+offsetPbPb));
		if(!gMassSysO) {
			Printf("%s%d not found", nameSys.Data(), ipt+offsetPbPb);
			listSystPbPb->ls();
			continue;
		}
		gMassSysO->SetName(TString::Format("%s%d", nameSys.Data(), ipt+offsetPbPb));
		gMassSysO->SetFillStyle(0);
		gMassSysO->SetLineWidth(2);
		
		TGraphErrors *gSysNew = (TGraphErrors*)listSystPbPb->FindObject(TString::Format("%s%d", hnewsysn.Data(), ipt+offsetPbPb));
		if(!gSysNew) {
			Printf("%s%d not found", hnewsysn.Data(), ipt+offsetPbPb);
			listSystPbPb->ls();
			continue;
		}
		gSysNew->SetName(TString::Format("%s%d", hnewsysn.Data(), ipt+offsetPbPb));
		
		TH1D *hMassPbPb = 0x0;
		TH1D *hSystrs   = TGraphToTH1D(gSysNew, TString::Format("hSystrs_PtBin%d", ipt+offsetPbPb));
		TH1D *hMasstrs  = TGraphToTH1D(gMassPbPb, TString::Format("hMasstrs_PtBin%d", ipt+offsetPbPb));
		
		TH1D *hSysNew = MinimumRange(hSystrs, hMasstrs, hMassPbPb);
		hSysNew->  SetName(TString::Format("hSysNew_PtBin%d", ipt+offsetPbPb));
		hMassPbPb->SetName(TString::Format("hMassPbPb_PtBin%d", ipt+offsetPbPb));
		
		Printf("New histograms %p, %p", hSysNew, hMassPbPb);
		hMassPbPb->SetMarkerStyle(gMassPbPb->GetMarkerStyle());
		hMassPbPb->SetMarkerColor(gMassPbPb->GetMarkerColor());
		hSysNew->SetFillStyle(1001);
		hSysNew->SetFillColor(kCyan);
		
		SetMassValueInSystematic(hSysNew, hMassPbPb);
		
		cMassPbPb->cd(ipt+1);
		//gMassPbPb->Draw("AP");
		hSysNew->GetYaxis()->SetRangeUser(0, 0.25);
		hSysNew->Draw("E2");
		hMassPbPb->Draw("sames");
		gMassSysO->Draw("E2sames");
		if(ipt == 0){
			legMPbPb->AddEntry(hSysNew, TString::Format("Syst %s", name.Data()));
			legMPbPb->AddEntry(gMassSysO, "Original Syst");
			
		}
		legMPbPb->Draw();
		
		//pPb
		
		TH1D* hMasspPb = GetpPb(finpPb, ipt);
		
		TH1D* hSyspPbNew = 0x0;
		TH1D* hSyspPb    = 0x0;
		GetpPbSystematics(rejectpPb, name, finpPbSy, ipt,  hSyspPbNew, hSyspPb);
		hSyspPbNew->SetFillStyle(0);
		
		if(!hMasspPb) Printf("Problem, ipt = %d hMasspPb %p", ipt, hMasspPb);
		if(!hSyspPb)  Printf("Problem, ipt = %d hSyspPb %p", ipt, hSyspPb);
		TH1D *hh = 0x0;
		TH1D *hhs = 0x0;
		
		//relative uncertainties, can uncomment and plot for debugging
		//cMasspPb->cd(ipt+1);
		//hSyspPbNew->Draw("E2");
		//hSyspPb->Draw("E2sames");
		SetMassValueInSystematic(hSyspPbNew, hMasspPb);
		SetMassValueInSystematic(hSyspPb, hMasspPb);
		hSyspPbNew->GetYaxis()->SetRangeUser(0., 0.25);
		hSyspPbNew->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		cMasspPb->cd(ipt+1);
		hSyspPbNew->Draw("E2");
		hSyspPb->Draw("E2sames");
		hMasspPb->Draw("sames");
		
		if(ipt == 0){
			legMpPb->AddEntry(hSyspPbNew, TString::Format("Syst %s", name.Data()));
			legMpPb->AddEntry(hSyspPb, "Original Syst");
			
		}
		legMpPb->Draw();
		
		// reduce the range to the minimum one
		hRatio[ipt] = MinimumRange(hMassPbPb, hMasspPb, hh);
		hRatioSys[ipt] = MinimumRange(hSysNew, hSyspPbNew, hhs);
		hRatio[ipt]->SetName(TString::Format("hMassPbPbOpPb_PtBin%d", ipt+1));
		hRatioSys[ipt]->SetName(TString::Format("h%sSysPbPbOpPb_PtBin%d", name.Data(), ipt+1));
		hRatio[ipt]->GetYaxis()->SetTitle("#Rgothic_{#sqrt{#it{s}}}");
		hRatioSys[ipt]->GetYaxis()->SetTitle("#Rgothic_{#sqrt{#it{s}}}");
		
		hRatio[ipt]->Divide(hh);
		hRatioSys[ipt]->Divide(hhs);
		hRatioSys[ipt]->SetFillStyle(0);
		hRatioSys[ipt]->SetFillColor(kRed);
		hRatioSys[ipt]->SetLineColor(kRed);
		hRatioSys[ipt]->SetLineWidth(2);
		
		TF1 *fpol3RatioDef = new TF1(TString::Format("fpol3RatioDef_Pt%d", ipt), "[0] + [1]*x + [2]*x*x + [3]*x*x", 0., hRatioSys[ipt]->GetBinLowEdge(hRatioSys[ipt]->GetNbinsX()+1));
		fpol3RatioDef->SetParameters(1, 1e-2, 1e-3, 1e-3);
		fpol3RatioDef->SetLineColor(hRatioSys[ipt]->GetLineColor());
		
		TH1D* hSysRatioDef = (TH1D*)fFinal->Get(TString::Format("%s%.0f_%.0f", nameRatioSysDef.Data(), ptlims[ipt], ptlims[ipt+1]));
		
		hSysRatioDef->Fit(fpol3RatioDef, "RL+0");
		TH1D* hSysPbPbRes = (TH1D*)fFinal->Get(TString::Format("%s%d", namePbPbSys.Data(), ipt+offsetPbPb)); //watch out the offset!!!!
		TH1D* hRatioPyt   = (TH1D*)fFinal->Get(TString::Format("hRatio276O502M_Pt%d", ipt));
		
		cMassRatio->cd(ipt+1);
		hSysRatioDef->GetYaxis()->SetRangeUser(0., 4.);
		hSysRatioDef->Draw("E2");
		hRatio[ipt]->Draw("sames");
		hRatioSys[ipt]->Draw("E2sames");
		fpol3RatioDef->Draw("sames");
		hRatioPyt->Draw("sames");
		
		if(ipt == 0){
			legRatio->AddEntry(hRatioSys[ipt], TString::Format("Syst %s", name.Data()));
			legRatio->AddEntry(hSysRatioDef, "Original Ratio");
			
		}
		legRatio->Draw();
		// results preliminary, uncomment to debug
		//cMassPbPb->cd(ipt+1);
		//hSysPbPbRes->Draw("E2sames");
		
	}
	SaveCv(cMassPbPb);
	SaveCv(cMasspPb );
	SaveCv(cMassRatio);
}


#endif
