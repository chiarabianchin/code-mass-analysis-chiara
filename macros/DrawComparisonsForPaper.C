#ifndef DrawComparisonsForPaper_C
#define DrawComparisonsForPaper_C
#include "/data/macros/LoadALICEFigures.C"
#include "/data/Work/MyCodeJetMass/utils/CommonTools.C"
#include "/data/Work/MyCodeJetMass/unfold/SystematicComparisonsUnf.C"
#include <TLegend.h>
//#include <TGaxis.h>

//global

Double_t maxM = 26.;
//definitions
TList* GetPbPbResults();
TList* GetPbPbResultsGraph();
TList* GetpPbResults(Bool_t kinecorr = kTRUE);
TList* GetpPbMarta();
TList* GetPythiapPb(Int_t input = 0);
TList* GetPythiapPbMarta();
TList* GetPythiaPbPb();
TList* GetJEWELPbPb();
TList* GetJEWELPbPbrecoilOff();
TList* GetJEWELpp();
TList* GetJEWEL(Int_t system);
TList* GetJEWELQPYTHIA(Int_t mod);


TH1D* CalculateMean(Double_t **xrange, Int_t nbins, Double_t lims[], TH1D** hdistr, TGraphErrors *&grmean, TH1D *&hrms, const char* namehmean);
void DrawMeanComparison(Double_t **xrange, const Int_t ninputs, TString inputDistr[], TString h1names[], const Int_t ninputsSy, TString inputSyst[], TString hsysnames[], Int_t offset[], TString leg[], Int_t mrk[], Int_t clr[], Int_t clrsy[], Int_t fillsy[]);
TH1D* CalculateSysMean(Double_t **xrange, Int_t nbins, Double_t lims[], TH1D** hSysdistr, TGraphErrors *&grSysmean, TH1D *&hSysrms, const char* namehsysmean);

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
		hMassSy->SetMarkerStyle(1);
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
		
		TGraphErrors* gMassSt = (TGraphErrors*)fResPbPbM->Get(Form("%s%d", nameMSt.Data(), ih+offset));
		TGraphErrors* gMassSy = (TGraphErrors*)fResPbPbM->Get(Form("%s%d", nameMSy.Data(),ih+offset));
		
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
		
		TH1D* hMassSt = new TH1D(Form("hUnfMass_PtBin%d", ih+offset), Form("hUnfMass_PtBin%d", ih+offset), npoints, binlimits);
		hMassSt->Sumw2();
		TH1D* hMassSy = new TH1D(Form("hUnfMassSyst_PtBin%d", ih+offset), Form("hUnfMassSyst_PtBin%d", ih+offset), npoints, binlimits);
		hMassSy->Sumw2();
		// set this nicely at some point
		hMassSt->SetMarkerColor(kBlue+2);
		hMassSt->SetMarkerStyle(20);
		hMassSt->SetLineColor(kBlue+2);
		
		hMassSy->SetFillColor(kBlue-10);
		hMassSy->SetMarkerStyle(1);
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
	const Int_t nhM = fin->GetListOfKeys()->GetEntries()/2;
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
		
		TH1D* hMassSt = new TH1D(Form("%sMax%.0f_%.0f", nameMSt.Data(), ptlims[ih],ptlims[ih+1]), Form("Mass pPb Pt %.0f - %.0f; #it{M} (GeV/c^2)", ptlims[ih],ptlims[ih+1]), nbinsM, 0, maxRangeMassFinal[ih]);
		hMassSt->Sumw2();
		hMassSt->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		
		TH1D* hMassSy = (TH1D*)hMassSt->Clone(Form("%sMax%.0f_%.0f", nameMSy.Data(), ptlims[ih],ptlims[ih+1]));
		hMassSy->Sumw2();
		
		TH1D* hMtmp = (TH1D*)fin->Get(Form("%s%.0f_%.0f", nameMSt.Data(), ptlims[ih],ptlims[ih+1]));
		
		TH1D* hSytmp = 0x0;
		if(kinecorr) hSytmp = (TH1D*)listKineEffSy->At(ih*nhM+1);
		else hSytmp = (TH1D*)fin->Get(Form("%s%.0f_%.0fb", nameMSy.Data(), ptlims[ih],ptlims[ih+1]));
		Printf("Low edge %f and %f (they have to be the same)", hMtmp->GetBinLowEdge(0), hSytmp->GetBinLowEdge(0));
		
		for(Int_t ib = 0; ib < nbinsM; ib++){
			Int_t bin = hMtmp->FindBin(hMassSt->GetBinCenter(ib+1));
			//Printf("BIN = %d", bin);
			if(bin == 0 || bin == -1) {
				//Printf("WAAARNIIIINGGGGGG");
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
		//hMassSt->GetYaxis()->SetTitleOffset(1.7);
		
		//hMassSt->GetXaxis()->SetRangeUser(0., maxRangeMassFinal[ih]);
		hMassSy->SetFillStyle(hSytmp->GetFillStyle());
		hMassSy->SetLineWidth(hSytmp->GetLineWidth());
		hMassSy->SetLineColor(hSytmp->GetLineColor());
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
		TH1D* hMassSt = (TH1D*)fRespPbM->Get(Form("%s%d", nameMSt.Data(), ih+1));
		TH1D* hMassSy = (TH1D*)fRespPbM->Get(Form("%s%d", nameMSy.Data(),ih+1));
		
		// set this nicely at some point
		hMassSt->SetLineColor(kOrange+1);
		hMassSt->SetMarkerColor(kOrange+1);
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
		hMassPythia[ih]->SetMarkerColor(kGreen-3);
		hMassPythia[ih]->SetLineColor(kGreen-3);
		listPythia->Add(hMassPythia[ih]);
	}
	return listPythia;
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
		hMassSt->SetMarkerColor(kBlue-7);
		hMassSt->SetLineColor(hMassSt->GetMarkerColor());
		hMassSt->SetTitle("; #it{M}_{ch jet} (GeV/#it{c}^{2}); #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{M}_{ch jet}} (#it{c}^{2}/GeV)");
		
		listPythia->Add(hMassSt);
	}
	
	return listPythia;
}

//________________________________________________________________________

TList* GetPythiapPbMarta(){
	
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
		TH1D* hMassSt = (TH1D*)fMPythia502->Get(Form("%s%d", nameMSt.Data(), ih+offset));
		
		// set this nicely at some point
		hMassSt->SetMarkerStyle(25);
		hMassSt->SetMarkerColor(kBlack);
		hMassSt->SetLineColor(kBlack);
		
		listPythia->Add(hMassSt);
	}
	
	return listPythia;
}

//________________________________________________________________________

TList* GetJEWELPbPb(){
	return GetJEWEL(1);
}

//________________________________________________________________________

TList* GetJEWELPbPbrecoilOff(){
	return GetJEWEL(2);
}

//________________________________________________________________________

TList* GetJEWELpp(){
	return GetJEWEL(0);
}

//________________________________________________________________________
TList* GetJEWEL(Int_t system){
	//system: 0 = pp, 1 = 0-10% Pb-Pb, 2 = 0-10% recoils off
	
	TString pathJewel = "/data/Work/jets/JetMass/PbPbJEWEL/alice_jetmass_histograms.root";
	if(system == 2) pathJewel =	"/data/Work/jets/JetMass/PbPbJEWEL/alice_jetmass_norecoils_histograms.root";
	TFile *fJewel = new TFile(pathJewel);
	if(!fJewel->IsOpen()){
		Printf("File %s not found", pathJewel.Data());
		return 0;
	}
	Int_t cent = 0;
	if(system > 0) cent = 1;
	TString nameMSt = Form("Mass_centbin%d_ptbin", cent);
	TString nameRat = "Ratio_centbin1_ptbin";
	const Int_t nhM = 3;
	Int_t offset = 0;
	
	TList *listJewel = new TList();
	listJewel->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fJewel->Get(Form("%s%d", nameMSt.Data(), ih+offset));
		TH1D* hMassRa = (TH1D*)fJewel->Get(Form("%s%d", nameRat.Data(), ih+offset));
		
		// set this nicely at some point
		hMassSt->SetMarkerStyle(28);
		if(system == 0){
			
			hMassSt->SetMarkerColor(kGreen+2);
		} 
		if(system == 1) {
			hMassSt->SetMarkerColor(kOrange+1);
		}
		if(system == 2) {
			hMassSt->SetMarkerColor(kRed);
		}
		hMassSt->SetLineColor(hMassSt->GetMarkerColor());
		hMassRa->SetMarkerStyle(hMassSt->GetMarkerStyle());
		hMassRa->SetMarkerColor(hMassSt->GetMarkerColor());
		hMassRa->SetLineColor(hMassSt->GetMarkerColor());
		
		listJewel->Add(hMassSt);
		listJewel->Add(hMassRa);
	}
	
	return listJewel;

}

//________________________________________________________________________

TList* GetJEWELQPYTHIA(Int_t mod){
	
	//mod: 1 = JEWEL recoil On, 2 = JEWEL recoil Off, 3 = QPYTHIA, 4 = Vacuum
	
	TString pathModels = "/data/Work/jets/JetMass/PbPbJEWEL/Marco/JetMassJEWELR040.root";
	
	TFile *fModels = new TFile(pathModels);
	if(!fModels->IsOpen()){
		Printf("File %s not found", pathModels.Data());
		return 0;
	}

	TString nameMSt = "";
	Int_t marker = 28;
	Int_t color  = kOrange+1;
	
	if(mod == 1) nameMSt = "hM_JEWEL RecOn_";
	if(mod == 2) {
		nameMSt = "hM_JEWEL RecOff_";
		marker = 28;
		color  = kRed;
	}
	if(mod == 3) {
		nameMSt = "hM_QPYTHIA_";
		marker = 33;
		color  = kRed+2;
	}
	if(mod == 4) {
		nameMSt = "hM_Vacuum_";
		marker = 21;
		color  = kViolet+5;
	}
	
	const Int_t nhM = 5;
	Int_t offset = 0;
	
	TList *listModels = new TList();
	listModels->SetOwner();
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		TH1D* hMassSt = (TH1D*)fModels->Get(Form("%s%d", nameMSt.Data(), ih+offset));
		
		hMassSt->SetMarkerStyle(marker);
		hMassSt->SetMarkerColor(color);
		hMassSt->SetLineColor(hMassSt->GetMarkerColor());
		
		listModels->Add(hMassSt);
	}
	
	return listModels;

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
void RatiopPbPbPb(Bool_t stylezero = kFALSE, Bool_t pPbpaperprop = kFALSE, Bool_t corrkine = kFALSE, Bool_t show1134 = kFALSE){
	TString suff = "";
	if(!corrkine) suff = "NoKineCor";
	//TGaxis::SetMaxDigits(2);
	if(stylezero) {
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		gStyle->SetTextFont(42);
		
		gStyle->SetLabelSize(0.05, "X");
		gStyle->SetTitleOffset(1.1, "X");
		gStyle->SetTitleSize(0.048, "X");
		
		gStyle->SetLabelSize(0.05, "Y");
		gStyle->SetTitleOffset(1.34, "Y");
		gStyle->SetTitleSize(0.048, "Y");
		
		gStyle->SetNdivisions(1004, "Y");
	}
	
	gStyle->SetPadBottomMargin(.13);
	gStyle->SetPadLeftMargin(.16);
	gStyle->SetPadRightMargin(.06);
	
	// text for the final plots
	TPaveText *pvgeneral = new TPaveText(0.58, 0.5, 0.83, 0.6, "NDC");
	pvgeneral->SetFillStyle(0);
	pvgeneral->SetBorderSize(0);
	pvgeneral->AddText("Anti-#it{k}_{T}, #it{R} = 0.4");
	
	
	TPaveText *pvSystPbPb = new TPaveText(0.22, 0.7, 0.78, 0.8, "NDC");
	pvSystPbPb->SetFillStyle(0);
	pvSystPbPb->SetBorderSize(0);
	pvSystPbPb->AddText("Pb-Pb 0-10% #sqrt{#it{s}}_{NN} = 2.76 TeV");
	TPaveText *pvSystpPb = new TPaveText(0.35, 0.6, 0.77, 0.7, "NDC");
	pvSystpPb->SetFillStyle(0);
	pvSystpPb->SetBorderSize(0);
	pvSystpPb->AddText("p-Pb #sqrt{#it{s}}_{NN} = 5.02 TeV");
	TPaveText *pvSystPyth = new TPaveText(0.35, 0.6, 0.77, 0.7, "NDC");
	pvSystPyth->SetFillStyle(0);
	pvSystPyth->SetBorderSize(0);
	pvSystPyth->AddText("PYTHIA Perugia 2011");
	
	//get the PbPb results
	TList* listPbPb = GetPbPbResultsGraph(); //GetPbPbResults();
		
	//get the pPb results
	TList* listpPb = GetpPbResults(corrkine);
	TList* listpPbM = GetpPbMarta();
	
	//get Pythia pPb
	TList* listPythia502 = GetPythiapPb();
	//get Pythia pPb
	TList* listPythia502M = GetPythiapPbMarta();
	
	//get Pythia PbPb
	TList* listPythia276 = GetPythiaPbPb();
	
	//get Jewel PbPb
	TList* listJewelPbPb = GetJEWELPbPb();
	
	TList* listJewelMRecoilOff = //GetJEWELPbPbrecoilOff(); 
	GetJEWELQPYTHIA(2);
	
	//get Jewel pPb
	TList* listJewelpp = GetJEWELpp();
	
	//get Q-PYTHIA
	TList* listqpythia = GetJEWELQPYTHIA(3);
	
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
	
	Int_t npPb    = listpPb ->GetEntries()/2;
	Int_t nPbPb   = listPbPb->GetEntries()/2;
	Int_t npPbM   = listpPbM->GetEntries()/2;
	Int_t nPy502  = listPythia502->GetEntries();
	Int_t nPy502M = listPythia502M->GetEntries();
	Int_t nPy276  = listPythia276->GetEntries();
	Int_t nJewpp  = listJewelpp ->GetEntries()/2;
	Int_t nJewPbPb= listJewelPbPb->GetEntries()/2;
	Int_t nJewPbPbMRoff = listJewelMRecoilOff->GetEntries();
	Int_t nqpythia= listqpythia->GetEntries();
	
	Int_t offsetPbPb   = nPbPb   - npPb;
	Int_t offsetpPbM   = npPbM   - npPb;
	Int_t offsetPy502  = nPy502  - npPb;
	Int_t offsetPy502M = nPy502M - npPb;
	Int_t offsetPy276  = nPy276  - npPb;
	Int_t offsetJewpp  = nJewpp  - npPb;
	Int_t offsetJewPbPb= nJewPbPb- npPb;
	Int_t offsetJewPbPbMRoff = nJewPbPbMRoff- npPb;
	Int_t offsetqpythia= nqpythia- npPb;
	
	Printf("offsets %d %d %d %d %d %d %d %d %d ", offsetPbPb, offsetpPbM , offsetPy502, offsetPy502M, offsetPy276, offsetJewpp, offsetJewPbPb, offsetJewPbPbMRoff, offsetqpythia);
	
	const Int_t nhM = npPb;
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nhM, nx, ny, dx, dy);
	TCanvas *cMass = new TCanvas(Form("cMass%s", suff.Data()), "Mass histograms", dx, dy);
	cMass->Divide(nx, ny);
	TLegend *legMass = new TLegend(0.3, 0.55, 0.85, 0.75);
	legMass->SetBorderSize(0);
	legMass->SetFillStyle(0);
	
	TCanvas *cRatioPbPbOpPb = new TCanvas(Form("cRatioPbPbOpPb%s", suff.Data()), "Ratio PbPb/pPb", dx, dy);
	cRatioPbPbOpPb->Divide(nx, ny);
	TLegend *legRatioPbPbOpPb = new TLegend(0.18, 0.55, 0.9, 0.8);
	legRatioPbPbOpPb->SetBorderSize(0);
	legRatioPbPbOpPb->SetFillStyle(0);
	
	TCanvas *cMasspPb = new TCanvas(Form("cMasspPb%s", suff.Data()), "Mass pPb compared with Pythia", dx, dy);
	cMasspPb->Divide(nx, ny);
	TLegend *legMasspPb = new TLegend(0.4, 0.5, 0.9, 0.8);
	legMasspPb->SetBorderSize(0);
	legMasspPb->SetFillStyle(0);
	
	TCanvas *cRelUnc = new TCanvas(Form("cRelUncSys%s", suff.Data()), "Relative systematic uncertainties mass spectra", dx, dy);
	cRelUnc->Divide(nx, ny);
	TLegend *legRelUnc = new TLegend(0.4, 0.5, 0.8, 0.9);
	legRelUnc->SetBorderSize(0);
	legRelUnc->SetFillStyle(0);
	
	TCanvas *cRatiopPbOPy = new TCanvas(Form("cRatiopPbOPy%s", suff.Data()), "Ratio pPb over Pythia 5.02 TeV", dx, dy);
	cRatiopPbOPy->Divide(nx, ny);
	TLegend *legRatiopPbOPy = new TLegend(0.22, 0.55, 0.9, 0.75);
	legRatiopPbOPy->SetBorderSize(0);
	legRatiopPbOPy->SetFillStyle(0);
	
	TCanvas *cMassPbPb = new TCanvas(Form("cMassPbPb%s", suff.Data()), "Mass PbPb compared with Pythia", dx, dy);
	cMassPbPb->Divide(nx, ny);
	TLegend *legMassPbPb = new TLegend(0.15, 0.5, 0.9, 0.8);
	legMassPbPb->SetBorderSize(0);
	legMassPbPb->SetFillStyle(0);
	
	TCanvas *cMassPbPbOnly = new TCanvas(Form("cMassPbPbOnly%s", suff.Data()), "Mass PbPb", dx, dy);
	cMassPbPbOnly->Divide(nx, ny);
	TLegend *legMassPbPbOnly = new TLegend(0.15, 0.55, 0.9, 0.75);
	legMassPbPbOnly->SetBorderSize(0);
	legMassPbPbOnly->SetFillStyle(0);
	
	TCanvas *cMassPbPbPy = new TCanvas(Form("cMassPbPbPy%s", suff.Data()), "Mass PbPb and Pythia", dx, dy);
	cMassPbPbPy->Divide(nx, ny);
	TLegend *legMassPbPbPy = new TLegend(0.15, 0.55, 0.9, 0.75);
	legMassPbPbPy->SetBorderSize(0);
	legMassPbPbPy->SetFillStyle(0);
	
	TCanvas *cMassPbPbModels = new TCanvas(Form("cMassPbPbModels%s", suff.Data()), "Mass PbPb in models", dx, dy);
	cMassPbPbModels->Divide(nx, ny);
	TLegend *legMassPbPbModels = new TLegend(0.15, 0.55, 0.9, 0.75);
	legMassPbPbModels->SetBorderSize(0);
	legMassPbPbModels->SetFillStyle(0);
	
	TCanvas *cRatioPbPbOPy = new TCanvas(Form("cRatioPbPbOPy%s", suff.Data()), "Ratio PbPb over Pythia 2.76 TeV", dx, dy);
	cRatioPbPbOPy->Divide(nx, ny);
	TLegend *legRatioPbPbOPy = new TLegend(0.22, 0.63, 0.9, 0.82);
	legRatioPbPbOPy->SetBorderSize(0);
	legRatioPbPbOPy->SetFillStyle(0);
	
	TCanvas *cRatioDataOPy = new TCanvas(Form("cRatioDataOPy%s", suff.Data()), "Ratio data over Pythia", dx, dy);
	cRatioDataOPy->Divide(nx, ny);
	TLegend *legRatioDataOPy = new TLegend(0.22, 0.58, 0.9, 0.78);
	legRatioDataOPy->SetBorderSize(0);
	legRatioDataOPy->SetFillStyle(0);
	
	Int_t n = 2;
	
	TPaveText **pvpt = new TPaveText*[nhM];
	
	for(Int_t ih = 0; ih<nhM ; ih++){
		maxM = maxRangeMassFinal[ih];
		pvpt[ih] = new TPaveText(0.33, 0.8, 0.83, 0.9, "NDC");
		pvpt[ih]->SetFillStyle(0);
		pvpt[ih]->SetBorderSize(0);
		pvpt[ih]->AddText(Form("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ih], ptlims[ih+1]));
		
		//mass
		TH1D* hMpPb     = (TH1D*)listpPb ->At(n*ih);
		TH1D* hMPbPb    = (TH1D*)listPbPb->At(n*(ih+offsetPbPb));
		TH1D* hMpPbM    = (TH1D*)listpPbM->At(n*(ih+offsetpPbM));
		                
		TH1D* hMPy502   = (TH1D*)listPythia502 ->At(ih+offsetPy502 );
		TH1D* hMPy502M  = (TH1D*)listPythia502M->At(ih+offsetPy502M);
		TH1D* hMPy276   = (TH1D*)listPythia276 ->At(ih+offsetPy276 );
		
		TH1D* hMJewpp   = (TH1D*)listJewelpp  ->At(n*(ih+offsetJewpp));
		TH1D* hMJewPbPb = (TH1D*)listJewelPbPb->At(n*(ih+offsetJewPbPb));
		
		TH1D* hMJewPbPbMRoff= (TH1D*)listJewelMRecoilOff ->At(ih+offsetJewPbPbMRoff);
		TH1D* hMQPy     = (TH1D*)listqpythia->At(ih+offsetqpythia);
		
		hMPy502->Scale(1./hMPy502->Integral("width"));
		hMPy502M->Scale(1./hMPy502M->Integral("width"));
		
		hMPy276->Scale(1./hMPy276->Integral("width"));
		
		//syst
		TH1D* hSypPb  = (TH1D*)listpPb ->At(n*ih+1);
		TH1D* hSypPbM = (TH1D*)listpPbM->At(n*(ih+offsetpPbM)+1);
		TH1D* hSyPbPb = (TH1D*)listPbPb->At(n*(ih+offsetPbPb)+1);
		//ratio Jewel
		TH1D* hPbPboppJewel = (TH1D*)listJewelPbPb->At(n*(ih+offsetJewPbPb)+1);
		
		//Draw comparison pPb, PbPb
		
		cMass->cd(ih+1);
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
		if(ih == 1) pvgeneral->Draw();
		pvpt[ih]->Draw();
		if(ih == nhM-1){
			legMass->AddEntry(hMPbPb , "Pb-Pb 0-10% #sqrt{#it{s}}_{NN} = 2.76 TeV", "PL");
			legMass->AddEntry(hSyPbPb, "Systematic Pb-Pb" , "F");
			legMass->AddEntry(hMpPb ,  "pPb #sqrt{#it{s}}_{NN} = 5.02 TeV", "PL");
			legMass->AddEntry(hSypPb,  "Systematic p-Pb" , "F");
			if(pPbpaperprop){
				legMass->AddEntry(hMpPbM , "Mass p-Pb paper prop", "PL");
				legMass->AddEntry(hSypPbM,  "Sys p-Pb paper prop" , "F");
			}
			legMass->Draw();
		}
		
		//relative uncertainties
		TH1D* hRelSysUncpPb  = (TH1D*)hMpPb->Clone(Form("hRelSysUncMpPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRelSysUncpPb->SetLineColor(kPink-9);
		hRelSysUncpPb->SetLineWidth(2);
		
		TH1D* hRelSysUncpPbM = (TH1D*)hMpPbM->Clone(Form("hRelSysUncMpPbM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		TH1D* hRelSysUncPbPb = (TH1D*)hMPbPb->Clone(Form("hRelSysUncMPbPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
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
		
		hSypPb ->GetYaxis()->SetRangeUser(0., 0.25);
		hSypPb ->Draw("E2");
		hMpPb  ->Draw("sames");
		if(show1134) hMPy502->Draw("sames");
		hMPy502M->Draw("sames");
		//hMJewpp-> Draw("sames");
		if(ih == nhM-1){
			legMasspPb->AddEntry(hMpPb, "p-Pb #sqrt{s_{NN}} = 5.02 TeV", "LP");
			legMasspPb->AddEntry(hSypPb, "Systematic p-Pb", "F");
			if(show1134) legMasspPb->AddEntry(hMPy502, "PYTHIA Perugia 2011 #sqrt{s_{NN}} = 5.02 TeV", "LP");
			legMasspPb->AddEntry(hMPy502M, "PYTHIA  Perugia 2011", "LP"); //  #sqrt{s_{NN}} = 5.02 TeV // paper prop
			//legMasspPb->AddEntry(hMJewpp, "JEWEL+PYTHIA pp", "LP");
			
			legMasspPb->Draw();
		}
		if(ih == 1) pvgeneral->Draw();
		pvpt[ih]->Draw();
		
		
		cMassPbPbOnly->cd(ih+1);
		hSyPbPb  ->Draw("E2");
		hMPbPb   ->Draw("sames");
		if(ih == 1) pvgeneral->Draw();
		pvpt[ih]->Draw();
		
		cMassPbPbPy->cd(ih+1);
		hSyPbPb  ->Draw("E2");
		hMPbPb   ->Draw("sames");
		hMPy276  ->Draw("sames");
		if(ih == 1) pvgeneral->Draw();
		pvpt[ih] ->Draw();
		//Draw comparison PbPb Pythia 2.76 TeV
		cMassPbPb->cd(ih+1);
		
		hSyPbPb  ->Draw("E2");
		hMPbPb   ->Draw("sames");
		hMPy276  ->Draw("sames");
		hMJewPbPb->Draw("sames");
		hMJewPbPbMRoff->Draw("sames");
		hMQPy    ->Draw("sames");
		pvpt[ih] ->Draw();
		
		cMassPbPbModels->cd(ih+1);
		hMPy276  ->GetYaxis()->SetNdivisions(1004, "Y");
		hMPy276  ->GetYaxis()->SetRangeUser(0, 0.25);
		hMPy276  ->Draw("sames");
		hMJewPbPb->Draw("sames");
		hMJewPbPbMRoff->Draw("sames");
		hMQPy    ->Draw("sames");
		pvpt[ih] ->Draw();
		
		if(ih == nhM-1){
			legMassPbPbOnly->AddEntry(hMPbPb,    "Pb-Pb 0-10% #sqrt{s_{NN}} = 2.76 TeV", "LP");
			legMassPbPbOnly->AddEntry(hSyPbPb,   "Systematic Pb-Pb" , "F");
			cMassPbPbOnly->cd(ih+1); legMassPbPbOnly->Draw();
			
			
			legMassPbPbModels->AddEntry(hMPy276,   "PYTHIA  Perugia 2011 #sqrt{s_{NN}} = 2.76 TeV", "LP");
			legMassPbPbModels->AddEntry(hMJewPbPb, "JEWEL+PYTHIA Recoil on 0-10% PbPb", "LP");
			legMassPbPbModels->AddEntry(hMJewPbPbMRoff, "JEWEL+PYTHIA Recoil off 0-10% PbPb", "LP");
			legMassPbPbModels->AddEntry(hMQPy,      "Q-PYTHIA", "LP");
			cMassPbPbModels->cd(ih+1); legMassPbPbModels->Draw();
			
			legMassPbPbPy->AddEntry(hMPbPb,    "Pb-Pb 0-10% #sqrt{s_{NN}} = 2.76 TeV", "LP");
			legMassPbPbPy->AddEntry(hSyPbPb,   "Systematic Pb-Pb" , "F");
			legMassPbPbPy->AddEntry(hMPy276,   "PYTHIA  Perugia 2011", "LP");
			cMassPbPbPy->cd(ih+1); legMassPbPbPy->Draw();
			
			legMassPbPb->AddEntry(hMPbPb,    "Pb-Pb 0-10% #sqrt{s_{NN}} = 2.76 TeV", "LP");
			legMassPbPb->AddEntry(hSyPbPb,   "Systematic Pb-Pb" , "F");
			legMassPbPb->AddEntry(hMPy276,   "PYTHIA  Perugia 2011 #sqrt{s_{NN}} = 2.76 TeV", "LP");
			legMassPbPb->AddEntry(hMJewPbPb, "JEWEL+PYTHIA Recoil on 0-10% PbPb", "LP");
			legMassPbPb->AddEntry(hMJewPbPbMRoff, "JEWEL+PYTHIA Recoil off 0-10% PbPb", "LP");
			legMassPbPb->AddEntry(hMQPy,      "Q-PYTHIA", "LP");
			cMassPbPb->cd(ih+1); legMassPbPb->Draw();
		}
		
		//ratio PbPb/pPb
		TH1* hRatioPbPbOpPb;
		TH1* hDividepPb;
		UniformTH1FForDivide(hMPbPb, hMpPb, hRatioPbPbOpPb, hDividepPb, "TH1D"); hRatioPbPbOpPb->SetName(Form("hRatioPbPbOpPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
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
		hRatioPbPbOpPbSys->SetName(Form("hRatioPbPbOpPbSys_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatioPbPbOpPbSys->GetYaxis()->SetTitle("#Rgothic_{#sqrt{#it{s}}}");
		hRatioPbPbOpPbSys->GetYaxis()->SetTitleOffset(1.15);
		//hRatioPbPbOpPbSys->GetYaxis()->SetTitleOffset(1.7);
		hRatioPbPbOpPbSys->Divide(hRDivpPb);
		hRatioPbPbOpPbSys->SetFillStyle(1001);
		hRatioPbPbOpPbSys->SetFillColor(kMagenta-10);
		hRatioPbPbOpPbSys->SetMarkerStyle(1);
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
		hRatioPbPbOpPbSysC->SetName(Form("hRatioPbPbOpPbSysC_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
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
		hRatioPbPbOpPbM->SetName(Form("hRatioPbPbOpPb_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
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
		hRatioPbPbOpPbSysM->SetName(Form("hRatioPbPbOpPbSysM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRDivpPbSyM->SetName(Form("hRDivpPbSyM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
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
		hRatio276O502M->GetYaxis()->SetTitleOffset(1.15);
		hRatio276O502M->Divide(hDividePyPy502M);
		hRatio276O502M->SetLineColor(kRed);
		hRatio276O502M->SetMarkerColor(kRed);
		hRatio276O502M->SetMarkerStyle(20);
		
		cRatioPbPbOpPb->cd(ih+1);
		
		
		
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
		hRatio276O502M->Draw("Lsames");
		// pythia ratio using 1134 output for pPb
		if(show1134) hRatio276O502->Draw("sames");
		//hPbPboppJewel->Draw("Lsames");
		//reference line
		lineOne->Draw();
		if(ih == 1) pvgeneral->Draw();
		pvpt[ih]->Draw();
		
		if(ih == nhM-1){
			legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPb, "Pb-Pb / p-Pb", "LP");// #sqrt{s_{NN}}
			legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPbSys, "Sys Pb-Pb/p-Pb", "F");
			if(pPbpaperprop){
				legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPbM, "Ratio PbPb/pPb paper prop", "LP");
				legRatioPbPbOpPb->AddEntry(hRatioPbPbOpPbSysM, "Sys PbPb/pPb paper prop", "F");
			}
			if(show1134) legRatioPbPbOpPb->AddEntry(hRatio276O502, "PYTHIA(2.76)/PYTHIA(5.02)Matched", "LP");
			
			legRatioPbPbOpPb->AddEntry(hRatio276O502M, "PYTHIA(2.76Tev)/PYTHIA(5.02TeV)", "LP"); // paper prop
			//legRatioPbPbOpPb->AddEntry(hPbPboppJewel, "Ratio JEWEL", "LP");
			legRatioPbPbOpPb->Draw();
			
		}
		
		//ratio pPb/PYTHIA
		TH1* hRatiopPbOPyM;
		TH1* hDividePy502M;
		UniformTH1FForDivide(hMpPb, hMPy502M, hRatiopPbOPyM, hDividePy502M, "TH1D"); hRatiopPbOPyM->SetName(Form("hRatiopPbOPyM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
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
		TH1* hDividePy502MNoErr = (TH1*)hDivideSysPy502M->Clone(Form("hDividePy502MNoErr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		for(Int_t ib = 0; ib<hDividePy502MNoErr->GetNbinsX(); ib++) hDividePy502MNoErr->SetBinError(ib, 0);
		
		hRatiopPbOPy502SysM->SetName(Form("hRatiopPbOPy502SysM_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOPy502SysM->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatiopPbOPy502SysM->GetYaxis()->SetTitleOffset(1.15);
		
		//hRatiopPbOPy502SysM->GetYaxis()->SetTitleOffset(1.7);
		hRatiopPbOPy502SysM->Divide(hDividePy502MNoErr);
		hRatiopPbOPy502SysM->SetFillStyle(hSypPb->GetFillStyle());
		hRatiopPbOPy502SysM->SetLineWidth(hSypPb->GetLineWidth());
		hRatiopPbOPy502SysM->SetLineColor(hSypPb->GetLineColor());
		
		//ratio pPb paper prop /PYTHIA
		TH1* hRatiopPbOPypappr;
		TH1* hDividePy502Mbis;
		UniformTH1FForDivide(hMpPbM, hMPy502M, hRatiopPbOPypappr, hDividePy502Mbis, "TH1D"); hRatiopPbOPypappr->SetName(Form("hRatiopPbMOPypappr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
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
		TH1* hDividePy502MbisNoErr = (TH1*)hDivideSysPy502Mbis->Clone(Form("hDividePy502MbisNoErr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		for(Int_t ib = 0; ib<hDividePy502MbisNoErr->GetNbinsX(); ib++) hDividePy502MbisNoErr->SetBinError(ib, 0);
		
		hRatiopPbOPy502Syspappr->SetName(Form("hRatiopPbOPy502Syspappr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatiopPbOPy502Syspappr->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatiopPbOPy502Syspappr->GetYaxis()->SetTitleOffset(1.15);
		
		
		//hRatiopPbOPy502Syspappr->GetYaxis()->SetTitleOffset(1.7);
		hRatiopPbOPy502Syspappr->Divide(hDividePy502MbisNoErr);
		hRatiopPbOPy502Syspappr->SetFillStyle(hSypPbM->GetFillStyle());
		hRatiopPbOPy502Syspappr->SetLineWidth(hSypPbM->GetLineWidth());
		hRatiopPbOPy502Syspappr->SetLineColor(hSypPbM->GetLineColor());
		hRatiopPbOPy502Syspappr->SetFillColor(hSypPbM->GetLineColor());
		hRatiopPbOPy502Syspappr->SetMarkerStyle(1);
		
		cRatiopPbOPy->cd(ih+1);
		hRatiopPbOPy502SysM->Draw("E2");
		hRatiopPbOPyM->Draw("sames");
		if(pPbpaperprop) {
			hRatiopPbOPy502Syspappr->Draw("E2sames");
			hRatiopPbOPypappr->Draw("sames");
			
		}
		lineOne->Draw();
		if(ih == 1) pvgeneral->Draw();
		pvpt[ih]->Draw();
		
		if(ih == nhM-1){
			legRatiopPbOPy->AddEntry(hRatiopPbOPyM, "p-Pb / PYTHIA(5.02TeV)", "PL");
			legRatiopPbOPy->AddEntry(hRatiopPbOPy502SysM, "(Sys p-Pb) / PYTHIA(5.02TeV)", "F");
			if(pPbpaperprop) {
				legRatiopPbOPy->AddEntry(hRatiopPbOPypappr, "pPb paper pr / PYTHIA(5.02TeV)", "PL");
				legRatiopPbOPy->AddEntry(hRatiopPbOPy502Syspappr, "(Sys pPb paper pr) / PYTHIA(5.02TeV)", "F");
			}
			legRatiopPbOPy->Draw();
		}
		if(show1134) {}
		//ratio PbPb/PYTHIA
		TH1* hRatioPbPbOPy;
		TH1* hDividePy276;
		UniformTH1FForDivide(hMPbPb, hMPy276, hRatioPbPbOPy, hDividePy276, "TH1D"); hRatioPbPbOPy->SetName(Form("hRatiopPbOPy_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		
		hRatioPbPbOPy->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatioPbPbOPy->GetYaxis()->SetTitleOffset(1.15);
		hRatioPbPbOPy->Divide(hDividePy276);
		hRatioPbPbOPy->GetXaxis()->SetRangeUser(0., maxM);
		hRatioPbPbOPy->GetYaxis()->SetRangeUser(0., 4);
		
		
		// sys ratio PbPb/PYTHIA
		TH1* hRatioPbPbOPy276Sys;
		UniformTH1FForDivide(hSyPbPb, hMPy276, hRatioPbPbOPy276Sys, hDividePy276, "TH1D");
		TH1* hDividePy276NoErr = (TH1*)hDividePy276->Clone(Form("hDividePy276NoErr_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		for(Int_t ib = 0; ib<hDividePy276NoErr->GetNbinsX(); ib++) hDividePy276NoErr->SetBinError(ib, 0);
		
		hRatioPbPbOPy276Sys->SetName(Form("hRatioPbPbOPy276Sys_Pt%.0f_%.0f", ptlims[ih],ptlims[ih+1]));
		hRatioPbPbOPy276Sys->GetYaxis()->SetTitle("Data/PYTHIA");
		hRatioPbPbOPy276Sys->GetYaxis()->SetTitleOffset(1.15);
		
		//hRatioPbPbOPy276Sys->GetYaxis()->SetTitleOffset(1.7);
		hRatioPbPbOPy276Sys->Divide(hDividePy276NoErr);
		hRatioPbPbOPy276Sys->SetFillStyle(1001);
		hRatioPbPbOPy276Sys->SetLineWidth(2);
		hRatioPbPbOPy276Sys->SetLineColor(hSyPbPb->GetLineColor());
		hRatioPbPbOPy276Sys->SetFillColor(hSyPbPb->GetFillColor());
				
		cRatioPbPbOPy->cd(ih+1);
		
		hRatioPbPbOPy276Sys->GetYaxis()->SetRangeUser(0., 4);
		hRatioPbPbOPy276Sys->Draw("E2");
		hRatioPbPbOPy->Draw("sames");
		lineOne->Draw();
		if(ih == nhM-1){
			legRatioPbPbOPy->AddEntry(hRatioPbPbOPy, "Pb-Pb / PYTHIA(2.76TeV)", "PL");
			legRatioPbPbOPy->AddEntry(hRatioPbPbOPy276Sys, "Sys Pb-Pb / PYTHIA(2.76TeV)", "F");
			legRatioPbPbOPy->Draw();
		}
		pvpt[ih]->Draw();
		if(ih == 1) pvgeneral->Draw();
		// data/ PYTHIA overlapped
		cRatioDataOPy->cd(ih+1);
		
		hRatiopPbOPy502SysM->GetYaxis()->SetRangeUser(0., 4.);
		hRatiopPbOPy502SysM->GetXaxis()->SetRangeUser(0., maxM);
		hRatioPbPbOPy276Sys->GetYaxis()->SetRangeUser(0., 4.);
		hRatioPbPbOPy276Sys->GetXaxis()->SetRangeUser(0., maxM);
		hRatioPbPbOPy276Sys->Draw("samesE2");
		hRatioPbPbOPy->Draw("sames");
		hRatiopPbOPy502SysM->Draw("samesE2");
		hRatiopPbOPyM->Draw("sames");
		pvpt[ih]->Draw();
		lineOne->Draw();
		if(ih == 1) pvgeneral->Draw();
		if(ih == nhM-1){
			legRatioDataOPy->AddEntry(hRatiopPbOPyM, "p-Pb / PYTHIA(5.02TeV)", "PL");
			legRatioDataOPy->AddEntry(hRatiopPbOPy502SysM, "Sys p-Pb / PYTHIA(5.02TeV)", "F");
			legRatioDataOPy->AddEntry(hRatioPbPbOPy, "Pb-Pb / PYTHIA(2.76TeV)", "PL");
			legRatioDataOPy->AddEntry(hRatioPbPbOPy276Sys, "Sys Pb-Pb / PYTHIA(2.76TeV)", "F");
			legRatioDataOPy->Draw();
		}
		
	}
	cMass->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cMasspPb->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cMassPbPb->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cMassPbPbOnly->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cMassPbPbPy->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cRatioPbPbOpPb->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cRatioDataOPy->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cRatioPbPbOPy->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	cRatiopPbOPy->cd(1);
	DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
	
	cRatiopPbOPy->cd(2);
	pvSystpPb->Draw();
	cRatiopPbOPy->cd(1);
	pvSystPyth->Draw();
	
	cRatioDataOPy->cd(2);
	pvSystpPb->Draw();
	cRatioDataOPy->cd(1);
	pvSystPyth->Draw();
	
	cRatioPbPbOpPb->cd(2);
	pvSystpPb->Draw();
	
	cRatioPbPbOPy->cd(2);
	pvSystPbPb->Draw();
	pvSystPyth->Draw();
	
	cRatioDataOPy->cd(2);
	pvSystPbPb->Draw();
	cRatioDataOPy->cd(1);
	pvSystPyth->Draw();
	cRatioPbPbOpPb->cd(2);
	pvSystPbPb->Draw();
	
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
	
	TH1D* hmean = new TH1D(namehmean, "Mean;#it{p}_{T}(GeV/#it{c}); <M_{ch jet}>(GeV/#it{c}^{2})", nbins, lims);
	hmean->Sumw2();
	
	hrms = new TH1D(Form("%sStdDev", namehmean), "StdDev;#it{p}_{T}(GeV/#it{c}); Std Dev(GeV/#it{c}^{2})", nbins, lims);
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
	
	hSysrms = new TH1D(Form("%sStdDev", namehsysmean), "StdDev;#it{p}_{T}(GeV/#it{c}); Std Dev(GeV/#it{c}^{2})", nbins, lims);
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

void DrawMeanComparison(Bool_t stylezero = kTRUE){
	
	if(stylezero) {
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		gStyle->SetTextFont(42);
	}
	
	const Int_t ninputs = 4;
	TString inputDistr[ninputs] = {
		"/data/Work/jets/JetMass/PbPbResults/TranformedIntoTH1.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160911/MasspPbResults.root",
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
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160915/TotalSystematicUnc.root"
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
	// -> then, somewhere else (e.g. CreateRooUnfoldResponse) do the prior variation ad described in the note: fill the response with a random number extracted from the projection on the reco level per each particle level bin
	
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
		cInput[ifi] = new TCanvas(Form("c%s", h1names[ifi].Data()), Form("%s", h1names[ifi].Data()), 800, 800);
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
				hnamethisloop =  Form("%s%d", hsysnames[ifi].Data(),  ib + offset[ifi]);
				hnamedoublethisloop = Form("%s%.0f_%.0fb", hsysnames[ifi].Data() ,  ptlims[ib + offset[ifi]],  ptlims[ib + 1 + offset[ifi]]);
				
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
			
			hSysmean[ifi] = CalculateSysMean(xrange, nptbins, ptlims, hSysdistr, grSysmean[ifi], hSysrms[ifi], Form("hSysMeanMass_%d_", ifi));
			
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
		legSy->AddEntry(hSysmean[ifi], Form("Systematics"), "F");
		
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
				hnamethisloop = Form("%s%d%s", hnamepref.Data(),  ib + offset[ifi], hnamesuff.Data());
				
				// for the comparison of the unfolded new
				
			}
			else {
				hnamethisloop = Form("%s%d", h1names[ifi].Data(),  ib + offset[ifi]);
				hnamedoublethisloop = Form("%s%.0f_%.0f", h1names[ifi].Data() ,  ptlims[ib + offset[ifi]],  ptlims[ib + 1 + offset[ifi]]);
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
		
		hmean[ifi] = CalculateMean(xrange, nptbins, ptlims, h1Mass, grmean[ifi], hrms[ifi], Form("hMeanMass_%d_", ifi));
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
#endif
