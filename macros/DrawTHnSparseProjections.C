#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>
#include <THnSparse.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TList.h>

void DrawTHnSparseProjections(Int_t axrange[], Int_t axproj[], const Int_t nbins, Double_t binlims[], const Int_t nfiles, TString finame[], TString hname[], TString legs[], TString analysisresultslistname[], TString fullname = "", Int_t base = -1, Int_t rebin = 1);

TH2D* Projections2DWithCuts(THnSparseF* hnsp, Int_t *axproj, Double_t *mincut, Double_t *maxcut, TString projname);

TH1D* Projections1DWithCuts(THnSparseF* hnsp, Int_t axproj, Double_t *mincut, Double_t *maxcut, TString projname);

THnSparseF* ReadInput(TString file, TString listname, TString hname);

void ClearSelections(THnSparseF *hsp);

//-------------------------------------------------------------------

// running the method on several files/axis projections/axis ranges

void RunDrawTHnSparseProjectionsdMOrdPt(){
	Int_t axrange1 = 3;
	Int_t axproj1 = 1;
	
	const Int_t nbins = 4;
	Double_t binlims[nbins+1] = {40., 60., 80., 100., 120.};
	const Int_t nfiles = 2;
	// these files contain dM dpT Mp ptp
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/analysis/DeltaMDeltaPt/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/DetFluc/DeltaMDeltaPt/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root"};
	
	TString hname1 = "fhResponseFinal";
	TString hname[nfiles] = {hname1, hname1};
	Int_t axrange[nfiles] = {axrange1, axrange1};
	Int_t axproj[nfiles]  = {axproj1, axproj1};
	TString legs[nfiles] = {"Deriv", "Const"};
	TString fullname = "dM";
	if(axproj1 == 1) fullname = "dPt";
	TString listname[nfiles] = {"", ""};
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

void RunDrawTHnSparseProjectionsRespSlices(){
	//Int_t axrange1 = 3; // pt par
	//Int_t axproj1 = 2;  // pt det
	Int_t axrange1 = 1; 
	Int_t axproj1 = 0;  
	
	const Int_t nbins = 4;
	Double_t binlims[nbins+1] = {40., 60., 80., 100., 120.};
	if(axrange1 == 1) {
		binlims[0] = 0.;
		binlims[1] = 5.;
		binlims[2] = 10.;
		binlims[3] = 30.;
	}
	const Int_t nfiles = 2;
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/analysis/WeightNormPerBin/DetFluc/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/DetFluc/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root"};
	//{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo/00resp_pt20_ptT10/responseDetFlBkgDeriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo/00resp_pt20_ptT10/responseDetFlBkgConst.root"};
	TString hname1 = "fhResponseFinal"; //"fhnReduced";
	TString hname[nfiles] = {hname1, hname1}; 
	Int_t axrange[nfiles] = {axrange1, axrange1};
	Int_t axproj[nfiles]  = {axproj1, axproj1};
	TString legs[nfiles] = {"Deriv", "Const"};
	TString fullname = "PtRec";
	TString listname[nfiles] = {"", ""};
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

void RunDrawTHnSparseProjectionsDerivOutCmp(){
	//Int_t axrange1 = 3; // pt par
	//Int_t axproj1 = 2;  // pt det
	Int_t axrange1 = 1; 
	Int_t axproj1 = 0;  
	
	const Int_t nbins = 4;
	Double_t binlims[nbins+1] = {40., 60., 80., 100., 120.};
	if(axrange1 == 1) {
		binlims[0] = 0.;
		binlims[1] = 5.;
		binlims[2] = 10.;
		binlims[3] = 30.;
	}
	const Int_t nfiles = 2;
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo/00resp_pt20_ptT10/responseDetFlBkgDeriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/DetFluc/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root"};
	TString hname[nfiles] = {"fhnReduced", "fhResponseFinal"};
	Int_t axrange[nfiles] = {axrange1, axrange1};
	Int_t axproj[nfiles]  = {axproj1, axproj1};
	TString legs[nfiles] = {"UnfInput", "TaskShape"};
	TString fullname = "MDerivSub";
	TString listname[nfiles] = {"", ""};
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

void RunDrawTHnSparseProjectionsConstOutCmp(){
	//Int_t axrange1 = 3; // pt par
	//Int_t axproj1 = 2;  // pt det
	Int_t axrange1 = 1; 
	Int_t axproj1 = 0;  
	
	const Int_t nbins = 4;
	Double_t binlims[nbins+1] = {40., 60., 80., 100., 120.};
	if(axrange1 == 1) {
		binlims[0] = 0.;
		binlims[1] = 5.;
		binlims[2] = 10.;
		binlims[3] = 30.;
	}
	const Int_t nfiles = 2;
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo/00resp_pt20_ptT10/responseDetFlBkgDeriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/DetFluc/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root"};
	TString hname[nfiles] = {"fhnReduced", "fhResponseFinal"}; 
	Int_t axrange[nfiles] = {axrange1, axrange1};
	Int_t axproj[nfiles]  = {axproj1, axproj1};
	TString legs[nfiles] = {"UnfInput", "TaskShape"};
	TString fullname = "MDerivSub";
	TString listname[nfiles] = {"", ""};
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

void RunDrawTHnSparseProjectionsRhoOrRhoM(){
	Int_t axrange1 = 3; //pt
	Int_t axproj1 = 1; //0 = rhom, 1 =rho
	const Int_t nbins = 2;
	Double_t binlims[nbins+1] = {0., 40., 120.}; // 60., 80., 100.,
	const Int_t nfiles = 2;
	// these files contain dM dpT Mp ptp
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/analysis/Rho/vsPtMSub/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/Rho/vsPtMSub/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root"};
	
	TString hname1 = "fhResponseFinal";
	TString hname[nfiles] = {hname1, hname1};
	Int_t axrange[nfiles] = {axrange1, axrange1};
	Int_t axproj[nfiles]  = {axproj1, axproj1};
	TString legs[nfiles] = {"Deriv", "Const"};
	TString fullname = "Rho";
	if(axproj1 == 1) fullname = "RhoM";
	TString listname[nfiles] = {"", ""};
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

void RunDrawTHnSparseProjectionsRhoOrRhoMDataResp(){
	//response
	Int_t axrange1 = 2; //pt
	Int_t axproj1 = 1; //0 = rhom, 1 =rho
	//data
	Int_t axrange2 = 2; //pt
	Int_t axproj2 = 0; // 0 = rho or 1 = rhom
	const Int_t nbins = 5;
	Double_t binlims[nbins+1] = {0., 40., 60., 80., 100.,120.}; // 
	
	const Int_t nfiles = 2;
	// the data file contains rho, pt_lead, centrality
	// 
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/RhoWeighted/analysis/BkgRho/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/MB/output1203-1204/cmpnew/AnalysisResults.root"};
	TString listname[nfiles] = {"", "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv"};
	
	Int_t axrange[nfiles] = {axrange1, axrange2};
	Int_t axproj[nfiles]  = {axproj1, axproj2};
	TString hname1 = "fhResponseFinal";
	TString hname2 = "fhnRhoVsRhoMVsLeadJetPtVsMassVsCent";//"fh2RhoVsLeadJetPtVsCent";// 
	TString hname[nfiles] = {hname1, hname2};
	
	TString legs[nfiles] = {"Resp", "Data"};
	TString fullname = "Rho";//fullname = "Rho";
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

void RunDrawTHnSparseProjectionsPtDataResp(){
	Int_t axrange1 = 0; //0 = rhom, 1 =rho
	Int_t axproj1 = 2;  //pt
	Int_t axrange2 = 0; // rhom or rho
	Int_t axproj2 = 1; //pt
	const Int_t nbins = 1;
	Double_t binlims[nbins+1] = {0., 20.}; // 
	
	const Int_t nfiles = 2;
	// the data file contains rho, pt_lead, centrality
	// 
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/analysis/Rho/vsPtMSub/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1191/AnalysisResults.root"};
	TString listname[nfiles] = {"", "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv"};
	
	Int_t axrange[nfiles] = {axrange1, axrange2};
	Int_t axproj[nfiles]  = {axproj1, axproj2};
	TString hname1 = "fhResponseFinal";
	TString hname2 = "fh2RhoMVsLeadJetPtVsCent";//"fh2RhoVsLeadJetPtVsCent";// 
	TString hname[nfiles] = {hname1, hname2};
	
	TString legs[nfiles] = {"Resp", "Data"};
	TString fullname = "Pt";//fullname = "Rho";
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

void DrawRhoVsCentralityData(){
	Int_t axrange2 = 2; //centrality
	//Int_t axproj2 = 1; //pt
	Int_t axproj2 = 0; // rhom or rho
	const Int_t nbins = 5;
	Double_t binlims[nbins+1] = {0., 20., 40., 60., 80., 100.}; // 
	
	const Int_t nfiles = 1;
	// the data file contains rho, pt_lead, centrality
	// 
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1191/AnalysisResults.root"};
	TString listname[nfiles] = {"JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv"};
	
	Int_t axrange[nfiles] = {axrange2};
	Int_t axproj[nfiles]  = {axproj2};
	
	TString hname2 = "fh2RhoVsLeadJetPtVsCent";// "fh2RhoMVsLeadJetPtVsCent";//
	TString hname[nfiles] = {hname2};
	
	TString legs[nfiles] = {"Data"};
	TString fullname = "Pt"; //"Rho";//fullname = "Rho";
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);

}

void DrawCentralityVsPtData(){
	Int_t axrange2 = 1; //pt
	Int_t axproj2 = 2; //centrality
	const Int_t nbins = 3;
	Double_t binlims[nbins+1] = {0., 40., 60., 80.}; // 
	
	const Int_t nfiles = 1;
	// the data file contains rho, pt_lead, centrality
	// 
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1191/AnalysisResults.root"};
	TString listname[nfiles] = {"JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv"};
	
	Int_t axrange[nfiles] = {axrange2};
	Int_t axproj[nfiles]  = {axproj2};
	
	TString hname2 = "fh2RhoVsLeadJetPtVsCent";// "fh2RhoMVsLeadJetPtVsCent";//
	TString hname[nfiles] = {hname2};
	
	TString legs[nfiles] = {"Data"};
	TString fullname = "PtCent"; 
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);

}

void RunDrawTHnSparseProjectionsPt(){
	//data
	Int_t axrange2 = 1; // M
	Int_t axproj2 = 0;  // pt
	//response
	Int_t axrange1 = 0;
	Int_t axproj1 = 2;	
	const Int_t nbins = 1;
	Double_t binlims[nbins+1] = {-20., 50.};
	const Int_t nfiles = 2;
	// these files contain dM dpT Mp ptp
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/InputDataEJEMB/dataRerunMacroMBEJETrigComb1/MassOutput.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160419/analysis/Deriv3Runs/DetFluc/FineBinning/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root"};
	
	TString hname1 = "fhResponseFinal";
	TString hname2 = "h2MPtTagged_BkgSub2_TrgCmb1";
	TString hname[nfiles] = {hname2, hname1};
	Int_t axrange[nfiles] = {axrange2, axrange1};
	Int_t axproj[nfiles]  = {axproj2, axproj1};
	TString legs[nfiles] = {"Data", "Resp"};
	TString fullname = "Pt";
	if(axproj1 == 1) fullname = "M";
	TString listname[nfiles] = {"", ""};
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

//-------------------------------------------------------------------

void RunDrawTHnSparseProjectionsM(){
	//data
	Int_t axrange2 = 0;  // pt
	Int_t axproj2  = 1; // M
	//response
	Int_t axrange1 = 2;
	Int_t axproj1 = 0;	
	const Int_t nbins = 4;
	Double_t binlims[nbins+1] = {40., 60., 80., 100., 120.};
	const Int_t nfiles = 2;
	// these files contain dM dpT Mp ptp
	TString finame[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/InputDataEJEMB/dataRerunMacroMBEJETrigComb1/MassOutput.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/analysis/NoSub/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root"};
	
	TString hname1 = "fhResponseFinal";
	TString hname2 = "h2MPtTagged_BkgSub1_TrgCmb1";
	TString hname[nfiles] = {hname2, hname1};
	Int_t axrange[nfiles] = {axrange2, axrange1};
	Int_t axproj[nfiles]  = {axproj2, axproj1};
	TString legs[nfiles] = {"Data", "Resp"};
	TString fullname = "Pt";
	if(axproj1 == 0) fullname = "M";
	TString listname[nfiles] = {"", ""};
	DrawTHnSparseProjections(axrange, axproj, nbins, binlims, nfiles, finame, hname, legs, listname);
}

//-------------------------------------------------------------------
// where things are done:

void DrawTHnSparseProjections(Int_t axrange[], Int_t axproj[], const Int_t nbins, Double_t binlims[], const Int_t nfiles, TString finame[], TString hname[], TString legs[], TString analysisresultslistname[], TString fullname, Int_t base, Int_t rebin){
	/// rebin: set when needed, but look into where is used, only in the TH3F
	
	THnSparseF *hsph[nfiles];
	TH3F  *h3[nfiles];
	TH2D  *h2[nfiles];
	TFile *fins[nfiles];
	Int_t nx, ny, dx, dy;
	CalculatePads(nbins, nx, ny, dx, dy);
	
	for(Int_t i = 0; i< nfiles; i++) fullname+=legs[i];
	TCanvas *cAnvas = new TCanvas(Form("c%s", fullname.Data()), "", 900, 700);
	cAnvas->Divide(nx, ny);
	TCanvas *cRatios = new TCanvas(Form("cRatios%s", fullname.Data()), "", 900, 700);
	cRatios->Divide(nx, ny);
	
	TLegend *leg = new TLegend(0.6, 0.4, 0.8, 0.6);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	
	if(base == -1) base = nfiles-1; //if default, use the last file as deniminator
	TH1D*** hRatioProj = new TH1D**[nfiles-1];
	TH1D* hRatioDenom[nbins];
	
	
	for(Int_t ifile = 0; ifile < nfiles; ifile++){
		hRatioProj[ifile] = new TH1D*[nbins];
		
		Printf("%d", ifile);
		for(Int_t ipt = 0 ; ipt < nbins; ipt++) hRatioProj[ifile][ipt] = 0x0;
		
		if(!analysisresultslistname[ifile].IsNull()) {
			h3[ifile] = 0x0;
			hsph[ifile] = 0x0;
			h2[ifile]   = 0x0;
			Printf("Reading %s", analysisresultslistname[ifile].Data());
			TList *list = ReadFile(finame[ifile], analysisresultslistname[ifile]);
			
			
			Printf("list = %p",list);
			
			if(!list){
				continue;
			}
			TObject* obj = list->FindObject(hname[ifile].Data());
			if(!obj) {
				Printf("No obj %s found", hname[ifile].Data());
				
			}
			TClass *objclas = obj->IsA();
			TString clname = objclas->GetName();
			Printf("TClass is %s", clname.Data());
			if(clname == "TH3D"){
				Printf("Reading TH3D");
				h3[ifile] = (TH3F*)list->FindObject(hname[ifile].Data());
				if(!h3[ifile]){
					list->ls();
					continue;
				} else {
					Printf("Info, TH3F %s found", hname[ifile].Data());
				}
			}
			if(clname == "THnSparseT<TArrayF>"){
				Printf("Reading THnSparse");
				
				hsph[ifile] = (THnSparseF*)list->FindObject(hname[ifile].Data());
				if(!hsph[ifile]) {
					Printf("%s not found", hname[ifile].Data());
					list->ls();
					continue;
				}
				
			}
		} else {
			h3[ifile] = 0x0;
			hsph[ifile] = 0x0;
			h2[ifile] = 0x0;
			//fins[ifile] = new TFile(finame[ifile]);
			TFile *fin = new TFile(finame[ifile]);
			if(!fin->IsOpen()){
				Printf("File %s not found", finame[ifile].Data());
				continue;
			}
			TClass *objclas = fin->Get(hname[ifile].Data())->IsA();
			TString clname = objclas->GetName();
			Printf("TClass is %s", clname.Data());
			if(clname == "TH2D"){
				Printf("Reading TH2D");
				h2[ifile] = (TH2D*)fin->Get(hname[ifile].Data());
				if(!h2[ifile]){
					Printf("TH2D %s not found, looking for THnSparse", hname[ifile].Data());
					fin->ls();
					continue;
				}
				
			}
			if(clname == "THnSparseT<TArrayF>"){
				Printf("Reading THnSparse");
				hsph[ifile] = (THnSparseF*)fin->Get(hname[ifile].Data());
				if(!hsph[ifile]) {
					Printf("%s not found", hname[ifile].Data());
					fin->ls();
					continue;
				}
				
			}
		}
		
		// first projection integrated
		// these histograms are internal to the file loop:
		// - M integrated in the other axes
		// - M in the full range defined by limits
		
		TCanvas *cFile = new TCanvas(Form("cDistributions%d", ifile), Form("Distributions %s", legs[ifile].Data()), 600, 600);
		//TCanvas *cRatioFile = new TCanvas(Form("cDistributionRatios%d", ifile), Form("Ratios %s", legs[ifile].Data()), 600, 600);
		
		TH1D *hM = 0x0;
		TH1D *hMInLimRange = 0x0;
		
		if(hsph[ifile]){
			TLegend *legpt = new TLegend(0.5, 0.4, 0.8, 0.8, legs[ifile]);
			legpt->SetBorderSize(0);
			legpt->SetFillStyle(0);
			
			Int_t ndim = hsph[ifile]->GetNdimensions();
			
			if(axrange[ifile] < 0 || axrange[ifile] > ndim-1){
				Printf("Check axis for binning %d (ndims %d) failed", axrange[ifile], ndim);
				continue;
			}
			
			if(axproj[ifile] < 0 || axproj[ifile] > ndim-1){
				Printf("Check axis for projection %d (ndims %d) failed", axproj[ifile], ndim);
				continue;
			}
			Printf("Binning on %s, projecting on %s", hsph[ifile]->GetAxis(axrange[ifile])->GetTitle(), hsph[ifile]->GetAxis(axproj[ifile])->GetTitle());
			
			hM = hsph[ifile]->Projection(axproj[ifile]);
			hM->SetName(Form("hM%d", ifile));
			hM->SetLineStyle(2);
			hM->SetLineWidth(2);
			hM->SetLineColor(kBlack);
			hM->Scale(1./hM->Integral("width"));
			
			cFile->cd();
			gPad->SetLogy();
			hM->Draw("");
			Int_t binrangeall[2] = {hsph[ifile]->GetAxis(axrange[ifile])->FindBin(binlims[0]), hsph[ifile]->GetAxis(axrange[ifile])->FindBin(binlims[nbins])};
			hsph[ifile]->GetAxis(axrange[ifile])->SetRange(binrangeall[0], binrangeall[1]);
			hMInLimRange = hsph[ifile]->Projection(axproj[ifile]);
			hMInLimRange->SetName(Form("hMInLimRange%d", ifile));
			hMInLimRange->SetLineWidth(2);
			hMInLimRange->SetLineColor(kGray+2);
			hMInLimRange->Scale(1./hMInLimRange->Integral("width"));
			cFile->cd();
			hMInLimRange->Draw("sames");
			legpt->AddEntry(hMInLimRange,  Form("%.0f < %s < %.0f", binlims[0], hsph[ifile]->GetAxis(axrange[ifile])->GetTitle(), binlims[nbins]), "LP");
			for(Int_t ipt = 0; ipt < nbins; ipt++){
	
				
				TPaveText *pvt = 0x0;
				if(ifile == 1) {
					pvt = new TPaveText(0.4, 0.6, 0.9, 0.8, "NDC");
					pvt->SetFillStyle(0);
					pvt->SetBorderSize(0);
					pvt->AddText(Form("%.0f < %s < %.0f", binlims[ipt], hsph[ifile]->GetAxis(axrange[ifile])->GetTitle(), binlims[ipt+1]));
				}
				Int_t binrange[2] = {hsph[ifile]->GetAxis(axrange[ifile])->FindBin(binlims[ipt]), hsph[ifile]->GetAxis(axrange[ifile])->FindBin(binlims[ipt+1] - 0.1)};
				
				Printf("PtBin = %d, Bin range %f-%f -> %d - %d", ipt, binlims[ipt], binlims[ipt+1] - 0.1, binrange[0], binrange[1]);
				
				hsph[ifile]->GetAxis(axrange[ifile])->SetRange(binrange[0], binrange[1]);
			
				TH1D* hproj = hsph[ifile]->Projection(axproj[ifile]);
				hproj->SetName(Form("hdM%s_%d", legs[ifile].Data(), ipt));
				hproj->SetLineColor(colors[ifile]);
				hproj->SetMarkerColor(colors[ifile]);
				hproj->SetMarkerStyle(20+ifile);
				
				if(ipt == 0) leg->AddEntry(hproj, legs[ifile], "LP");
				
				cAnvas->cd(ipt+1);
				gPad->SetLogy();
				hproj->Scale(1./hproj->Integral("width"));
				if(ifile == 0) hproj->Draw();
				else hproj->Draw("sames");
				if(ifile == 1) pvt->Draw();
				
				TH1D* hprojcopy = (TH1D*)hproj->Clone(Form("%s_color", hproj->GetName()));
				hprojcopy->SetMarkerColor(colors[ipt]);
				legpt->AddEntry(hprojcopy, Form("%.0f < %s < %.0f", binlims[ipt],  hsph[ifile]->GetAxis(axrange[ifile])->GetTitle(),  binlims[ipt+1]), "LP");
				cFile->cd();
				hprojcopy->Draw("sames");
				
				if(ifile == base) {
					hRatioDenom[ipt] = (TH1D*)hprojcopy->Clone(Form("hDenom%s_%d", legs[base].Data(), ipt));
					Printf("Filling denominator with file %d, thnsparse input", ifile);
				}
				else{
					if(!hRatioProj[ifile][ipt]) hRatioProj[ifile][ipt] = (TH1D*)hprojcopy->Clone(Form("h%sOver%s_%d", legs[ifile].Data(), legs[base].Data(), ipt));
					Printf("Filling numerator with file %d, thnsparse input", ifile);
				}
			}
			cFile->cd();
			legpt->Draw();
			SaveCv(cFile);
		}
		if(h3[ifile]){

			
			Printf("%p", h3[ifile]);
			TLegend *legpt = new TLegend(0.5, 0.4, 0.8, 0.8, legs[ifile]);
			legpt->SetBorderSize(0);
			legpt->SetFillStyle(0);
			Printf("Rebinning TH3 from");
			if(axproj[ifile] == 0) {
				printf("%f to ", h3[ifile]->GetXaxis()->GetBinWidth(1));
				if(rebin > 1) h3[ifile]->RebinX(rebin);
				Printf("%f", h3[ifile]->GetXaxis()->GetBinWidth(1));
				
				hM = h3[ifile]->ProjectionX();

			}
			if(axproj[ifile] == 1)  {
				printf("%f to ", h3[ifile]->GetYaxis()->GetBinWidth(1));
				if(rebin > 1) h3[ifile]->RebinY(rebin);
				Printf("%f", h3[ifile]->GetYaxis()->GetBinWidth(1));
				hM = h3[ifile]->ProjectionY();
			}
			if(axproj[ifile] == 2) {
				printf("%f to ", h3[ifile]->GetZaxis()->GetBinWidth(1));
				if(rebin > 1) h3[ifile]->RebinZ(rebin);
				Printf("%f", h3[ifile]->GetZaxis()->GetBinWidth(1));
				hM = h3[ifile]->ProjectionZ();
			}
			legpt->AddEntry(hM, "Full range", "LP");
			hM->SetName(Form("hM%d", ifile));
			hM->SetLineStyle(2);
			hM->SetLineWidth(2);
			hM->SetLineColor(kBlack);
			hM->Scale(1./hM->Integral("width"));
			
			cFile->cd();
			gPad->SetLogy();
			hM->Draw("hist");
			//Int_t binrangeall[2] = {hsph[ifile]->GetAxis(axrange[ifile])->FindBin(binlims[0]), hsph[ifile]->GetAxis(axrange[ifile])->FindBin(binlims[nbins])};
			//
			//if(axproj[ifile] == 0) hMInLimRange = h3[ifile]->ProjectionX("_p", binrangeall[0], binrangeall[1]);
			//if(axproj[ifile] == 0) hMInLimRange = h3[ifile]->ProjectionX("_p", binrangeall[0], binrangeall[1]);
			//if(axproj[ifile] == 0) hMInLimRange = h3[ifile]->ProjectionX("_p", binrangeall[0], binrangeall[1]);
			//
			//hMInLimRange->SetName(Form("hMInLimRange%d", ifile));
			//hMInLimRange->SetLineWidth(2);
			//hMInLimRange->SetLineColor(kGray+2);
			//hMInLimRange->Scale(1./hMInLimRange->Integral("width"));
			//cFile->cd();
			//hMInLimRange->Draw("sames");
			//legpt->AddEntry(hMInLimRange,  Form("%.0f < %s < %.0f GeV/#it{c}", binlims[0], hsph[ifile]->GetAxis(axrange[ifile])->GetTitle(), binlims[nbins]), "LP");
			TString titleleg = "";
			for(Int_t ipt = 0; ipt < nbins; ipt++){
				
				TPaveText *pvt = 0x0;
				if(ifile == 1) {
					pvt = new TPaveText(0.4, 0.6, 0.9, 0.8, "NDC");
					pvt->SetFillStyle(0);
					pvt->SetBorderSize(0);
				}
				
				TH1D* hproj = 0x0;
				Int_t binrange[2] = {1,1};
				if(axrange[ifile] == 0) {
					binrange[0] = h3[ifile]->GetXaxis()->FindBin(binlims[ipt]);
					binrange[1] = h3[ifile]->GetXaxis()->FindBin(binlims[ipt+1] - 0.01);
					titleleg = h3[ifile]->GetXaxis()->GetTitle();
				}
				
				if(axrange[ifile] == 1) {
					binrange[0] = h3[ifile]->GetYaxis()->FindBin(binlims[ipt]);
					binrange[1] = h3[ifile]->GetYaxis()->FindBin(binlims[ipt+1] - 0.01);
					titleleg = h3[ifile]->GetYaxis()->GetTitle();
				}
				
				if(axrange[ifile] == 2) {
					binrange[0] = h3[ifile]->GetZaxis()->FindBin(binlims[ipt]);
					binrange[1] = h3[ifile]->GetZaxis()->FindBin(binlims[ipt+1] - 0.01);
					titleleg = h3[ifile]->GetZaxis()->GetTitle();
				}
				
				TString namepj = Form("hdM%s_%d", legs[ifile].Data(), ipt);
				Printf("Ranges ax  %d (%s) %d-%d", axrange[ifile], h3[ifile]->GetYaxis()->GetTitle(), binrange[0], binrange[1]);
				
				if(axproj[ifile] == 0)  {
					if(axrange[ifile] == 1) hproj = h3[ifile]->ProjectionX(namepj, binrange[0], binrange[1], 0, -1);
					if(axrange[ifile] == 2) hproj = h3[ifile]->ProjectionX(namepj, 0, -1, binrange[0], binrange[1]);
				}
				Printf("0 ) %p", hproj);
				if(axproj[ifile] == 1)  {
					if(axrange[ifile] == 0) hproj = h3[ifile]->ProjectionY(namepj, binrange[0], binrange[1], 0, -1);
					if(axrange[ifile] == 2) hproj = h3[ifile]->ProjectionY(namepj, 0, -1, binrange[0], binrange[1]);
				}
				Printf("1 ) %p", hproj);
				if(axproj[ifile] == 2)  {
					if(axrange[ifile] == 0) hproj = h3[ifile]->ProjectionZ(namepj, binrange[0], binrange[1], 0, -1);
					if(axrange[ifile] == 1) hproj = h3[ifile]->ProjectionZ(namepj, 0, -1, binrange[0], binrange[1]);
				}
				if(!hproj) {
					Printf("ERROR! %p", hproj);
					continue;	
				}
				hproj->Scale(1./hproj->Integral("width"));
				hproj->SetLineColor(colors[ifile]);
				hproj->SetMarkerColor(colors[ifile]);
				hproj->SetMarkerStyle(20+ifile);
				
				if(ipt == 0) leg->AddEntry(hproj, legs[ifile], "LP");
				if(ifile == 1) pvt->AddText(Form("%.0f < %s < %.0f", binlims[ipt], titleleg.Data(), binlims[ipt+1]));
				
				cAnvas->cd(ipt+1);
				gPad->SetLogy();
				if(ifile == 0) hproj->Draw("P");
				else hproj->Draw("Psames");
				if(ifile == 1) pvt->DrawClone();
				
				TH1D* hprojcopy = (TH1D*)hproj->Clone(Form("%s_color", hproj->GetName()));
				hprojcopy->SetMarkerColor(colors[ipt]);
				hprojcopy->SetFillColor(colors[ipt]);
				legpt->AddEntry(hprojcopy, Form("%.0f < %s < %.0f", binlims[ipt],  titleleg.Data(), binlims[ipt+1]) , "LP");
				cFile->cd();
				hprojcopy->Draw("sames");
				
				if(ifile == base) {
					hRatioDenom[ipt] = (TH1D*)hprojcopy->Clone(Form("hDenom%s_%d", legs[base].Data(), ipt));
					Printf("Filling denominator with file %d, th3f input", ifile);
				}
				else{
					if(!hRatioProj[ifile][ipt]) hRatioProj[ifile][ipt] = (TH1D*)hprojcopy->Clone(Form("h%sOver%s_%d", legs[ifile].Data(), legs[base].Data(), ipt));
					Printf("Filling numerator with file %d, th3f input", ifile);
				}
			}
			cFile->cd();
			legpt->Draw();
			SaveCv(cFile);
		}
		
		if(h2[ifile]){
			
			Printf("%p", h2[ifile]);
			TLegend *legpt = new TLegend(0.5, 0.4, 0.8, 0.8, legs[ifile]);
			legpt->SetBorderSize(0);
			legpt->SetFillStyle(0);
			Printf("Rebinning TH2 from");
			if(axproj[ifile] == 0) {
				printf("%f to ", h2[ifile]->GetXaxis()->GetBinWidth(1));
				if(rebin > 1) h2[ifile]->RebinX(rebin);
				Printf("%f", h2[ifile]->GetXaxis()->GetBinWidth(1));
				
				hM = h2[ifile]->ProjectionX();

			}
			
			if(axproj[ifile] == 1)  {
				printf("%f to ", h2[ifile]->GetYaxis()->GetBinWidth(1));
				if(rebin > 1) h2[ifile]->RebinY(rebin);
				Printf("%f", h2[ifile]->GetYaxis()->GetBinWidth(1));
				hM = h2[ifile]->ProjectionY();
			}
			
			legpt->AddEntry(hM, "Full range", "LP");
			hM->SetName(Form("hM%d", ifile));
			hM->SetLineStyle(2);
			hM->SetLineWidth(2);
			hM->SetLineColor(kBlack);
			hM->Scale(1./hM->Integral("width"));
			
			cFile->cd();
			gPad->SetLogy();
			hM->Draw("hist");
			
			TString titleleg = "";
			for(Int_t ipt = 0; ipt < nbins; ipt++){
				
				TPaveText *pvt = 0x0;
				if(ifile == 1) {
					pvt = new TPaveText(0.4, 0.6, 0.9, 0.8, "NDC");
					pvt->SetFillStyle(0);
					pvt->SetBorderSize(0);
				}
				
				TH1D* hproj = 0x0;
				Int_t binrange[2] = {1,1};
				if(axrange[ifile] == 0) {
					binrange[0] = h2[ifile]->GetXaxis()->FindBin(binlims[ipt]);
					binrange[1] = h2[ifile]->GetXaxis()->FindBin(binlims[ipt+1] - 0.01);
					titleleg = h2[ifile]->GetXaxis()->GetTitle();
				}
				
				if(axrange[ifile] == 1) {
					binrange[0] = h2[ifile]->GetYaxis()->FindBin(binlims[ipt]);
					binrange[1] = h2[ifile]->GetYaxis()->FindBin(binlims[ipt+1] - 0.01);
					titleleg = h2[ifile]->GetYaxis()->GetTitle();
				}
				
				TString namepj = Form("hdM%s_%d", legs[ifile].Data(), ipt);
				Printf("Ranges ax  %d (%s) %d-%d", axrange[ifile], h2[ifile]->GetYaxis()->GetTitle(), binrange[0], binrange[1]);
				
				if(axproj[ifile] == 0)  {
					if(axrange[ifile] == 1) hproj = h2[ifile]->ProjectionX(namepj, binrange[0], binrange[1]);
				}
				Printf("0 ) %p", hproj);
				if(axproj[ifile] == 1)  {
					if(axrange[ifile] == 0) hproj = h2[ifile]->ProjectionY(namepj, binrange[0], binrange[1]);
				}
				
				if(!hproj) {
					Printf("ERROR! %p", hproj);
					continue;	
				}
				hproj->Scale(1./hproj->Integral("width"));
				hproj->SetLineColor(colors[ifile]);
				hproj->SetMarkerColor(colors[ifile]);
				hproj->SetMarkerStyle(20+ifile);
				
				if(ipt == 0) leg->AddEntry(hproj, legs[ifile], "LP");
				if(ifile == 1) pvt->AddText(Form("%.0f < %s < %.0f", binlims[ipt], titleleg.Data(), binlims[ipt+1]));
				
				cAnvas->cd(ipt+1);
				gPad->SetLogy();
				if(ifile == 0) hproj->Draw("P");
				else hproj->Draw("Psames");
				if(ifile == 1) pvt->DrawClone();
				
				TH1D* hprojcopy = (TH1D*)hproj->Clone(Form("%s_color", hproj->GetName()));
				hprojcopy->SetMarkerColor(colors[ipt]);
				hprojcopy->SetFillColor(colors[ipt]);
				legpt->AddEntry(hprojcopy, Form("%.0f < %s < %.0f", binlims[ipt],  titleleg.Data(), binlims[ipt+1]) , "LP");
				cFile->cd();
				hprojcopy->Draw("sames");
				
				if(ifile == base) {
					hRatioDenom[ipt] = (TH1D*)hprojcopy->Clone(Form("hDenom%s_%d", legs[base].Data(), ipt));
					Printf("Filling denominator with file %d, th2f input", ifile);
				}
				else{
					if(!hRatioProj[ifile][ipt]) hRatioProj[ifile][ipt] = (TH1D*)hprojcopy->Clone(Form("h%sOver%s_%d", legs[ifile].Data(), legs[base].Data(), ipt));
					Printf("Filling numerator with file %d, th2f input", ifile);
				}
			}
			cFile->cd();
			legpt->Draw();
			SaveCv(cFile);
		}
		
	}
	cAnvas->cd(1);
	leg->Draw();
	
	for(Int_t ifile = 0; ifile < nfiles; ifile++){
		for(Int_t ipt = 0 ; ipt < nbins; ipt++) {
			cRatios->cd(ipt+1);
			Printf("Ifile %d , pt %d", ifile, ipt);
			if(ifile != base && hRatioProj[ifile][ipt] && hRatioDenom[ipt]){
				TH1 *hNum = 0x0;
				TH1 *hDen = 0x0;
				Printf("%p/%p", hRatioProj[ifile][ipt], hRatioDenom[ipt]);
				Int_t exit = UniformTH1FForDivide(hRatioProj[ifile][ipt], hRatioDenom[ipt], hNum, hDen, "TH1D", kFALSE);
				if(exit < 0) {
					Printf("Error");
					continue;	
				}
				hNum->Divide(hDen); 
				hNum->Draw("sames");
				//hRatioProj[ifile][ipt]->Divide(hRatioDenom[ipt]);
				//hRatioProj[ifile][ipt]->Draw("sames");
				//hRatioDenom[ipt]->Draw("sames");
			}
		}
	}
	SaveCv(cAnvas);
}


void ReadAndProject(Int_t bkgSubMethod = 0, TString fileEmbFlucOnly = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1088/output/AnalysisResults.root", Int_t rhoOrrhoM = 1){

	// study the correlation between the dM and dPt from fluctuation and rho
	// bkgSubMethod: 0 = Deriv, 1 = Const Sub, -1 = no sub
	//rhoOrrhoM: rho = 1, rho_m =2
	
	TString listDeriv = "JetShapeDeriv_JetEmb_AKTChargedR040_PicoTracksEmb_pT0150_E_scheme_TC";
	TString listConst = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracksEmb_pT0150_E_scheme_TC";
	
	TString hnameFlucVsRho = "fhnDeltaMassAndBkgInfo";
	
	TString bkgsub = "Deriv";
	THnSparseF *hSparse = 0x0;
	if(bkgSubMethod == 0) {
		hSparse = ReadInput(fileEmbFlucOnly, listDeriv, hnameFlucVsRho);
	}
	if(TMath::Abs(bkgSubMethod) == 1) { // includes bkgSubMethod == -1. The dM and dpT read is the one without background subtraction (see below)
		hSparse = ReadInput(fileEmbFlucOnly, listConst, hnameFlucVsRho);
		bkgsub = "Const";
	}
	if(bkgSubMethod == -1) {
		bkgsub = "NoBkg";
	}
	
	Int_t axprojM[2]  = {0, 8}; //Delta mass, rho
	Int_t axprojPt[2] = {1, 8}; //Delta Pt, rho
	if(bkgSubMethod<0) {
		Printf("Getting the no sub axes");
		axprojM[0] = 2;
		axprojPt[0] = 3;
	}
	if(rhoOrrhoM == 2) {
		axprojM[1] = 9;
		axprojPt[1] = 9;
	}
	const Int_t ndim = 10;
	// dM, dpt, Munsub- Mpart, pTunsub - pTpar #it{M}_{det}, #it{M}_{unsub}, #it{p}_{T,det}, #it{p}_{T,unsub}, #rho, #rho_m
	Int_t nbinsM = 4;
	Double_t binsM[nbinsM] = {0, 5, 10, 30};
	//Int_t nbinsPt = 9;
	//Double_t binsPt[nbinsPt] = {0, 10, 20, 30, 40, 60, 80, 100., 120};
	Int_t nbinsPt = 4;
	Double_t binsPt[nbinsPt] = {20, 40, 60, 80};
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nbinsPt-1, nx, ny, dx, dy);
	
	TCanvas *cdMvsRhoPtBins = new TCanvas(Form("cdMvsRho%sPtBins%s", rhoOrrhoM == 2 ? "M" : "", bkgsub.Data()), "dM vs Rho (PtBins)", 900, 800);
	cdMvsRhoPtBins->Divide(nx, ny);
	
	TCanvas *cdPtvsRhoPtBins = new TCanvas(Form("cdPtvsRho%sPtBins%s", rhoOrrhoM == 2 ? "M" : "", bkgsub.Data()), "dPt vs Rho (PtBins)", 900, 800);
	cdPtvsRhoPtBins->Divide(nx, ny);
	
	
	CalculatePads(nbinsM-1, nx, ny, dx, dy);
	TCanvas *cdMvsRhoMBins = new TCanvas(Form("cdMvsRho%sMBins%s", rhoOrrhoM == 2 ? "M" : "", bkgsub.Data()), "dM vs Rho (MBins)", 900, 800);
	cdMvsRhoMBins->Divide(nx, ny);
	
	TCanvas *cdPtvsRhoMBins = new TCanvas(Form("cdPtvsRho%sMBins%s", rhoOrrhoM == 2 ? "M" : "", bkgsub.Data()), "dPt vs Rho (MBins)", 900, 800);
	cdPtvsRhoMBins->Divide(nx, ny);
	
	for(Int_t ipt = 0; ipt<nbinsPt-1; ipt++){
		
		TPaveText *pave = new TPaveText(0.35, 0.8, 0.75, 0.9, "NDC");
		pave->SetFillStyle(0);
		pave->SetBorderSize(0);
		pave->AddText(Form("%.0f < #it{p}_{T} (GeV/#it{c}) < %.0f", binsPt[ipt], binsPt[ipt+1]));
		//start with low mass region and bins of pT
		//Double_t mincut[ndim] = {-99, -99, binsM[0], -99, binsPt[ipt], -99, -99, -99};
		//Double_t maxcut[ndim] = {-99, -99, binsM[nbinsM-1], -99, binsPt[ipt+1]-0.01, -99, -99, -99};
		Double_t mincut[ndim] = {-99, -99, -99, -99, -99, binsM[0], -99, binsPt[ipt], -99, -99};
		Double_t maxcut[ndim] = {-99, -99, -99, -99, -99, binsM[nbinsM-1], -99, binsPt[ipt+1]-0.01, -99, -99};
		
		TH2D* hPjdMDeriv = Projections2DWithCuts(hSparse, axprojM, mincut, maxcut, Form("hdMvsrho_PtUnsub%d", ipt));
		hPjdMDeriv->GetZaxis()->SetRangeUser(10, 1e5);
		
		TH2D* hPjdPtDeriv = Projections2DWithCuts(hSparse, axprojPt, mincut, maxcut, Form("hdPtvsrho_PtUnsub%d", ipt));
		hPjdPtDeriv->GetZaxis()->SetRangeUser(10, 1e5);
		
		cdMvsRhoPtBins->cd(ipt+1);
		gPad->SetLogz();
		gPad->SetGridy();
		hPjdMDeriv->Draw("colz");
		pave->Draw();
		
		cdPtvsRhoPtBins->cd(ipt+1);
		gPad->SetLogz();
		gPad->SetGridy();
		hPjdPtDeriv->Draw("colz");
		pave->Draw();
		
		
	} //the output of this loop shows no dependence on pT
	
	for(Int_t ipt = 0; ipt<nbinsM-1; ipt++){
		
		TPaveText *pave = new TPaveText(0.35, 0.8, 0.75, 0.9, "NDC");
		pave->SetFillStyle(0);
		pave->SetBorderSize(0);
		pave->AddText(Form("%.0f < #it{M} (GeV) < %.0f", binsM[ipt], binsM[ipt+1]));
		//start with low mass region and bins of pT
		//Double_t mincut[ndim] = {-99, -99, binsM[ipt], -99, binsPt[0], -99, -99, -99};
		//Double_t maxcut[ndim] = {-99, -99, binsM[ipt+1]-0.01, -99, binsPt[nbinsPt-1], -99, -99, -99};
		
		Double_t mincut[ndim] = {-99, -99, -99, -99, -99, binsM[ipt], -99, binsPt[0], -99, -99};
		Double_t maxcut[ndim] = {-99, -99, -99, -99, -99, binsM[ipt+1]-0.01, -99, binsPt[nbinsPt-1], -99, -99};
		TH2D* hPjdMDeriv = Projections2DWithCuts(hSparse, axprojM, mincut, maxcut, Form("hdMvsrho_MUnsub%d", ipt));
		hPjdMDeriv->GetZaxis()->SetRangeUser(10, 1e5);
		
		TH2D* hPjdPtDeriv = Projections2DWithCuts(hSparse, axprojPt, mincut, maxcut, Form("hdPtvsrho_MUnsub%d", ipt));
		hPjdPtDeriv->GetZaxis()->SetRangeUser(10, 1e5);
		
		cdMvsRhoMBins->cd(ipt+1);
		gPad->SetLogz();
		gPad->SetGridy();
		hPjdMDeriv->Draw("colz");
		pave->Draw();
		
		cdPtvsRhoMBins->cd(ipt+1);
		gPad->SetLogz();
		gPad->SetGridy();
		hPjdPtDeriv->Draw("colz");
		pave->Draw();
	} //the output of this loop shows no dependence on M
	
	SaveCv(cdMvsRhoPtBins);
	SaveCv(cdPtvsRhoPtBins);
	SaveCv(cdMvsRhoMBins);
	SaveCv(cdPtvsRhoMBins);
}

void ComparisonDerivConst1D(TString fileEmbFlucOnly = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1088/output/AnalysisResults.root", Int_t rhoOrrhoM = 1, Bool_t stylezero = kTRUE, Int_t binFig = 2){

	// study the correlation between the dM and dPt from fluctuation and rho
	// bkgSubMethod: Deriv, Const Sub, no Sub
	//rhoOrrhoM: rho = 1, rho_m =2
	
	if(stylezero) {
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		gStyle->SetTextFont(42);
		
		gStyle->SetLabelSize(0.05, "X");
		gStyle->SetTitleSize(0.048, "X");
		
		gStyle->SetLabelSize(0.05, "Y");
		gStyle->SetTitleSize(0.048, "Y");
		
		gStyle->SetPadBottomMargin(.13);
		gStyle->SetPadLeftMargin(.13);
		gStyle->SetPadRightMargin(.11);
	}
	
	TString listDeriv = "JetShapeDeriv_JetEmb_AKTChargedR040_PicoTracksEmb_pT0150_E_scheme_TC";
	TString listConst = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracksEmb_pT0150_E_scheme_TC";
	
	TString hnameFlucVsRho = "fhnDeltaMassAndBkgInfo";
	
	TString bkgsub[2] = {"Deriv", "Const"};
	THnSparseF *hSparseDeriv = ReadInput(fileEmbFlucOnly, listDeriv, hnameFlucVsRho);
	
	THnSparseF *hSparseConst = ReadInput(fileEmbFlucOnly, listConst, hnameFlucVsRho);
	
	Int_t axprojM[2]  = {0, 8}; //Delta mass, rho
	Int_t axprojPt[2] = {1, 8}; //Delta Pt, rho
	if(rhoOrrhoM == 2) {
		axprojM[1] = 9;
		axprojPt[1] = 9;
	}
	Int_t axNoBkgdM = 2, axNoBkgdpT = 3;
	Int_t axprojFluc[2] = {0, 1};
	Int_t axprojFlucNoBkg[2] = {axNoBkgdM, axNoBkgdpT};//{0, 1};
	
	const Int_t ndim = 10;
	// dM, dpt, #it{M}_{det}, #it{M}_{unsub}, #it{p}_{T,det}, #it{p}_{T,unsub}, #rho, #rho_m
	Int_t nbinsM = 4;
	Double_t binsM[nbinsM] = {0, 5, 10, 30};
	Int_t nbinsPt = 9;
	Double_t binsPt[nbinsPt] = {0, 10, 20, 30, 40, 60, 80, 100., 120};
	Int_t nbinsrho = 6;
	Double_t binsrho[nbinsrho] = {0, 1.5, 3, 4.5, 10, 20};
	TString projtxt = "#rho";
	if(rhoOrrhoM == 2){ //bins for rhom
		binsrho[0] = 0;
		binsrho[1] = 0.2;
		binsrho[2] = 0.4;
		binsrho[3] = 0.6;
		binsrho[4] = 0.8;
		binsrho[5] = 1;
		projtxt = "#rho_{M}";
	}
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nbinsrho-1, nx, ny, dx, dy, 2);
	
	TCanvas *cdMRhoBins = new TCanvas(Form("cdMRho%sBins", rhoOrrhoM==2 ? "M" : ""), Form("Delta M in rho %s bins", rhoOrrhoM==2 ? "M" : ""), 900, 800);
	cdMRhoBins->Divide(nx, ny);
	TCanvas *cdPtRhoBins = new TCanvas(Form("cdPtRho%sBins", rhoOrrhoM==2 ? "M" : ""), Form("Delta Pt in rho %s bins", rhoOrrhoM==2 ? "M" : ""), 900, 800);
	cdPtRhoBins->Divide(nx, ny);
	
	TCanvas *cdMRatiosRhoBins = new TCanvas(Form("cdMRatiosRho%sBins", rhoOrrhoM==2 ? "M" : ""), Form("Delta M Ratio in rho %s bins", rhoOrrhoM==2 ? "M" : ""), 900, 800);
	cdMRatiosRhoBins->Divide(nx, ny);
	TCanvas *cdPtRatiosRhoBins = new TCanvas(Form("cdPtRatiosRho%sBins", rhoOrrhoM==2 ? "M" : ""), Form("Delta M Ratio in rho %s bins", rhoOrrhoM==2 ? "M" : ""), 900, 800);
	cdPtRatiosRhoBins->Divide(nx, ny);
	
	TCanvas *cdMdPtRhoBinsDeriv = new TCanvas(Form("cdMdPtRho%sBinsDeriv", rhoOrrhoM==2 ? "M" : ""), Form("Delta M vs delta Ptin rho %s bins (Deriv)", rhoOrrhoM==2 ? "M" : ""), 900, 800);
	cdMdPtRhoBinsDeriv->Divide(nx, ny);
	
	TCanvas *cdMdPtRhoBinsConst = new TCanvas(Form("cdMdPtRho%sBinsConst", rhoOrrhoM==2 ? "M" : ""), Form("Delta M vs delta Ptin rho %s bins (Const)", rhoOrrhoM==2 ? "M" : ""), 900, 800);
	cdMdPtRhoBinsConst->Divide(nx, ny);
	TCanvas *cdMdPtRhoBinsNoBkg = new TCanvas(Form("cdMdPtRho%sBinsNoBkg", rhoOrrhoM==2 ? "M" : ""), Form("Delta M vs delta Ptin rho %s bins (NoBkg)", rhoOrrhoM==2 ? "M" : ""), 900, 800);
	cdMdPtRhoBinsNoBkg->Divide(nx, ny);
	
	TCanvas *cdMdPtRhoBinsDerivFig = new TCanvas(Form("cdMdPtRho%sBinsDerivFig", rhoOrrhoM==2 ? "M" : ""), Form("Delta M vs delta Ptin rho %s bins (Deriv)", rhoOrrhoM==2 ? "M" : ""), 400, 400);
	
	TCanvas *cdMdPtRhoBinsConstFig = new TCanvas(Form("cdMdPtRho%sBinsConstFig", rhoOrrhoM==2 ? "M" : ""), Form("Delta M vs delta Ptin rho %s bins (Const)", rhoOrrhoM==2 ? "M" : ""), 400, 400);
	      
	TCanvas *cdMdPtRhoBinsNoBkgFig = new TCanvas(Form("cdMdPtRho%sBinsNoBkgFig", rhoOrrhoM==2 ? "M" : ""), Form("Delta M vs delta Ptin rho %s bins (NoBkg)", rhoOrrhoM==2 ? "M" : ""), 400, 400);
	
	TLegend *legBkg = new TLegend(0.6, 0.4, 0.8, 0.6);
	legBkg->SetFillStyle(0);
	legBkg->SetBorderSize(0);
	
	for(Int_t ipt = 0; ipt<nbinsrho-1; ipt++){
		
		TPaveText *pave = new TPaveText(0.35, 0.8, 0.75, 0.9, "NDC");
		pave->SetFillStyle(0);
		pave->SetBorderSize(0);
		pave->AddText(Form("%.1f < %s (GeV) < %.1f", binsrho[ipt], projtxt.Data(), binsrho[ipt+1]));
		//start with low mass region and bins of pT
		//Double_t mincut[ndim] = {-99, -99, binsM[0], -99, binsPt[0], -99, binsrho[ipt], -99};
		//Double_t maxcut[ndim] = {-99, -99, binsM[nbinsM-1], -99, binsPt[nbinsPt-1], -99, binsrho[ipt+1]-0.01, -99};
	
		Double_t mincut[ndim] = {-99, -99, -99, -99, binsM[0], -99, binsPt[0], -99, binsrho[ipt], -99};
		Double_t maxcut[ndim] = {-99, -99, -99, -99, binsM[nbinsM-1], -99, binsPt[nbinsPt-1], -99, binsrho[ipt+1]-0.01, -99};

		//Double_t mincut[ndim] = {-99, -99, -99, -99, -99, binsM[0], -99, binsPt[0], binsrho[ipt], -99};
		//Double_t maxcut[ndim] = {-99, -99, -99, -99, -99, binsM[nbinsM-1], -99, binsPt[nbinsPt-1], binsrho[ipt+1]-0.01, -99};
		ClearSelections(hSparseDeriv);
		TH1D* hPjdMDeriv = Projections1DWithCuts(hSparseDeriv, axprojM[0], mincut, maxcut, Form("hdM_Deriv_Rho%d", ipt));
		ClearSelections(hSparseConst);
		TH1D* hPjdMConst = Projections1DWithCuts(hSparseConst, axprojM[0], mincut, maxcut, Form("hdM_Const_Rho%d", ipt));
		ClearSelections(hSparseDeriv);
		TH1D* hPjdMNoBkg = Projections1DWithCuts(hSparseDeriv, axNoBkgdM, mincut, maxcut, Form("hdM_NoBkg_Rho%d", ipt));
		
		
		hPjdMDeriv->SetLineColor(colors[0]);
		hPjdMConst->SetLineColor(colors[1]);
		hPjdMNoBkg->SetLineColor(colors[2]);
		if(ipt == 0){
			legBkg->AddEntry(hPjdMDeriv, bkgsub[0], "l");
			legBkg->AddEntry(hPjdMConst, bkgsub[1], "l");
			legBkg->AddEntry(hPjdMNoBkg, "NoBkg"  , "l");
			
		}
		ClearSelections(hSparseDeriv);
		TH1D* hPjdPtDeriv = Projections1DWithCuts(hSparseDeriv, axprojPt[0], mincut, maxcut, Form("hdPt_Deriv_Rho%d", ipt));
		ClearSelections(hSparseConst);
		TH1D* hPjdPtConst = Projections1DWithCuts(hSparseConst, axprojPt[0], mincut, maxcut, Form("hdPt_Const_Rho%d", ipt));
		ClearSelections(hSparseConst);
		TH1D* hPjdPtNoBkg = Projections1DWithCuts(hSparseConst, axNoBkgdpT, mincut, maxcut, Form("hdPt_NoBkg_Rho%d", ipt));
		
		hPjdPtDeriv->SetLineColor(colors[0]);
		hPjdPtConst->SetLineColor(colors[1]);
		hPjdPtNoBkg->SetLineColor(colors[2]);
		
		cdMRhoBins->cd(ipt+1);
		gPad->SetLogy();
		//gPad->SetGridy();
		hPjdMDeriv->Draw();
		hPjdMConst->Draw("sames");
		hPjdMNoBkg->Draw("sames");
		pave->Draw();
		legBkg->Draw();
		
		
		cdPtRhoBins->cd(ipt+1);
		gPad->SetLogy();
		//gPad->SetGridy();
		hPjdPtDeriv->Draw();
		hPjdPtConst->Draw("sames");
		hPjdPtNoBkg->Draw("sames");
		pave->Draw();
		legBkg->Draw();
		
		//Ratios
		TH1D* hPjRatioPtDeriv = (TH1D*)hPjdPtDeriv->Clone("hPjRatioPtDeriv");
		TH1D* hPjRatioPtConst = (TH1D*)hPjdPtConst->Clone("hPjRatioPtConst");
		hPjRatioPtDeriv->Add(hPjdPtNoBkg, -1);
		hPjRatioPtConst->Add(hPjdPtNoBkg, -1);
		
		TH1D* hPjRatioMDeriv = (TH1D*)hPjdMDeriv->Clone("hPjRatioMDeriv");
		TH1D* hPjRatioMConst = (TH1D*)hPjdMConst->Clone("hPjRatioMConst");
		hPjRatioMDeriv->Add(hPjdMNoBkg, -1);
		hPjRatioMConst->Add(hPjdMNoBkg, -1);
		
		hPjRatioPtDeriv->GetXaxis()->SetRangeUser(-10, 10);
		hPjRatioPtConst->GetXaxis()->SetRangeUser(-10, 10);
		hPjRatioMDeriv-> GetXaxis()->SetRangeUser(-10, 10);
		hPjRatioMConst-> GetXaxis()->SetRangeUser(-10, 10);
		
		hPjRatioPtDeriv->GetYaxis()->SetTitle("BkgSub-NoSub");
		hPjRatioPtConst->GetYaxis()->SetTitle("BkgSub-NoSub");
		hPjRatioMDeriv-> GetYaxis()->SetTitle("BkgSub-NoSub");
		hPjRatioMConst-> GetYaxis()->SetTitle("BkgSub-NoSub");
		
		cdPtRatiosRhoBins->cd(ipt+1);
		hPjRatioPtDeriv->Draw();
		hPjRatioPtConst->Draw("sames");
		pave->Draw();
		legBkg->Draw();
		
		
		cdMRatiosRhoBins->cd(ipt+1);
		hPjRatioMDeriv->Draw();
		hPjRatioMConst->Draw("sames");
		pave->Draw();
		legBkg->Draw();
		
		//2D correlation dMdpT
		
		ClearSelections(hSparseDeriv);
		TH2D *hdMdPtDeriv = Projections2DWithCuts(hSparseDeriv, axprojFluc, mincut, maxcut, Form("hdMdPtDerivRho%d", ipt));
		
		ClearSelections(hSparseConst);
		TH2D *hdMdPtConst = Projections2DWithCuts(hSparseConst, axprojFluc, mincut, maxcut, Form("hdMdPtConstRho%d", ipt));
		
		ClearSelections(hSparseConst);
		TH2D *hdMdPtNoBkg = Projections2DWithCuts(hSparseConst, axprojFlucNoBkg, mincut, maxcut, Form("hdMdPtoBkgRho%d", ipt));
		
		hdMdPtDeriv->GetZaxis()->SetRangeUser(10, 1e5);
		hdMdPtConst->GetZaxis()->SetRangeUser(10, 1e5);
		hdMdPtNoBkg->GetZaxis()->SetRangeUser(10, 1e5);
		
		hdMdPtDeriv->GetXaxis()->SetRangeUser(-10, 20);
		hdMdPtConst->GetXaxis()->SetRangeUser(-10, 20);
		hdMdPtNoBkg->GetXaxis()->SetRangeUser(-10, 20);
		hdMdPtDeriv->GetYaxis()->SetRangeUser(-10, 20);
		hdMdPtConst->GetYaxis()->SetRangeUser(-10, 20);
		hdMdPtNoBkg->GetYaxis()->SetRangeUser(-10, 20);
		
		cdMdPtRhoBinsDeriv->cd(ipt+1);
		gPad->SetLogz();
		hdMdPtDeriv-> Draw("colz");
		pave->Draw();
		cdMdPtRhoBinsConst->cd(ipt+1);
		gPad->SetLogz();
		hdMdPtConst-> Draw("colz");
		pave->Draw();
		cdMdPtRhoBinsNoBkg->cd(ipt+1);
		gPad->SetLogz();
		hdMdPtNoBkg-> Draw("colz");
		pave->Draw();
		
		if(ipt == binFig){
			TPaveText *pvDer = new TPaveText(0.48, 0.2, 0.89, 0.3, "NDC");
			pvDer->SetFillStyle(0);
			pvDer->SetBorderSize(0);
			TPaveText *pvCons = (TPaveText*)pvDer->Clone("Const");
			TPaveText *pvNoSub = (TPaveText*)pvDer->Clone("NoSub");
			pvDer->AddText("Derivative subtraction");
			pvCons->AddText("Constituent subtraction");
			pvNoSub->AddText("No background subtraction");
			
			cdMdPtRhoBinsDerivFig->cd();
			gPad->SetLogz();
			hdMdPtDeriv-> Draw("colz");
			pave->Draw();
			pvDer->Draw();
			DrawLogo(3, 0.45, 0.7, 0.85, 0.8, "", 42, "");
			
			cdMdPtRhoBinsConstFig->cd();
			gPad->SetLogz();
			hdMdPtConst-> Draw("colz");
			pave->Draw();
			pvCons->Draw();
			DrawLogo(3, 0.45, 0.7, 0.85, 0.8, "", 42, "");
			
			cdMdPtRhoBinsNoBkgFig->cd();
			gPad->SetLogz();
			hdMdPtNoBkg-> Draw("colz");
			pave->Draw();
			pvNoSub->Draw();
			DrawLogo(3, 0.45, 0.7, 0.85, 0.8, "", 42, "");
		}
	}
	
	
	
	SaveCv(cdPtRhoBins);
	SaveCv(cdMRhoBins);
	SaveCv(cdPtRatiosRhoBins);
	SaveCv(cdMRatiosRhoBins);
	SaveCv(cdMdPtRhoBinsDeriv);
	SaveCv(cdMdPtRhoBinsConst);
	SaveCv(cdMdPtRhoBinsNoBkg);
	SaveCv(cdMdPtRhoBinsDerivFig);
	SaveCv(cdMdPtRhoBinsConstFig);
	SaveCv(cdMdPtRhoBinsNoBkgFig);
	
}

// return 2D projection of a ThnSparse
TH2D* Projections2DWithCuts(THnSparseF* hnsp, Int_t *axproj, Double_t *mincut, Double_t *maxcut, TString projname){
	
	//assuming 2D projections, so axproj is a two 2 array
	// mincut and maxcut have dimension naxes. The position corresponding to the projection axes are not considered
	const Int_t n = 2;
	Int_t naxes = hnsp->GetNdimensions();
	for(Int_t iax = 0; iax<naxes; iax++){
		
		if(iax == axproj[0] || iax == axproj[1]) continue;
		if(mincut[iax] == -99) continue;
		Int_t binrange[2] = {hnsp->GetAxis(iax)->FindBin(mincut[iax]), hnsp->GetAxis(iax)->FindBin(maxcut[iax])};
		Printf("(%f, %f) = (%d, %d)", mincut[iax], maxcut[iax], binrange[0], binrange[1]);
		
		hnsp->GetAxis(iax)->SetRange(binrange[0], binrange[1]);
		
	}
	
	TH2D *hproj = hnsp->Projection(axproj[0], axproj[1]);
	hproj->SetName(projname);

	return hproj;
	
}

// return 1D projection of a ThnSparse
TH1D* Projections1DWithCuts(THnSparseF* hnsp, Int_t axproj, Double_t *mincut, Double_t *maxcut, TString projname){
	
	//assuming 2D projections, so axproj is a two 2 array
	// mincut and maxcut have dimension naxes. The position corresponding to the projection axes are not considered
	const Int_t n = 1;
	Int_t naxes = hnsp->GetNdimensions();
	for(Int_t iax = 0; iax<naxes; iax++){
		
		if(iax == axproj) continue;
		if(mincut[iax] == -99) continue;
		Int_t binrange[2] = {hnsp->GetAxis(iax)->FindBin(mincut[iax]), hnsp->GetAxis(iax)->FindBin(maxcut[iax])};
		Printf("(%f, %f) = (%d, %d)", mincut[iax], maxcut[iax], binrange[0], binrange[1]);
		
		hnsp->GetAxis(iax)->SetRange(binrange[0], binrange[1]);
		
	}
	
	TH1D *hproj = hnsp->Projection(axproj);
	hproj->SetName(projname);
	hproj->SetLineWidth(2);
	
	return hproj;
	
}

//reset all selections
void ClearSelections(THnSparseF *hnsp){
	Printf("Cleaning up range selectios!");
	Int_t naxes = hnsp->GetNdimensions();
	for(Int_t iax = 0; iax<naxes; iax++){
		hnsp->GetAxis(iax)->SetRange(0, -1);
	}

}
// read THnSparse from input file
THnSparseF* ReadInput(TString file, TString listname, TString hname){
	TList *list = ReadFile(file, listname);
	
	if(!list){
		Printf("List not found, return 0");
		return 0x0;
	}
	
	THnSparseF *hsp = (THnSparseF*)list->FindObject(hname);
	return hsp;
}


