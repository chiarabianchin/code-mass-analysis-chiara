#ifndef CreateRooUnfoldResponse_C
#define CreateRooUnfoldResponse_C
#include <TH1D.h>
#include <THnSparse.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TParameter.h>
#include <TTree.h>
#include <TRandom3.h>

#include "RooUnfoldBayes.h"

#include "/data/Work/code-mass-analysis-chiara/utils/CommonTools.C"
#include "/data/macros/LoadALICEFigures.C"

void CreateRooUnfoldResponseVarWidth(const Int_t nbinsPt, Double_t ptlims[], const Int_t nbinsM, Double_t mlims[], Double_t mTwidth, Double_t minMT, Double_t maxMT, Double_t ptTwidth, Double_t minPtT, Double_t maxPtT, TString strIn = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root",  TString hname = "fhResponseFinal", TString suff = "", Bool_t fillmiss = kTRUE);

void CreateRooUnfoldResponseVarWidthFromFile(TString fileData = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrMB.root", Double_t massbw = 2, Double_t massminT = 0, Double_t massmaxT = 20, Double_t ptbw = 20, Double_t ptminT = 20., Double_t ptmaxT = 100., TString strIn = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", TString hname = "fhResponseFinal", TString suff = "", Int_t idx = 1, Bool_t fillmiss = kTRUE);

void CreateRooUnfoldResponse(TString strIn = "AnalysisResults.root", TString strL = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", TString tag = "DerivPart", Int_t binWidthPt = 5., Int_t skipBins = 0, Double_t pt_min = -40., Double_t m_min = -20., Double_t pt_minT = 0., Double_t m_minT = 0., Double_t pt_max = 150., Double_t m_max = 40., Double_t pt_maxT = 150., Double_t m_maxT = 40., Int_t colType = -1, Bool_t fillmiss = kTRUE);

TTree* ResponseToTree(Int_t nevents, THnSparse *hresp);
void CompareTreeAndResp(TTree *tree, THnSparse *hresp);
TTree* ReadTree(Long64_t& nEv,TString filetree = "TreeResponse.root", TString trname = "treeresp");

// variable width for the rec level and fixed width for the gen level
void CreateRooUnfoldResponseVarWidthMB(TString fileData = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrMB.root"){
	
	/*
	0, ibm = 0 -> lim 0.000000
1, ibm = 1 -> lim 1.000000
2, ibm = 2 -> lim 2.000000
3, ibm = 3 -> lim 4.000000
4, ibm = 4 -> lim 6.000000
5, ibm = 5 -> lim 8.000000
6, ibm = 6 -> lim 10.000000
7, ibm = 7 -> lim 12.000000
Mass Systematic down
0, ibm = 0 -> lim 0.000000
1, ibm = 1 -> lim 1.000000
2, ibm = 2 -> lim 2.000000
3, ibm = 3 -> lim 4.000000
4, ibm = 4 -> lim 6.000000
5, ibm = 5 -> lim 8.000000
6, ibm = 6 -> lim 10.000000
Mass Systematic up
Here binSyUMado = 0
0, ibm = 0 -> lim 0.000000
1, ibm = 1 -> lim 1.000000
2, ibm = 2 -> lim 2.000000
3, ibm = 3 -> lim 4.000000
4, ibm = 4 -> lim 6.000000
5, ibm = 5 -> lim 8.000000
6, ibm = 6 -> lim 10.000000
7, ibm = 7 -> lim 12.000000
8, ibm = 8 -> lim 14.000000
pT default
4, ibm = 0 -> lim 20.000000
5, ibm = 1 -> lim 25.000000
6, ibm = 2 -> lim 30.000000
7, ibm = 3 -> lim 35.000000
8, ibm = 4 -> lim 40.000000
9, ibm = 5 -> lim 45.000000
10, ibm = 6 -> lim 50.000000
11, ibm = 7 -> lim 60.000000
12, ibm = 8 -> lim 70.000000
13, ibm = 9 -> lim 80.000000
pT Systematic down
3, ibm = 0 -> lim 15.000000
4, ibm = 1 -> lim 20.000000
5, ibm = 2 -> lim 25.000000
6, ibm = 3 -> lim 30.000000
7, ibm = 4 -> lim 35.000000
8, ibm = 5 -> lim 40.000000
9, ibm = 6 -> lim 45.000000
10, ibm = 7 -> lim 50.000000
11, ibm = 8 -> lim 60.000000
12, ibm = 9 -> lim 70.000000
pT Systematic up
5, ibm = 0 -> lim 25.000000
6, ibm = 1 -> lim 30.000000
7, ibm = 2 -> lim 35.000000
8, ibm = 3 -> lim 40.000000
9, ibm = 4 -> lim 45.000000
10, ibm = 5 -> lim 50.000000
11, ibm = 6 -> lim 60.000000
12, ibm = 7 -> lim 70.000000
13, ibm = 8 -> lim 80.000000
14, ibm = 9 -> lim 100.000000

	*/
	TFile *finData = new TFile(fileData);
	if(!finData){
		Printf("File %s not found", fileData.Data());
		return;
	}
	
	const Int_t nbinsXmb = ((TParameter<Int_t>*)finData->Get("nbinsPt"))->GetVal();
	Double_t binlimsXmb[nbinsXmb+1];
	const Int_t nbinsYmb = ((TParameter<Int_t>*)finData->Get("nbinsM")) ->GetVal();
	Double_t binlimsYmb[nbinsYmb+1];
	printf("Lims pt [%d] = {", nbinsXmb);
	for(Int_t i = 0; i<nbinsXmb+1; i++){
		binlimsXmb[i] = ((TParameter<Double_t>*)finData->Get(Form("limPt%d", i)))->GetVal();
		printf("%f, ", binlimsXmb[i]);
	}
	Printf("}");
	printf("Lims M [%d] = {", nbinsYmb);
	for(Int_t i = 0; i<nbinsYmb+1; i++){
		binlimsYmb[i] = ((TParameter<Double_t>*)finData->Get(Form("limM%d", i)))->GetVal();
		printf("%f, ", binlimsYmb[i]);
	}
	Printf("}");
	
	CreateRooUnfoldResponseVarWidth(nbinsXmb, binlimsXmb, nbinsYmb, binlimsYmb, 2, 0, 20, 20, 20., 100.);
}

//_________________________________________________________________

// Variable width for the rec level and fixed width for the gen level. 
// Customizable settings for MB and EJE samples
// the fileData has to be created with DefineRangeUnfolding for MB or EJE sample (MassVsPtVarWBkgNoSubTrMB.root or MassVsPtVarWBkgNoSubTrJ1.root)
// Give the bin width and range for the generated level mass and pT, which will be with fixed binning

void CreateRooUnfoldResponseVarWidthFromFile(TString fileData, Double_t massbw, Double_t massminT, Double_t massmaxT, Double_t ptbw, Double_t ptminT, Double_t ptmaxT, TString strIn, TString hname, TString suff, Int_t idx, Bool_t fillmiss){
	
	//idx: 1 = default, 2 = SysDown; 3 = SysUp; 4 = default range fixed bins
	
	TString namenptbins = "nbinsPt";
	TString namenmabins = "nbinsM";
	TString nameptlims  = "limPt";
	TString namemalims  = "limM";
	TString sufname     = "";
	if(idx == 2) sufname= "SyD";
	if(idx == 3) sufname= "SyU";
	if(idx == 4) sufname= "Fix";
	
	TFile *finData = new TFile(fileData);
	if(!finData){
		Printf("File %s not found", fileData.Data());
		return;
	}
	Printf("Reading %s%s and %s%s ", namenptbins.Data(), sufname.Data(), namenmabins.Data(), sufname.Data());
	
	const Int_t nbinsXmb = ((TParameter<Int_t>*)finData->Get(Form("%s%s", namenptbins.Data(), sufname.Data())))->GetVal();
	Double_t binlimsXmb[nbinsXmb+1];
	const Int_t nbinsYmb = ((TParameter<Int_t>*)finData->Get(Form("%s%s", namenmabins.Data(), sufname.Data()))) ->GetVal();
	Double_t binlimsYmb[nbinsYmb+1];
	printf("Lims pt [%d] = {", nbinsXmb);
	
	Double_t deltaPt = 0, deltaM = 0;
	
	if(idx == 4){
		binlimsXmb[0] = ((TParameter<Double_t>*)finData->Get(Form("%sMin", nameptlims.Data())))->GetVal();
		deltaPt = (((TParameter<Double_t>*)finData->Get(Form("%sMax", nameptlims.Data())))->GetVal() - binlimsXmb[0])/(Double_t)nbinsXmb;
		
		binlimsYmb[0] = ((TParameter<Double_t>*)finData->Get(Form("%sMin", namemalims.Data())))->GetVal();
		deltaM = (((TParameter<Double_t>*)finData->Get(Form("%sMax", namemalims.Data())))->GetVal() - binlimsYmb[0])/(Double_t)nbinsYmb;
		
	}
	for(Int_t i = 0; i<nbinsXmb+1; i++){
		if(idx != 4){
			binlimsXmb[i] = ((TParameter<Double_t>*)finData->Get(Form("%s%s%d", nameptlims.Data(), sufname.Data(), i)))->GetVal();
		} else{
			if(i>0) binlimsXmb[i] = binlimsXmb[0]+i*deltaPt;
		}
		printf("%f, ", binlimsXmb[i]);
	}
	Printf("}");
	printf("Lims M [%d] = {", nbinsYmb);
	for(Int_t i = 0; i<nbinsYmb+1; i++){
		if(idx != 4){
			binlimsYmb[i] = ((TParameter<Double_t>*)finData->Get(Form("%s%s%d", namemalims.Data(), sufname.Data(), i)))->GetVal();
		}else {
			if(i>0) binlimsYmb[i] = binlimsYmb[0]+i*deltaM;
		}
		printf("%f, ", binlimsYmb[i]);
	}
	Printf("}");
	
	CreateRooUnfoldResponseVarWidth(nbinsXmb, binlimsXmb, nbinsYmb, binlimsYmb, massbw, massminT, massmaxT, ptbw, ptminT, ptmaxT, strIn, hname, suff, fillmiss);
	
}
/*
void GetLimits(Int_t& nbinsXmb, Double_t& binlimsXmb[], Int_t& nbinsYmb, Double_t& binlimsYmb[]){
TString namenptbins = "nbinsPt";
	TString namenmabins = "nbinsM";
	TString nameptlims  = "limPt";
	TString namemalims  = "limM";
	TString sufname     = "";
	if(idx == 2) sufname= "SyD";
	if(idx == 3) sufname= "SyU";
	if(idx == 4) sufname= "Fix";
	
	TFile *finData = new TFile(fileData);
	if(!finData){
		Printf("File %s not found", fileData.Data());
		return;
	}
	Printf("Reading %s%s and %s%s ", namenptbins.Data(), sufname.Data(), namenmabins.Data(), sufname.Data());
	
	const Int_t nbinsXmb = ((TParameter<Int_t>*)finData->Get(Form("%s%s", namenptbins.Data(), sufname.Data())))->GetVal();
	Double_t binlimsXmb[nbinsXmb+1];
	const Int_t nbinsYmb = ((TParameter<Int_t>*)finData->Get(Form("%s%s", namenmabins.Data(), sufname.Data()))) ->GetVal();
	Double_t binlimsYmb[nbinsYmb+1];
	printf("Lims pt [%d] = {", nbinsXmb);
	
	Double_t deltaPt = 0, deltaM = 0;
	
	if(idx == 4){
		binlimsXmb[0] = ((TParameter<Double_t>*)finData->Get(Form("%sMin", nameptlims.Data())))->GetVal();
		deltaPt = (((TParameter<Double_t>*)finData->Get(Form("%sMax", nameptlims.Data())))->GetVal() - binlimsXmb[0])/(Double_t)nbinsXmb;
		
		binlimsYmb[0] = ((TParameter<Double_t>*)finData->Get(Form("%sMin", namemalims.Data())))->GetVal();
		deltaM = (((TParameter<Double_t>*)finData->Get(Form("%sMax", namemalims.Data())))->GetVal() - binlimsYmb[0])/(Double_t)nbinsYmb;
		
	}
	for(Int_t i = 0; i<nbinsXmb+1; i++){
		if(idx != 4){
			binlimsXmb[i] = ((TParameter<Double_t>*)finData->Get(Form("%s%s%d", nameptlims.Data(), sufname.Data(), i)))->GetVal();
		} else{
			if(i>0) binlimsXmb[i] = binlimsXmb[0]+i*deltaPt;
		}
		printf("%f, ", binlimsXmb[i]);
	}
	Printf("}");
	printf("Lims M [%d] = {", nbinsYmb);
	for(Int_t i = 0; i<nbinsYmb+1; i++){
		if(idx != 4){
			binlimsYmb[i] = ((TParameter<Double_t>*)finData->Get(Form("%s%s%d", namemalims.Data(), sufname.Data(), i)))->GetVal();
		}else {
			if(i>0) binlimsYmb[i] = binlimsYmb[0]+i*deltaM;
		}
		printf("%f, ", binlimsYmb[i]);
	}
	Printf("}");

}
*/
//_________________________________________________________________

void CreateRooUnfoldResponseVarWidth(const Int_t nbinsPt, Double_t ptlims[], const Int_t nbinsM, Double_t mlims[], Double_t mTwidth, Double_t minMT, Double_t maxMT, Double_t ptTwidth, Double_t minPtT, Double_t maxPtT, TString strIn, TString hname, TString suff, Bool_t fillmiss){
	// Int_t triggerType = 1 MB, 2 = EJE
	TStopwatch watch;
	watch.Start();
	
	THnSparseF *hn = 0x0;
	
	
	TFile *f = new TFile(strIn.Data());
	if(!strIn.Contains("AnalysisResults")){
		Printf("Get response weighted");
		hn = static_cast<THnSparseF*>(f->Get(hname.Data()));
		
	} else {
		TString nameResponse = "fhnMassResponse";
		Int_t   ic           = 0; //centrality bin
		TList *lst = static_cast<TList*>(f->Get(hname.Data()));
		if(!lst){
			Printf("List not found");
		}
		lst->Print();
		//Get response
		hn = static_cast<THnSparseF*>(lst->FindObject(Form("%s%s_%d",nameResponse.Data(),"",ic)));
		if(!hn) hn = static_cast<THnSparseF*>(lst->FindObject(nameResponse.Data()));
		if(!hn) {
			Printf("Could not find fhnMassResponse_%d",ic);
			return;
		}
	}
	
	if(!hn) {
     	Printf("Could not find %s", hname.Data());
     	return;
    }
    TCanvas *cDet = new TCanvas("cDet", "Detector level", 800, 400);
    cDet->Divide(2,1);
    
    TCanvas *cPar = new TCanvas("cPar", "Particle level", 800, 400);
    cPar->Divide(2,1);
    
    TCanvas *cMis = new TCanvas("cMis", "Non reconstructed", 400, 400);
    
    Int_t axMDet = 0, axMPar = 1, axPtDet = 2, axPtPar = 3;
    Printf("Number of dimensions = %d, proj det level %d, %d", hn->GetNdimensions(), axMDet, axPtDet);
    TH2D* hDetLev = (TH2D*)hn->Projection(axMDet, axPtDet, "E");
    hDetLev->SetName("fh2Smear");
    TH2D* hDetRbVarW = RebinVariableWidth2D(nbinsPt, ptlims, nbinsM, mlims, hDetLev);
    hDetRbVarW->SetName("fh2RespDimM");
    hDetRbVarW->Scale(1./hDetRbVarW->Integral("width"));
    cDet->cd(1);
    hDetLev->Draw("colz");
    cDet->cd(2);
    hDetRbVarW->Draw("colz");
    
    Double_t nbinsD = (maxMT - minMT)/mTwidth;
    Printf("nbinsM = %f", nbinsD);
    
    const Int_t nbinsMT = (Int_t)nbinsD;
    Double_t limsMT[nbinsMT];
    
    nbinsD = (maxPtT - minPtT)/ptTwidth;
    Printf("nbinsPt  = %f", nbinsD);
    
    const Int_t nbinsPtT = (Int_t)nbinsD;
    Double_t limsPtT[nbinsPtT];
    
    for(Int_t ib = 0; ib < nbinsPtT+1; ib++){
    	limsPtT[ib] = minPtT + ib*ptTwidth;
    	Printf("ptT %d = %f", ib, limsPtT[ib]);
    
    }
    
    for(Int_t ib = 0; ib < nbinsMT+1; ib++){
    	limsMT[ib] = minMT + ib*mTwidth;
    
    }
    
    TH2D* hParLev = (TH2D*)hn->Projection(axMPar, axPtPar, "E");
    hParLev->SetName("fh2Prior");
    TH2D* hParRbVarW = RebinVariableWidth2D(nbinsPtT, limsPtT, nbinsMT, limsMT, hParLev);
    hParRbVarW->SetName("fh2RespDimT");
    hParRbVarW->Scale(1./hParRbVarW->Integral("width"));
    cPar->cd(1);
    hParLev->Draw("colz");
    cPar->cd(2);
    hParRbVarW->Draw("colz");

    TH2D *fh2Miss = (TH2D*)hParRbVarW->Clone("fh2Miss");
    fh2Miss->Reset();
    
    RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp","resp");
    fResponse->Setup(hDetRbVarW, hParRbVarW);
    
    //loop over the unrebinned ThnSparse, sum the content of the bins corresponding to the range I want and fill the response
    //Int_t bins[4] = {nbinsM, nbinsMT, nbinsPt, nbinsPtT };
    //Double_t limsmin[4] = {mlims[0], limsMT[0], ptlims[0], limsPtT[0]};
    //Double_t limsmax[4] = {mlims[0], limsMT[0], ptlims[0], limsPtT[0]};
    //THnSparseF *hresp = new THnSparseF("hresp", 4, bins, limsmin, limsmax);
    Int_t nDim=4;
    Int_t nbin = hn->GetNbins();
    Int_t* coord = new Int_t[nDim];
    for(Int_t bin=0; bin<nbin; bin++) {
    	Double_t w = hn->GetBinContent(bin,coord);
    	Double_t msub = hn->GetAxis(0)->GetBinCenter(coord[0]);
    	Double_t mtru = hn->GetAxis(1)->GetBinCenter(coord[1]); 
    	Double_t ptsub = hn->GetAxis(2)->GetBinCenter(coord[2]); 
    	Double_t pttru = hn->GetAxis(3)->GetBinCenter(coord[3]);
    	if(msub>=mlims[0] && msub<mlims[nbinsM]
    		&& mtru>=limsMT[0] && mtru<limsMT[nbinsMT]
    		&& ptsub>=ptlims[0] && ptsub<ptlims[nbinsPt]
    		&& pttru>=limsPtT[0] && pttru<limsPtT[nbinsPtT]
    		)
    	fResponse->Fill(ptsub,msub,pttru,mtru,w); // filling the response in the requested ranges
    	else {
    		//outside the ranges fill the missed events in the true matrix
    		if(fillmiss){
    			fResponse->Miss(pttru,mtru,w); 
    			fh2Miss->Fill(pttru,mtru,w);
    		}
    	}
    }
    /*
    for(Int_t i = 0; i < nbinsM; i++){ //axMDet
    	Int_t binrangeOrigM[2] = {hn->GetAxis(axMDet)->FindBin(mlims[i]), hn->GetAxis(axMDet)->FindBin(mlims[i+1])-1};
    	
    	for(Int_t j = 0; j < nbinsMT; j++){ //axMPar
    		Int_t binrangeOrigMT[2] = {hn->GetAxis(axMPar)->FindBin(limsMT[j]), hn->GetAxis(axMPar)->FindBin(limsMT[j+1])-1};
    		
    		for(Int_t k = 0; k < nbinsPt; k++){ //axPtDet
    			Int_t binrangeOrigPt[2] = {hn->GetAxis(axPtDet)->FindBin(ptlims[k]), hn->GetAxis(axPtDet)->FindBin(ptlims[k+1])-1};
    			
    			for(Int_t l = 0; l < nbinsPtT; l++){ //axPtPar
    				Int_t binrangeOrigPtT[2] = {hn->GetAxis(axPtPar)->FindBin(limsPtT[l]), hn->GetAxis(axPtPar)->FindBin(limsPtT[l+1])-1};
    				
    				Int_t binMin[4] = {binrangeOrigM[0], binrangeOrigMT[0], binrangeOrigPt[0], binrangeOrigPtT[0]};
    				Int_t binMax[4] = {binrangeOrigM[1], binrangeOrigMT[1], binrangeOrigPt[1], binrangeOrigPtT[1]};
    				
    				//Printf("Bin M %d-%d, binMT %d-%d, Bin Pt %d-%d, binPtT %d-%d", binrangeOrigM[0], binrangeOrigM[1], binrangeOrigMT[0], binrangeOrigMT[1], binrangeOrigPt[0], binrangeOrigPt[1], binrangeOrigPtT[0], binrangeOrigPtT[1]);
    				Long64_t globBinMin = hn->GetBin(binMin), globBinMax = hn->GetBin(binMax);
    				//Printf("Integrating gl bins %lld  %lld", globBinMin, globBinMax);
    				Double_t integralInRange = 0;
    				for (Long64_t ib = globBinMin; ib <= globBinMax; ib++) {
    					Double_t v = hn->GetBinContent(ib);
    					integralInRange += v;
    				}
    				//Printf(" = %f", integralInRange);
    				Printf("fResponse->Fill(%f, %f, %f, %f, %e)", (mlims[i+1] + mlims[i])*0.5, (limsMT[i+1] + limsMT[i])*0.5, (ptlims[i+1] + ptlims[i])*0.5, (limsPtT[i+1] + limsPtT[i])*0.5, integralInRange);
    				fResponse->Fill((ptlims[i+1] + ptlims[i])*0.5, (mlims[i+1] + mlims[i])*0.5,  (limsPtT[i+1] + limsPtT[i])*0.5, (limsMT[i+1] + limsMT[i])*0.5, integralInRange);
    				
    			}
    		}
    	}
	}
	
	Long64_t nglbins = hn->GetNbins();
	Int_t* coord = new Int_t[4];
	for(Long64_t ib = 0; ib<nglbins; ib++){
		Double_t w = hn->GetBinContent(ib,coord);
		Double_t msub  = hn->GetAxis(axMDet) ->GetBinCenter(coord[0]);
		Double_t mtru  = hn->GetAxis(axMPar) ->GetBinCenter(coord[1]); 
		Double_t ptsub = hn->GetAxis(axPtDet)->GetBinCenter(coord[2]); 
		Double_t pttru = hn->GetAxis(axPtPar)->GetBinCenter(coord[3]);
		if(msub<mlims[0] || msub>mlims[nbinsM] || mtru < limsMT[0] || mtru > limsMT[nbinsMT] || ptsub < ptlims[0] || ptsub > ptlims[nbinsPt] || pttru < limsPtT[0] || pttru > limsPtT[nbinsPtT]){
			fResponse->Miss(pttru,mtru,w);
			fh2Miss->Fill(pttru,mtru,w);
		}
	}
	
	*/
	cMis->cd();
	if(fillmiss) fh2Miss->Draw("colz");
	else delete cMis;
	
	//Write response + 2D histos to file
	TFile *fout = new TFile(Form("response%s.root",suff.Data()),"RECREATE");
	hn->Write();
	fResponse->Write("resp");
	hParRbVarW->Write();
	hDetRbVarW->Write();
	hParLev->Write();
	hDetLev->Write();
	if(fillmiss) fh2Miss->Write();
	fout->Write();
	fout->Close();
	
	
	//delete hParRbVarW ;
	//delete hDetRbVarW ;
	//delete hParLev    ;
	//delete hDetLev    ;
	//delete fResponse  ;
	
	watch.Stop();
	watch.Print();
	Printf("N bins reco: %d x %d = %d", nbinsM, nbinsPt, nbinsM*nbinsPt);
	Printf("N bins gene: %d x %d = %d", nbinsMT, nbinsPtT, nbinsMT*nbinsPtT);
}

//______________________________________________________________________________

void CreateRooUnfoldResponseReadingLimits(TString fileData = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrMB.root", TString strIn = "AnalysisResults.root", TString strL = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", TString tag = "DerivPart"){
	
	TFile *finData = new TFile(fileData);
	if(!finData){
		Printf("File %s not found", fileData.Data());
		return;
	}
	
	
	TH2D *h2dData = (TH2D*)finData->Get("hMeasStdRange");
	if(!h2dData){
		Printf("File of data not found");
		finData->ls();
		return;
	}
	//X = pt, Y = M	
	const Int_t nbinsXmb = h2dData->GetNbinsX();
	Double_t binlimsXmb[nbinsXmb+1];
	const Int_t nbinsYmb = h2dData->GetNbinsY();
	Double_t binlimsYmb[nbinsYmb+1];

	Double_t minX = h2dData->GetXaxis()->GetBinLowEdge(1), maxX = h2dData->GetXaxis()->GetBinLowEdge(nbinsXmb+1);
	Double_t minY = h2dData->GetYaxis()->GetBinLowEdge(1), maxY = h2dData->GetYaxis()->GetBinLowEdge(nbinsYmb+1);
	
	Printf("Lims pt [%d] = {%.0f, %.0f}", nbinsXmb, minX, maxX);
	
	Printf("Lims M [%d] = {%.0f, %.0f}", nbinsYmb, minY, maxY);
	
	Double_t pt_minT = 0., pt_maxT = 150.;
	Double_t m_minT = 0., m_maxT = 40.;
	CreateRooUnfoldResponse(strIn, strL, tag, h2dData->GetXaxis()->GetBinWidth(2), 0, minX, minY, pt_minT, m_minT, maxX, maxY, pt_maxT, m_maxT, -1);
	
	return;
}

//______________________________________________________________________________

void CreateRooUnfoldResponse(TString strIn, TString strL, TString tag, Int_t binWidthPt, Int_t skipBins, Double_t pt_min, Double_t m_min, Double_t pt_minT, Double_t m_minT, Double_t pt_max, Double_t m_max, Double_t pt_maxT, Double_t m_maxT, Int_t colType, Bool_t fillmiss) {
	
	Printf("Miss is %sfilled", fillmiss ? "" : "NOT");
  //The response matrix will determine in which binning the data will be evaluated
  //RooUnfold needs to be installed
  // pt_min = minimum data pT
  // m_min = minimum data M
  // pt_minT = minimum response pT
  // m_minT = minimum response M
  //colType
  //-1: use default settings (MB data sets)
  // 1: triggered p-Pb data set ptMinMeas=60
  // 2: triggered p-Pb data set ptMinMeas=40

  //RooUnfold library: set path to RooUnfold as environment variable
  //gSystem->Load("$ROOUNFOLD/libRooUnfold.so");

  //strIn: root file containing response matrix
  //strL: name of list from which to take response matrix. Default name is fhnMassResponse*
  //tag: suffix for root file which will be created
  //binWidthPt: default 5 GeV. Can be changed to other value but needs to be possible with binning of original response in THnSparse
  //skipBins: exclude first pt bin of response
  //pt_min: minimum jet pT on detector-level axis
  //m_min: minimum jet mass on detector-level axis

  TString nameResponse = "fhnMassResponse";
  Int_t   ic           = 0; //centrality bin

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->LoadMacro("$gitJetMass/utils/plotUtils.C");
  //gROOT->LoadMacro("$gitJetMass/utils/style.C");
  //gROOT->LoadMacro("$gitJetMass/utils/utilsMV.C");
  //SetStyle(0);
  //gStyle->SetMarkerSize(gStyle->GetMarkerSize()*1.25);
  TStopwatch watch;
  watch.Start();
  
  THnSparseF *hn = 0x0;
  
  TFile *f = new TFile(strIn.Data());
  f->ls();
  

  if(!strIn.Contains("AnalysisResults")){
     Printf("Get response weighted");
     hn = static_cast<THnSparseF*>(f->Get(strL.Data()));
     if(!hn) {
     	Printf("Could not find %s", strL.Data());
     	return;
     }
  }else {
     TList *lst = static_cast<TList*>(f->Get(strL.Data()));
       if(!lst){
       	  Printf("List not found");
       }
     lst->Print();
     //Get response
     hn = static_cast<THnSparseF*>(lst->FindObject(Form("%s%s_%d",nameResponse.Data(),"",ic)));
     if(!hn) hn = static_cast<THnSparseF*>(lst->FindObject(nameResponse.Data()));
     if(!hn) {
     	Printf("Could not find fhnMassResponse_%d",ic);
     	return;
     }
  }
  
  hn->Sumw2();
  Double_t scale = 1e8;
  hn->Scale(scale);

  const Int_t ndim = 4;
  Int_t dim[ndim];
  for(Int_t i = 0; i<ndim; i++)
    dim[i] = i;

  THnSparse *hnRebin = static_cast<THnSparse*>(hn->Rebin(1));
  THnSparse *fhnSparseReduced = static_cast<THnSparse*>(hnRebin->Projection(ndim,dim,"E"));

  Int_t nDim = fhnSparseReduced->GetNdimensions();
  Int_t iMSub   = 0;
  Int_t iMTrue  = 1;
  Int_t iPtSub  = 2;
  Int_t iPtTrue = 3;

  TH2D *fh2Smear = dynamic_cast<TH2D*>(fhnSparseReduced->Projection(iMSub,iPtSub,"E")); //projection on pTrec (x), Mrec (y)
  TH2D *fh2Prior = dynamic_cast<TH2D*>(fhnSparseReduced->Projection(iMTrue,iPtTrue,"E")); //projection on pTpar (x), Mpar (y)

  // 2 is the number of 2D histograms, 0 is det, 1 is true
  Int_t nBinPt[2];// = {36,30};//meas(5),true
  Double_t ptmin[2] = {pt_min,pt_minT};
  Double_t ptmax[2] = {pt_max,pt_maxT}; //{140.,150.};
  Int_t nBinM[2];// = {30,20};//meas(2),true(2)
  Double_t mmin[2] = {m_min,m_minT};
  Double_t mmax[2] = {m_max,m_maxT}; // {40.,40.};

  //change binning for p-Pb triggered data set
  if(colType==1) {
    ptmin[0] = 60.; //meas
    ptmin[1] = 40.; //gen
  }
  else if(colType==2) {
    ptmin[0] = 40.; //meas
    ptmin[1] = 40.; //gen
  }

  Double_t binWidthM = 2.;
  for(Int_t i = 0; i<2; i++) {
    nBinPt[i] = TMath::Nint((ptmax[i]-ptmin[i])/binWidthPt);
    nBinM[i] = TMath::Nint((mmax[i]-mmin[i])/binWidthM);
  }

  if(skipBins>0) {
    nBinPt[1]--;
    ptmin[1]+=binWidthPt;
  }

  //dimensions of measured axis
  TH2D *fh2RespDimM = new TH2D("fh2RespDimM","fh2RespDimM",nBinPt[0],ptmin[0],ptmax[0],nBinM[0],mmin[0],mmax[0]);
  //dimensions of true axis
  TH2D *fh2RespDimT = new TH2D("fh2RespDimT","fh2RespDimT",nBinPt[1],ptmin[1],ptmax[1],nBinM[1],mmin[1],mmax[1]);
  //feed-out of response
  TH2D *fh2Miss     = new TH2D("fh2Miss","fh2Miss",nBinPt[1],ptmin[1],ptmax[1],nBinM[1],mmin[1],mmax[1]);
  Printf( "fh2Smear->GetEntries() %f ", fh2Smear->GetEntries());
  
  //fill detector-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimM->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimM->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimM->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Smear->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Smear->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimM->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimM->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimM->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Smear->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Smear->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Smear->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      fh2RespDimM->SetBinContent(ix,iy,con);
      fh2RespDimM->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimM");

  //fill particle-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimT->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimT->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimT->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Prior->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Prior->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimT->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimT->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimT->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Prior->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Prior->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Prior->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      fh2RespDimT->SetBinContent(ix,iy,con);
      fh2RespDimT->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimT");

  //response object for RooUnfold
  // the two histograms are like hSmear ( = rec) and hPrior ( = part), but limited at the axis ranges wanted 
  RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp","resp");
  fResponse->Setup(fh2RespDimM,fh2RespDimT);

  //Fill RooUnfoldResponse object
  Int_t* coord = new Int_t[nDim];
  Int_t nbin = fhnSparseReduced->GetNbins();
  for(Int_t bin=0; bin<nbin; bin++) {
  	  Double_t w = fhnSparseReduced->GetBinContent(bin,coord);
  	  Double_t msub = fhnSparseReduced->GetAxis(0)->GetBinCenter(coord[0]);
  	  Double_t mtru = fhnSparseReduced->GetAxis(1)->GetBinCenter(coord[1]); 
  	  Double_t ptsub = fhnSparseReduced->GetAxis(2)->GetBinCenter(coord[2]); 
  	  Double_t pttru = fhnSparseReduced->GetAxis(3)->GetBinCenter(coord[3]);
  	  if(msub>=mmin[0] && msub<mmax[0]
  	  	  && mtru>=mmin[1] && mtru<mmax[1]
  	  	  && ptsub>=ptmin[0] && ptsub<ptmax[0]
  	  	  && pttru>=ptmin[1] && pttru<ptmax[1]
  	  	  )
      fResponse->Fill(ptsub,msub,pttru,mtru,w); // filling the response in the requested ranges
      else {
      	  //outside the ranges fill the missed events in the true matrix
      	  if(fillmiss){
      	  	  fResponse->Miss(pttru,mtru,w); 
      	  	  fh2Miss->Fill(pttru,mtru,w);
      	  }
      }
  }
  
  delete [] coord;

  //Write response + 2D histos to file
  TFile *fout = new TFile(Form("response%s.root",tag.Data()),"RECREATE");
  hn->Write("fhn");
  fhnSparseReduced->Write("fhnReduced");
  fResponse->Write("resp");
  fh2Smear->Write("fh2Smear");
  fh2Prior->Write("fh2Prior");
  fh2RespDimM->Write();
  fh2RespDimT->Write();
  if(fillmiss) fh2Miss->Write();
  fout->Write();
  fout->Close();
  
  
  delete fh2RespDimM ;
  delete fh2RespDimT ;
  delete fh2Miss     ;
  delete fResponse;
  
  watch.Stop();
  watch.Print();
}

//_________________________________________________________________

void DefineRangeWithEnoughStatistics(TString strDat = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/InputDataMB/MassOutput.root", TString hname2D = "h2MPtTagged_BkgSub2_TrgCmb0", Double_t finBWx = 5., Double_t finBWy = 2., Double_t mincounts = 10){
	
	// hname2D for the MC response fhResponseFinal_proj_3_1
	
	TFile *fin = new TFile(strDat);
	if(!fin->IsOpen()){
		Printf("%s not found", strDat.Data());
		return;
	}
	
	TH2D* hMassPt = (TH2D*)fin->Get(hname2D);
	
	Int_t rebinx = finBWx/hMassPt->GetXaxis()->GetBinWidth(1);
	Int_t rebiny = finBWy/hMassPt->GetYaxis()->GetBinWidth(1);
	Printf("Rebin x -> %d, Rebin y -> %d", rebinx, rebiny);
	TH2* hMassPtRb = hMassPt->Rebin2D(rebinx, rebiny);
	TH2D* hInvertedContent = (TH2D*)hMassPtRb->Clone("hInvertedContent");
	hInvertedContent->Reset();
	
	for(Int_t ibx = 0; ibx < hMassPtRb->GetNbinsX(); ibx++){
		for(Int_t iby = 0; iby < hMassPtRb->GetNbinsY(); iby++){
			Double_t content = hMassPtRb->GetBinContent(ibx, iby);
			if(content < mincounts){
				hInvertedContent->SetBinContent(ibx+1, iby+1, mincounts);
				
			}
		}
	}
	
	TCanvas *cRanges = new TCanvas(Form("c%sRanges_min%.0f", hname2D.Data(), mincounts) , Form("Usable range (content > %.0f)", mincounts), 600, 600);
	cRanges->cd();
	hInvertedContent->Draw("colz");
	
	SaveCv(cRanges);

}

//_________________________________________________________________

void ResponseToTree(Int_t nevents, TString filename, TString hname = "fhResponseFinal"){
	
	TFile *fin = new TFile(filename);
	if(!fin->IsOpen()){
		Printf("File not found");
		return;
	}
	
	THnSparse *hresp = (THnSparse*) fin->Get(hname);
	
	TFile* fout = new TFile("TreeResponse.root", "recreate");
	
	TTree* tree = ResponseToTree(nevents, hresp);
	if(!tree) {
		Printf("No output tree");
		return;
	}

	CompareTreeAndResp(tree, hresp);
	
	tree->Write();
	fout->Close();
	Printf("Tree written in root file");
	

	return;
}
//_________________________________________________________________

void CompareTreeAndRespFromFile(TString filename, TString filetree, TString hname = "fhResponseFinal", TString trname = "treeresp"){
	TFile *fin = new TFile(filename);
	if(!fin->IsOpen()){
		Printf("File not found");
		return;
	}
	
	THnSparse *hresp = (THnSparse*) fin->Get(hname);
	
	TFile *fin2 = new TFile(filetree);
	if(!fin2->IsOpen()){
		Printf("File not found");
		return;
	}
	
	TTree *tree = (TTree*) fin2->Get(trname);
	
	if(!tree) return;
	if(!hresp) return;

	CompareTreeAndResp(tree, hresp);

}
//_________________________________________________________________

TTree* ResponseToTree(Int_t nevents, THnSparse *hresp){
	TStopwatch watch;
	watch.Start();
	if(!hresp) return 0x0;
	gRandom = new TRandom3();
	TTree *tree = new TTree("treeresp", "Response tree");
	
	
	Double_t mt, pTt, mr, pTr, weight;
	
	tree->Branch("mr", &mr);
	tree->Branch("mt", &mt);
	tree->Branch("pTr", &pTr);
	tree->Branch("pTt", &pTt);
	tree->Branch("weight", &weight);
	Int_t ndim = hresp->GetNdimensions();
	if(ndim!= 4) {
		Printf("The dimentions considered are the first 4");
	}
	
	Double_t integral = hresp->ComputeIntegral();
	Printf("Integral of the response %f", integral);

	Double_t binw[4];
	for(Int_t iax = 0; iax<4; iax++) {
		binw[iax] = hresp->GetAxis(iax)->GetBinWidth(1);
		Printf("Bin width %d = %f (nbins = %d)", iax, binw[iax], hresp->GetAxis(iax)->GetNbins());
	}
	Double_t coord[4];
	for(Int_t i = 0; i<hresp->GetAxis(0)->GetNbins(); i++){ // mr
		coord[0] = hresp->GetAxis(0)->GetBinCenter(i+1);
		Printf("Mreco = %f", mr);
		for(Int_t j = 0; j<hresp->GetAxis(1)->GetNbins(); j++){ //mt
			coord[1] = hresp->GetAxis(1)->GetBinCenter(j+1);
			
			for(Int_t k = 0; k<hresp->GetAxis(2)->GetNbins(); k++){ //pTr
				coord[2] = hresp->GetAxis(2)->GetBinCenter(k+1);
				
				for(Int_t l = 0; l<hresp->GetAxis(3)->GetNbins(); l++){ //pTt
					coord[3] = hresp->GetAxis(3)->GetBinCenter(l+1);
					
					
					Int_t binaxes[4] = {i+1, j+1, k+1, l+1};
					
					Double_t content = hresp->GetBinContent(binaxes);
					if (content == 0) continue;
					Printf("Bin %d, %d, %d, %d", i+1, j+1, k+1, l+1);
					
					weight = (Double_t)content/integral;
					Printf("Filling with weight %e", weight);
					for(Int_t iev = 0; iev < nevents; iev++){
						//might add here each time a 
						mr  = gRandom->Uniform(coord[0]-binw[0]/2. , coord[0]+binw[0]/2.);
						mt  = gRandom->Uniform(coord[1]-binw[1]/2. , coord[1]+binw[1]/2.);
						pTr = gRandom->Uniform(coord[2]-binw[2]/2., coord[2]+binw[2]/2.);
						pTt = gRandom->Uniform(coord[3]-binw[3]/2., coord[3]+binw[3]/2.);
						
						tree->Fill();
					}
					
				}
				
			}
		}
	}
	watch.Stop();
	watch.Print();
	
	return tree;
}

//_________________________________________________________________
void CompareTreeAndResp(TTree *tree, THnSparse *hresp){
	
	Double_t mt, pTt, mr, pTr, weight;
	tree->SetBranchAddress("mr", &mr);
	tree->SetBranchAddress("mt", &mt);
	tree->SetBranchAddress("pTr", &pTr);
	tree->SetBranchAddress("pTt", &pTt);
	tree->SetBranchAddress("weight", &weight);
	TH1D** htree = new TH1D*[4];
	TH1D** hpjres= new TH1D*[4];
	TH1D** hrati = new TH1D*[4];
	TCanvas *cCompare = new TCanvas("cCompare", "Compare histo and tree", 800, 800);
	cCompare->Divide(2,2);
	TCanvas *cRatios = new TCanvas("cRatios", "Ratio tree over resp proj", 800, 800);
	cRatios->Divide(2,2);
	for(Int_t iax  = 0; iax<4; iax++){
		
		hpjres[iax] = (TH1D*) hresp->Projection(iax);
		hpjres[iax]->SetLineColor(kRed);
		hpjres[iax]->SetLineWidth(2);
		hpjres[iax]->Sumw2();
		//hpjres[iax]->Scale(1./hpjres[iax]->Integral());
		Printf("Integral pj axis %d = %f",iax, hpjres[iax]->Integral());
		cCompare->cd(iax+1);
		gPad->SetLogy();
		hpjres[iax]->DrawClone();
		htree[iax] = (TH1D*)hpjres[iax]->Clone(Form("htreeAx_%d", iax));
		htree[iax]->Reset();
		htree[iax]->SetLineColor(kBlack);
		
	}
	Long64_t nentries = tree->GetEntries();
	for(Int_t ie = 0; ie< nentries; ie++){
		tree->GetEntry(ie);
		htree[0]->Fill(mr,  weight);
		htree[1]->Fill(mt,  weight);
		htree[2]->Fill(pTr, weight);
		htree[3]->Fill(pTt, weight);
	}
	for(Int_t iax  = 0; iax<4; iax++){
		cCompare->cd(iax+1);
		//htree[iax]->Scale(1./htree[iax]->Integral());
		htree[iax]->Scale(1./nentries);
		Printf("Integral branch axis %d= %f",iax, htree[iax]->Integral());
		htree[iax]->DrawClone("sames");
		hrati[iax] = (TH1D*)htree[iax]->Clone(Form("htreeOresp_%d", iax));
		hrati[iax]->Divide(hpjres[iax]);
		cRatios->cd(iax+1);
		hrati[iax]->Draw();
	}
}

//_________________________________________________________________

TTree* ReadTree(Long64_t& nEv,TString filetree, TString trname){
	
	TFile *fin2 = new TFile(filetree);
	if(!fin2->IsOpen()){
		Printf("File not found");
		return 0x0;
	}
	
	TTree *tree = (TTree*) fin2->Get(trname);
	
	if(!tree) return 0x0;
	
	nEv = tree->GetEntries();
	
	return tree;
}
//_________________________________________________________________

void SmearedResp(TString fileData = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrMB.root", Double_t mTwidth = 2, Double_t minMT = 0, Double_t maxMT = 40, Double_t ptTwidth = 20, Double_t minPtT = 20., Double_t maxPtT = 140., TString strIn = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", TString suff = "", Int_t idx = 4, Bool_t fillmiss = kTRUE, Double_t frac = 0.05, TString filetree = "TreeResponse.root", TString trname = "treeresp"){
	
	TString namenptbins = "nbinsPt";
	TString namenmabins = "nbinsM";
	TString nameptlims  = "limPt";
	TString namemalims  = "limM";
	TString sufname     = "";
	if(idx == 2) sufname= "SyD";
	if(idx == 3) sufname= "SyU";
	if(idx == 4) sufname= "Fix";
	
	TFile *finData = new TFile(fileData);
	if(!finData){
		Printf("File %s not found", fileData.Data());
		return;
	}
	Printf("Reading %s%s and %s%s ", namenptbins.Data(), sufname.Data(), namenmabins.Data(), sufname.Data());
	
	const Int_t nbinsPt = ((TParameter<Int_t>*)finData->Get(Form("%s%s", namenptbins.Data(), sufname.Data())))->GetVal();
	Double_t ptlims[nbinsPt+1];
	const Int_t nbinsM = ((TParameter<Int_t>*)finData->Get(Form("%s%s", namenmabins.Data(), sufname.Data()))) ->GetVal();
	Double_t mlims[nbinsM+1];
	printf("Lims pt [%d] = {", nbinsPt);
	
	Double_t deltaPt = 0, deltaM = 0;
	
	if(idx == 4){
		ptlims[0] = ((TParameter<Double_t>*)finData->Get(Form("%sMin", nameptlims.Data())))->GetVal();
		deltaPt = (((TParameter<Double_t>*)finData->Get(Form("%sMax", nameptlims.Data())))->GetVal() - ptlims[0])/(Double_t)nbinsPt;
		
		mlims[0] = ((TParameter<Double_t>*)finData->Get(Form("%sMin", namemalims.Data())))->GetVal();
		deltaM = (((TParameter<Double_t>*)finData->Get(Form("%sMax", namemalims.Data())))->GetVal() - mlims[0])/(Double_t)nbinsM;
		
	}
	for(Int_t i = 0; i<nbinsPt+1; i++){
		if(idx != 4){
			ptlims[i] = ((TParameter<Double_t>*)finData->Get(Form("%s%s%d", nameptlims.Data(), sufname.Data(), i)))->GetVal();
		} else{
			if(i>0) ptlims[i] = ptlims[0]+i*deltaPt;
		}
		printf("%f, ", ptlims[i]);
	}
	Printf("}");
	printf("Lims M [%d] = {", nbinsM);
	for(Int_t i = 0; i<nbinsM+1; i++){
		if(idx != 4){
			mlims[i] = ((TParameter<Double_t>*)finData->Get(Form("%s%s%d", namemalims.Data(), sufname.Data(), i)))->GetVal();
		}else {
			if(i>0) mlims[i] = mlims[0]+i*deltaM;
		}
		printf("%f, ", mlims[i]);
	}
	Printf("}");
		
	Long64_t nEv;
	TTree* tree = ReadTree(nEv, filetree, trname);
	
	Double_t mt, pTt, mr, pTr, weight;
	tree->SetBranchAddress("mr",  &mr);
	tree->SetBranchAddress("mt",  &mt);
	tree->SetBranchAddress("pTr", &pTr);
	tree->SetBranchAddress("pTt", &pTt);
	tree->SetBranchAddress("weight", &weight);
	
	THnSparseF *hn = 0x0;
	TString hname = "fhResponseFinal";
	
	TFile *f = new TFile(strIn.Data());
	hn = static_cast<THnSparseF*>(f->Get(hname.Data()));
	if(!hn) {
     	Printf("Could not find %s", hname.Data());
     	return;
    }
    
	Int_t axMDet = 0, axMPar = 1, axPtDet = 2, axPtPar = 3;
    TH2D* hDetLev = (TH2D*)hn->Projection(axMDet, axPtDet, "E");
    hDetLev->SetName("fh2Smear");
    TH2D* hDetRbVarW = RebinVariableWidth2D(nbinsPt, ptlims, nbinsM, mlims, hDetLev);
    hDetRbVarW->SetName("fh2RespDimM");
    hDetRbVarW->Reset();
    
    Double_t nbinsD = (maxMT - minMT)/mTwidth;
    Printf("nbinsM = %f", nbinsD);
    
    const Int_t nbinsMT = (Int_t)nbinsD;
    Double_t limsMT[nbinsMT];
    
    nbinsD = (maxPtT - minPtT)/ptTwidth;
    Printf("nbinsPt  = %f", nbinsD);
    
    const Int_t nbinsPtT = (Int_t)nbinsD;
    Double_t limsPtT[nbinsPtT];
    
    for(Int_t ib = 0; ib < nbinsPtT+1; ib++){
    	limsPtT[ib] = minPtT + ib*ptTwidth;
    	Printf("ptT %d = %f", ib, limsPtT[ib]);
    
    }
    
    for(Int_t ib = 0; ib < nbinsMT+1; ib++){
    	limsMT[ib] = minMT + ib*mTwidth;
    
    }
    
    TCanvas *cDet = new TCanvas("cDet", "Detector level", 800, 400);
    cDet->Divide(2,1);
    
    TCanvas *cPar = new TCanvas("cPar", "Particle level", 800, 400);
    cPar->Divide(2,1);
    
    TCanvas *cMis = new TCanvas("cMis", "Non reconstructed", 400, 400);
    
    TH2D* hParLev = (TH2D*)hn->Projection(axMPar, axPtPar, "E");
    hParLev->SetName("fh2Prior");
    TH2D* hParRbVarW = RebinVariableWidth2D(nbinsPtT, limsPtT, nbinsMT, limsMT, hParLev);
    hParRbVarW->Reset(); //cleaning entries
    hParRbVarW->SetName("fh2RespDimT");
    
   
    
    TH2D *fh2Miss = (TH2D*)hParRbVarW->Clone("fh2Miss");
    fh2Miss->Reset();
    
    
	double xx,yy;
	gRandom = new TRandom3();
	
	RooUnfoldResponse *responsenew = new RooUnfoldResponse("resp","respnew");
	
	responsenew->Setup(hDetRbVarW, hParRbVarW);
	
	for(int iEntry=0; iEntry< nEv; iEntry++){
		tree->GetEntry(iEntry); 
		
		//I smear the true shape according to a gaussian of given width
		double mtsmear = gRandom->Gaus(mt,0.01*frac*mt);
		//then I get the corresponding bin numbers for the true pT, the smeared true shape and the detector level pT
		//int bin1=htrueptd->FindFixBin(TMath::Abs(smear));
		//int bin2=htruept->FindFixBin(ptJetMatch);
		//int bin3=hsmearpt->FindFixBin(ptJet);
		//
		//if(bin1==0 || bin2==0) continue;
		//if(ptJet>80 || ptJet<20) continue;
		//then you need to get the 2D detector level histogram for the corresponding true bins and project onto shape axis
		//for given detector level bin
		//TH1F *histomult1d=(TH1F*) histomult[bin1-1][bin2-1]->ProjectionX( "histomult1d",bin3,bin3); 
		
		//if(histomult1d->Integral(1,histomult1d->GetNbinsX()) ==0) cout<<ptJet<<" "<<ptJetMatch<<" "<<bin1<<" "<<bin2<<" "<<bin3<<endl;
		//double xx=histomult1d->GetRandom(); 
		////then you get the detector level shape randomly according to its distribution
		//if(xx<0.3) continue;
		//outside the ranges fill the missed events in the true matrix
		if(mr>=mlims[0] && mr<mlims[nbinsM]
			&& mtsmear>=limsMT[0] && mtsmear<limsMT[nbinsMT]
			&& pTr>=ptlims[0] && pTr<ptlims[nbinsPt]
			&& pTt>=limsPtT[0] && pTt<limsPtT[nbinsPtT]){
		
		responsenew->Fill(pTr,mr,pTt,mtsmear, weight);
		hDetRbVarW->Fill(pTr, mr, weight);
		hParRbVarW->Fill(pTt, mtsmear, weight);
		
			}
			else {
				if(fillmiss){
					responsenew->Miss(pTt,mtsmear, weight); 
					fh2Miss->Fill(pTt,mtsmear, weight);
				}
			}
			
    }
    
    cPar->cd(1);
    hParLev->Draw("colz");
    cPar->cd(2);
    hParRbVarW->Draw("colz");
    
    cDet->cd(1);
    hDetLev->Draw("colz");
    cDet->cd(2);
    hDetRbVarW->Draw("colz");
    
    cMis->cd();
	if(fillmiss) fh2Miss->Draw("colz");
	else delete cMis;
	
    //Write response + 2D histos to file
	TFile *fout = new TFile(Form("response%s.root",suff.Data()),"RECREATE");
	hn->Write();
	responsenew->Write("resp");
	hParRbVarW->Write();
	hDetRbVarW->Write();
	hParLev->Write();
	hDetLev->Write();
	if(fillmiss) fh2Miss->Write();
	fout->Write();
	fout->Close();
}

void ShowVariation(){
	const Int_t n = 4;
	TString files[n] = {
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/responseDetFlBkgNoSub.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/analysis/NoSub/testPriorVariation/response2percv2.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/analysis/NoSub/testPriorVariation/response5percv2.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/analysis/NoSub/testPriorVariation/response10percv2.root"
	};
	
	TString hnames[n] = {
		"fh2RespDimT",
		"fh2RespDimT",
		"fh2RespDimT",
		"fh2RespDimT"
	};
	Int_t nametype = 1;
	const Int_t nbinsrange = 3;
	const char* newname = "hMassTSmeared";
	
	Double_t rangelims[nbinsrange+1] = {60, 80, 100, 120};
	TH2D** hist = new TH2D*[n];
	for(Int_t i = 0; i<n; i++){
		
		hist[i] = 0x0;
		TFile *fin = new TFile(files[i]);
		if(!fin->IsOpen()){
			Printf("File %s not found", files[i].Data());
			continue;
		}
		Printf(" - Color %d", colors[i]);
		hist[i] = (TH2D*)fin->Get(hnames[i]);
		
	}
	
	TH1D*** hProj = DrawProjFrom2D(n, hist, 1, 0, nbinsrange, rangelims, newname, colors, markers, nametype);
}

void ShowVariationData(){
	const Int_t n = 2;
	TString files[n] = {
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrMB.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/responseDetFlBkgNoSub.root"
		
	};
	
	TString hnames[n] = {
		"hMeasStdRange",
		"fh2RespDimM"
	};
	Int_t nametype = 1;
	const Int_t nbinsrange = 2;
	const char* newname = "hMassDataAndResp";
	
	Double_t rangelims[nbinsrange+1] = {40, 60, 80};
	TH2D** hist = new TH2D*[n];
	for(Int_t i = 0; i<n; i++){
		
		hist[i] = 0x0;
		TFile *fin = new TFile(files[i]);
		if(!fin->IsOpen()){
			Printf("File %s not found", files[i].Data());
			continue;
		}
		Printf(" - Color %d", colors[i]);
		hist[i] = (TH2D*)fin->Get(hnames[i]);
		
	}
	
	TH1D*** hProj = DrawProjFrom2D(n, hist, 1, 0, nbinsrange, rangelims, newname, colors, markers, nametype);
}


#endif
