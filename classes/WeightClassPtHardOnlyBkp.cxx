#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TList.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TParameter.h>

#include </data/Work/MyCodeJetMass/classes/WeightClass.h>


ClassImp(WeightClass)

//-------------------------------------------------------------------------------------------------
WeightClass::WeightClass() : TNamed(),
fpathPythiaFile(""),
flistPythiaFile(""),
fnamehtrials("fHistTrials"),
fnamehxsec("fHistXsection"),
fnamehevents("fHistEvents"),
fResponsethnsp(""),
fDMdPtthnsp(""),
fNpthB(10),
fFirstpthB(0),
fAxesResp(),
fAxesMPt3D(),
fAxisMp(1),
fAxisMd(0),
fAxisPtp(3),
fAxisPtd(2),
fAxis2DPtp(1),
fAxis2DPtd(0),
fAxisMdEmb(-1),
fAxisPtdEmb(-1),
fxsec(),
ftrials(),
fnevts(),
fInputFileName(""),
fInputListName(""),
fInputList(0),
fhCrossSec(0),
fhTrials(0),
fhScaleF(0),
fhNormInt(0),
fhptMpar(0),
fhResponse(0),
fhDMdPt(0),
fhMPt3D(0),
fhResponseN(0),
fhDMdPtN(0),
fhptMparW(0),
fhptMdetW(0),
fhptMparN(0),
fhptMdetN(0),
fhResponseW(0),
fhDMdPtW(0),
fhMPt3DW(0),
fhResponseWFinal(0),
fhDMdPtWFinal(0),
fhMPt3DWFinal(0),
fDoNormPerInt(0),
fIs2D(0),
fh2dResponse(0),
fh2dResponseW(0),
fh2DRespFinal(0),
fCutAxis(-1),
fMinCut(-999),
fMaxCut(-999)

{
   //default constructor
   
   fhResponse      = new THnSparseF*[fNpthB];
   fhDMdPt        = new THnSparseF*[fNpthB];
   fhMPt3D         = new THnSparseF*[fNpthB];
   fhResponseW     = new THnSparseF*[fNpthB];
   fhDMdPtW       = new THnSparseF*[fNpthB];
   fhMPt3DW        = new THnSparseF*[fNpthB];
   fhResponseN     = new THnSparseF*[fNpthB];
   fhDMdPtN        = new THnSparseF*[fNpthB];
   fhptMpar        = new TH2*[fNpthB];
   fhptMparW       = new TH2*[fNpthB];
   fhptMdetW       = new TH2*[fNpthB];
   fhptMparN       = new TH2*[fNpthB];
   fhptMdetN       = new TH2*[fNpthB];
  
   fh2dResponse    = new TH2*[fNpthB];
   fh2dResponseW   = new TH2*[fNpthB];
   
   //QA
   
   fhCrossSec = new TH1F("fhCrossSec", "Cross sections; p_{T}-hard bin; cross section", fNpthB, 0.5, fNpthB+0.5);
   fhTrials   = new TH1F("fhTrials", "Number of trials; p_{T}-hard bin; N of trials"  , fNpthB, 0.5, fNpthB+0.5);
   
   Printf("Remember to set axes order with void SetOrderAxesRespose(Int_t axMp, Int_t axMd, Int_t axptp, Int_t axptd) before Initialise!! They are currently set to :  axMp = %d, axMd = %d, axptp = %d, axptd = %d", fAxisMp, fAxisMd, fAxisPtp, fAxisPtd);

}

//-------------------------------------------------------------------------------------------------
WeightClass::WeightClass(const char* name, const char* title, const Int_t npthbins) : TNamed(name, title),
fpathPythiaFile(""),
flistPythiaFile(""),
fnamehtrials("fHistTrials"),
fnamehxsec("fHistXsection"),
fnamehevents("fHistEvents"),
fResponsethnsp(""),
fDMdPtthnsp(""),
fNpthB(npthbins),
fFirstpthB(0),
fAxesResp(),
fAxesMPt3D(),
fAxisMp(1),
fAxisMd(0),
fAxisPtp(3),
fAxisPtd(2),
fAxis2DPtp(1),
fAxis2DPtd(0),
fAxisMdEmb(-1),
fAxisPtdEmb(-1),
fxsec(),
ftrials(),
fnevts(),
fInputFileName(""),
fInputListName(""),
fInputList(0),
fhCrossSec(0),
fhTrials(0),
fhScaleF(0),
fhNormInt(0),
fhptMpar(0),
fhResponse(0),
fhDMdPt(0),
fhMPt3D(0),
fhResponseN(0),
fhDMdPtN(0),
fhptMparW(0),
fhptMdetW(0),
fhptMparN(0),
fhptMdetN(0),
fhResponseW(0),
fhDMdPtW(0),
fhMPt3DW(0),
fhResponseWFinal(0),
fhDMdPtWFinal(0),
fhMPt3DWFinal(0),
fDoNormPerInt(0),
fIs2D(0),
fh2dResponse(0),
fh2dResponseW(0),
fh2DRespFinal(0),
fCutAxis(-1),
fMinCut(-999),
fMaxCut(-999)

{
   //standard constructor
   fhResponse      = new THnSparseF*[fNpthB];
   fhDMdPt        = new THnSparseF*[fNpthB];
   fhMPt3D         = new THnSparseF*[fNpthB];
   fhResponseW     = new THnSparseF*[fNpthB];
   fhDMdPtW        = new THnSparseF*[fNpthB];
   fhMPt3DW        = new THnSparseF*[fNpthB];
   fhResponseN     = new THnSparseF*[fNpthB];
   fhDMdPtN        = new THnSparseF*[fNpthB];
   fhptMpar        = new TH2*[fNpthB];
   fhptMparW       = new TH2*[fNpthB];
   fhptMdetW       = new TH2*[fNpthB];
   fhptMparN       = new TH2*[fNpthB];
   fhptMdetN       = new TH2*[fNpthB];
   
   fh2dResponse    = new TH2*[fNpthB];
   fh2dResponseW   = new TH2*[fNpthB];

   //QA
   
   fhCrossSec = new TH1F(Form("fhCrossSec%s", GetName()), "Cross sections; p_{T}-hard bin; cross section", fNpthB, 0.5, fNpthB+0.5);
   fhTrials   = new TH1F(Form("fhTrials%s", GetName()), "Number of trials; p_{T}-hard bin; N of trials"  , fNpthB, 0.5, fNpthB+0.5);
   fhScaleF   = new TH1F(Form("fhScaleF%s", GetName()), "Scale factor x-sec/nTrials; p_{T}-hard bin; x-sec/nTrials", fNpthB, 0.5, fNpthB+0.5);
   fhNormInt  = new TH1F(Form("fhNormInt%s", GetName()), "Normalization for the integral; p_{T}-hard bin; #int_{p_{T}-hard bin} d#sigma/d #it{p}_{T, par} d #it{M, par}", fNpthB, 0.5, fNpthB+0.5);
   
   Printf("Remember to set axes order with void SetOrderAxesRespose(Int_t axMp, Int_t axMd, Int_t axptp, Int_t axptd) before Initialise!! They are currently set to :  axMp = %d, axMd = %d, axptp = %d, axptd = %d", fAxisMp, fAxisMd, fAxisPtp, fAxisPtd);
}

//-------------------------------------------------------------------------------------------------

//WeightClass::WeightClass(const WeightClass &weightcl) : TNamed(),
//fpathPythiaFile(weightcl.fpathPythiaFile),
//flistPythiaFile(weightcl.flistPythiaFile),
//fnamehtrials(weightcl.fnamehtrials),
//fnamehxsec(weightcl.fnamehxsec),
//fnamehevents(weightcl.fnamehevents),
//fResponsethnsp(weightcl.fResponsethnsp),
//fNpthB(weightcl.fNpthB),
//fAxesResp(weightcl.fAxesResp),
//fAxisMp(weightcl.fAxisMp),
//fAxisMd(weightcl.fAxisMd),
//fAxisPtp(weightcl.fAxisPtp),
//fAxisPtd(weightcl.fAxisPtd),
//fxsec(weightcl.fxsec),
//ftrials(weightcl.ftrials),
//fInputFile(weightcl.fInputFile),
//fInputList(weightcl.fInputList),
//fhCrossSec(weightcl.fhCrossSec),
//fhTrials(weightcl.fhTrials),
//fhptMpar(weightcl.fhptMpar),
//fhResponse(weightcl.fhResponse),
//fhResponseN(weightcl.fhResponseN),
//fhptMparW(weightcl.fhptMparW),
//fhResponseW(weightcl.fhResponseW)
//{
//   //copy constructor
//   //,fhResponseWFinal(weightcl.fhResponseWFinal)
//}
//

//WeightClass& WeightClass::operator=(const WeightClass &weightcl){
//// to be implemented!!!!
//}
   
   
//-------------------------------------------------------------------------------------------------
void WeightClass::Initialise(TString path, TString list, TString inputFile, TString listnameInp){
   
   /// This has to be run before the loop on the pT hard bins
   /// Initialise what is needed for retrieving the weights
   /// path = base path to get the ntrial and xsec histograms to determine the weights. It will then append 0<n>/AnalysisResults.root in ReadWeight
   /// inputFile = name of the file where the response is. It is needed to use the method ReadInputDataFromList, otherwise default ""
   /// inputFile = name of the list within inputfile where the response is. It is needed to use the method ReadInputDataFromList, otherwise default ""
   
   //these refer to PYTHIA weights
   fpathPythiaFile = path;
   flistPythiaFile = list;
   
   //this refers to the input distributions
   fInputFileName = inputFile;
   fInputListName = listnameInp;
   // for the moment, fpathPythiaFile and flistPythiaFile are also used for the input
   Printf("Initialise:: Base input path: %s, Pythia weights in: %s; Response input: %s, %s", fpathPythiaFile.Data(), flistPythiaFile.Data(), fInputFileName.Data(), fInputListName.Data());
   
   /// axes in the following order : Mass_det, Mass_par, Pt_det, Pt_par
   if(!fIs2D){
      SetArrayOrderAxesResponse();
   }
   if(fAxisMdEmb != -1 && fAxisPtdEmb != -1){
      Printf("LAST TWO AXES DEFINED %d, %d", fAxisMdEmb, fAxisPtdEmb);
      SetArrayOrderAxesMPt3D();
   }
   return;
   
}
//-------------------------------------------------------------------------------------------------
Bool_t WeightClass::ReadWeight(Int_t ipthb){
   
   ///retrieve histograms with weights for this pT hard bin
   TList* listPythia = ReadFile(Form("%s/%02d/AnalysisResults.root", fpathPythiaFile.Data(), ipthb+1), flistPythiaFile);
   if(!listPythia){
      TString secondpart = "_pT0150_E_scheme_TC"; //to be changed
      //secondpart += "noOvl";
      
      listPythia = ReadFile(Form("%s/AnalysisResults.root", fpathPythiaFile.Data()), Form("%sPtH%d%s", flistPythiaFile.Data(), ipthb+1, secondpart.Data()));
      
      if(!listPythia){   
      	 Printf("Cannot proceed with pT-hard %d, weight list not found", ipthb+1);
      	 //I don't even want to continue, because I miss one bin
      	 return kFALSE;
      }
   }
   Printf("Reading BIN %d", ipthb+1);
   Int_t offset = 1;
   TH1F*     htrials = (TH1F*)listPythia->FindObject(fnamehtrials);
   TProfile* hxsec = (TProfile*)listPythia->FindObject(fnamehxsec);
   TH1F*     hevents = (TH1F*)listPythia->FindObject(fnamehevents); //not sure if needed
   if(!htrials || !hxsec || !hevents){
      Printf("Histograms for weights not found, please check file in %s%02d/, I exit", fpathPythiaFile.Data(), ipthb+1);
      listPythia->ls();
      return kFALSE;
   }
   //TCanvas *cWeightsInput = new TCanvas("cWeightsInput", "Weights Input");
   //cWeightsInput->cd();
   //hxsec->Draw();
   
   fxsec   = hxsec->GetBinContent(ipthb+offset+1);
   ftrials = htrials->GetBinContent(ipthb+offset+1);
   fnevts  = hevents->GetBinContent(ipthb+offset+1);
   if(fxsec==0 && ftrials==0){
      Printf("Warning! using integral, something is wrong with the pT hard bin number!!");
      fxsec   = hxsec->Integral();
      ftrials = htrials->Integral();
      if(fxsec==0 && ftrials==0) {
      	 Printf("Attempt failed... empty cross section and ntrials");
      	 return kFALSE;
      }
   }
   Printf("Pt-hard %d xsec = %e, trials = %f", ipthb+1, fxsec, ftrials);
   fhCrossSec->Fill(ipthb+1, fxsec);
   fhTrials  ->Fill(ipthb+1, ftrials);
   
   return kTRUE;
}

//-------------------------------------------------------------------------------------------------
Bool_t WeightClass::ReadInputDataFromList(Int_t ipthb){
   
   /// Read list of response to be weighted for this pT hard bin
   /// Input file is the output of the "JetShape" task or in general a TList that contains a the THnSparse of the response (here only the list is read)
   
   Bool_t bpThInSameFile = kFALSE;
   TString filename = Form("%s%02d/AnalysisResults.root", fInputFileName.Data(), ipthb+1);
   TFile *inputFile = new TFile(filename);
   if(!inputFile->IsOpen()){
      filename = Form("%s/AnalysisResults.root", fInputFileName.Data());
      inputFile = new TFile(filename);
      if(!inputFile->IsOpen()){
      	 Printf("Input file %s not found, return", filename.Data());
      	 return 0;
      } else bpThInSameFile = kTRUE;
   }
   //fInputList = new TList((TList*) inputFile->Get(fInputListName));
   TString namelist = fInputListName;
   
   if(bpThInSameFile){
      Printf("-> Reading lists in same file");
      TString secondpart = "_pT0150_E_scheme_TC"; //to be changed
      //secondpart += "noOvl";
      
      namelist = Form("%sPtH%d%s", fInputListName.Data(), ipthb+1, secondpart.Data());
   
   }
   
   fInputList = (TList*) inputFile->Get(namelist);
   
   if(!fInputList){
      Printf("Cannot proceed with pT-hard %d, response list %s not found", ipthb+1, namelist.Data());
      //I don't even want to continue, because I miss one bin
      inputFile->ls();
      inputFile->Close();
      return 0;
   }
   //Printf("Pointer list %s %p", fInputListName.Data(), fInputList);
   //inputFile->ls();
   //fInputList->ls();
   return 1;
}

//-------------------------------------------------------------------------------------------------
void WeightClass::ReadInputData(Int_t ipthb, THnSparseF* hspresp){

   /// read the THnSparse containing the response and projecting on the axes defined
   if(!hspresp) {
      Printf("Null input!!");
      return;
   }
   
   const Int_t nD = 4;
   fhResponse[ipthb] = (THnSparseF*)hspresp->ProjectionND(nD, fAxesResp, "E");
   fhResponse[ipthb]->SetName(Form("hResponse%s_pTh%d", GetName(), ipthb+1));
   fhResponse[ipthb]->Sumw2();
   
   for(Int_t iax = 0; iax < 4; iax++){
      Printf("axis titles %d = %s,", iax, fhResponse[ipthb]->GetAxis(iax)->GetTitle());
   }
   //after the ND projection, this is the order
   fAxisMd = 0; fAxisMp = 1; fAxisPtd = 2; fAxisPtp = 3;
   SetArrayOrderAxesResponse();
   // if the output is from embedding, also the detector level embedded particle is present
   if(fAxisMdEmb != -1 && fAxisPtdEmb != -1){
      SetAxesMPt3DWithArray();
      const Int_t n3D = 6;
      
      fhMPt3D[ipthb] = (THnSparseF*)hspresp->ProjectionND(n3D, fAxesMPt3D, "E");
      fhMPt3D[ipthb]->SetName(Form("hhMPt3D%s_pTh%d", GetName(), ipthb+1));
      fhMPt3D[ipthb]->Sumw2();
      
      for(Int_t iax = 0; iax < 4; iax++){
      	 Printf("%d", fAxesMPt3D[iax]);
      	 Printf("MPt x 3D axis titles %d = %s,", iax, fhMPt3D[ipthb]->GetAxis(iax)->GetTitle());
      }
      //after the ND projection, this is the order
      fAxisMd = 0; fAxisMp = 1; fAxisPtd = 2; fAxisPtp = 3; fAxisMdEmb = 4; fAxisPtdEmb = 5;
      SetArrayOrderAxesMPt3D();
      Printf("Here they are %d %d", fAxisMp, fAxisPtp);
   }
   
   
   return;
}

//-------------------------------------------------------------------------------------------------
void WeightClass::ReadInputData2D(Int_t ipthb, TH2* h2dresp){

   /// read the THnSparse containing the response and projecting on the axes defined
   if(!h2dresp) {
      Printf("Null input!!");
      return;
   }

   fh2dResponse[ipthb] = h2dresp;
   fh2dResponse[ipthb]->SetName(Form("hResponse2D_pTh%d", ipthb+1));
   fh2dResponse[ipthb]->Sumw2();
   return;
}

//-------------------------------------------------------------------------------------------------
Bool_t WeightClass::SetResponse(Int_t ipthb, TString hnspname){
   
   /// To be used for the output of JetShape task
   /// read the THnSparse containig the response and performs the projection on the 4D needed. The indices can be customize (tbi)
   
   ReadInputDataFromList(ipthb);
   //Printf("SetResponse Pointer list %p", fInputList);
   if(!fInputList) {
      Printf("Input list not found!");
      return 0 ;
   }
   fResponsethnsp = hnspname;
   
   
   if(fIs2D) return SetResponse2D(ipthb);
   
   
   THnSparseF *hnsp = (THnSparseF*)fInputList->FindObject(fResponsethnsp);
   if(!hnsp){
      Printf("Cannot proceed with pT-hard %d, %s not found", ipthb+1, fResponsethnsp.Data());
      //I don't even want to continue, because I miss one bin
      fInputList->ls();
      return 0;
   }
   
   if(fCutAxis > -1 && fCutAxis < 4) {
      Int_t bmin = 0, bmax = -1;
      if(fMinCut != -999) hnsp->GetAxis(fCutAxis)->FindBin(fMinCut);
      if(fMaxCut != -999) hnsp->GetAxis(fCutAxis)->FindBin(fMaxCut);
      
      hnsp->GetAxis(fCutAxis)->SetRange(bmin, bmax);
      
   }
   ReadInputData(ipthb, hnsp);
   if (fhResponse[ipthb]) return 1;
   else return 0;
}

//-------------------------------------------------------------------------------------------------
Bool_t WeightClass::SetResponse2D(Int_t ipthb){
   
   /// read the TH2 in the list with name fResponsethnsp containig the 2D response. The indices can be customize (tbi)
   
   if(!fInputList) {
      Printf("Input list not found!");
      return 0 ;
   }
   THnSparseF *hnsp = (THnSparseF*)fInputList->FindObject(fResponsethnsp);
   if(!hnsp){
      Printf("Cannot proceed with pT-hard %d, %s not found", ipthb+1, fResponsethnsp.Data());
      //I don't even want to continue, because I miss one bin
      fInputList->ls();
      return 0;
   }
   
   ReadInputData2D(ipthb, hnsp->Projection(fAxisPtp, fAxisPtd));
   if (fh2dResponse[ipthb]) return 1;
   else return 0;
}

//-------------------------------------------------------------------------------------------------
Bool_t WeightClass::WeightResponse(Int_t ipthb){
   
   /// weight response for this pT hard bin
   
   if(fIs2D) return WeightResponse2D(ipthb);
   
   if(!fhResponse[ipthb]){
      Printf("No response to be weighted");
      return 0;
   }
   
   // projections useful for QA
   Printf("What are %d and %d", fAxisMp, fAxisPtp);
   fhptMpar[ipthb]   = fhResponse[ipthb]->Projection(fAxisMp, fAxisPtp);
   fhptMpar[ipthb]->SetName(Form("fhptMpar_pTh%d", ipthb+1));
   Double_t integralptMpar = fhptMpar[ipthb]->Integral();
   Double_t integralptpar = 1;
   //TH1D *h = IntegralParPt(ipthb, 0, integralptpar);
   Printf("---->>>>>>>>> Integral 2D par = %f, integral 1D par = %f, fNevents = %f", integralptMpar, integralptpar, fnevts);
   // definition of the scale according to the input from ReadWeight
   Double_t scale = fxsec/ftrials;
   fhScaleF->Fill(ipthb+1, scale);
   
   fhResponseW[ipthb] = (THnSparseF*)fhResponse[ipthb]->Clone(Form("fhResponseW_pTh%d", ipthb+1));
   Printf("-----------> Original pointer (pth %d) %p", ipthb, fhResponseW[ipthb]);
   fhResponseW[ipthb]->Sumw2();
   //fhResponseW[ipthb]->Scale(1./integralptMpar);
   fhResponseW[ipthb]->Scale(scale);
   fhptMparW[ipthb] = fhResponseW[ipthb]->Projection(fAxisMd, fAxisPtd); // x = fAxisMd, y = fAxisPtd
   fhptMparW[ipthb]->SetName(Form("fhptMparW_pTh%d", ipthb+1));
   fhptMdetW[ipthb] = fhResponseW[ipthb]->Projection(fAxisMd, fAxisPtd); // x = fAxisMd, y = fAxisPtd
   fhptMdetW[ipthb]->SetName(Form("fhptMdetW_pTh%d", ipthb+1));
   
   //if (fDoNormPerInt) NormalizePerIntegral(ipthb);
   if (fDoNormPerInt) NormalizePerBin(ipthb);
   //if (fDoNormPerInt){
   //   if(!fListNJets || fListNJets->GetEntries() == 0){
   //   	 Printf("Set the list of N jets for the denominator with ReadNumberOfJetPerPtHBin. Will skip the normalization");
   //   } else {
   //   	 TParameter<double> *par = (TParameter<double>*)fListNJets->At(ipthb);
   //   	 if(!par){
   //   	    Printf("Par number %d not found", ipthb);
   //   	 } else {
   //   	    Double_t denominator = par->GetVal();
   //   	    //denominator = fnevts;
   //   	    fhResponseN[ipthb] = (THnSparseF*)fhResponseW[ipthb]->Clone(Form("fhResponseWNormPerJet_pthb%d", ipthb));
   //   	    fhResponseN[ipthb]->Scale(integralptpar/denominator);
   //   	    
   //   	    //particle level projection
   //   	    fhptMparN[ipthb] = fhResponseN[ipthb]->Projection(fAxisMp, fAxisPtp);
   //   	    fhptMparN[ipthb]->Sumw2();
   //   	    fhptMparN[ipthb]->SetName(Form("hptMpar_pthb%d", ipthb));
   //   	    Printf("999999999999 Scaling (%f) effect: %e to %e to %e", integralptpar/denominator, fhptMpar[ipthb]->Integral(), fhptMparW[ipthb]->Integral(), fhptMparN[ipthb]->Integral());
   //   	 }
   //   }
   //}
   else {
      // fill some histograms that for QA
      fhResponseN[ipthb] = 0;
      
   }
   //
   //fhptMparW[ipthb] = (TH2D*)fhptMpar[ipthb]->Clone(Form("fhptMparW_pTh%d", ipthb+1));
   //fhptMparW[ipthb]->Scale(scale);
   
      
   
   
   //fhNormInt->Fill(ipthb+1, integralptMparW);
   //if(integralptMparW > 0) fhResponseW[ipthb]->Scale(scale/integralptMparW);
   //if(integralptMparW > 0) fhResponseW[ipthb]->Scale(1./integralptMparW);
   //else  return 0;

   if(fAxisMdEmb != -1 && fAxisPtdEmb != -1){
      if(!fhMPt3D[ipthb]) {
      	 Printf("No MPt3D to be weighted");
      	 return 0;
      }
      
      fhMPt3DW[ipthb] = (THnSparseF*)fhMPt3D[ipthb]->Clone(Form("fhMPt3DW_pTh%d", ipthb+1));
      fhMPt3DW[ipthb]->Sumw2();
      fhMPt3DW[ipthb]->Scale(scale);
      
   }
   
   return 1;
}

//-------------------------------------------------------------------------------------------------

TH1D* WeightClass::IntegralParPt(Int_t ipthb, Int_t level, Double_t integralptpar) const {
   
   /// method that return the integral of the particle level pt distribution for the pt hard bin ipthb, creates and fills the histogram hNJetsPerPtPBinThisPtHBin with the content of each pT,par bin
   
   //prepare the histogram range and bins : this is (should be) the same for all pT-hard bins
   Printf("First bin = %d, axis ptp = %d", ipthb, fAxisPtp);
   TH1D *hPtPar = 0x0;
   if(level == 0) {
      if(!fhResponse[ipthb]) {
      	 Printf("Response not found");
      	 return 0x0;
      }
      hPtPar = fhResponse[ipthb]->Projection(fAxisPtp);
      //hPtPar = fhResponse[ipthb]->Projection(3);
   }
   if(level == 1) {
      if(!fhResponseW[ipthb]) {
      	 Printf("Response not found");
      	 return 0x0;
      }
      hPtPar = fhResponseW[ipthb]->Projection(fAxisPtp);
      //hPtPar = fhResponseW[ipthb]->Projection(3);
   }
   if(level > 1){
      Printf("Level %d not implemented!!", level);
      return 0x0;
   }
   if(!hPtPar){
      Printf("Ooops... Projection %d is null", fAxisPtp);
      return 0x0;
   }
   Int_t nbinsPtP = hPtPar->GetNbinsX();
   Double_t minPtP = hPtPar->GetXaxis()->GetBinLowEdge(1), maxPtP = hPtPar->GetXaxis()->GetBinLowEdge(nbinsPtP + 1);
   
   TH1D *hNJetsPerPtPBinThisPtHBin = new TH1D(Form("hNJetsPerPtPBin_PtH%d", ipthb), Form("Content per #it{p}_{T, par} bin, pTH bin %d; #it{p}_{T, par}; Content", ipthb), nbinsPtP, minPtP, maxPtP);
   
   for(Int_t iptp = 0; iptp<nbinsPtP; iptp++){
      hNJetsPerPtPBinThisPtHBin->SetBinContent(iptp+1, hPtPar->GetBinContent(iptp+1));
   }
   
   integralptpar = hPtPar->Integral();
   return hNJetsPerPtPBinThisPtHBin;
   
}
//-------------------------------------------------------------------------------------------------
void WeightClass::NormalizePerIntegral(Int_t ipthb){
   
   
   Int_t rebinall[4] = {2,2,2,2}; //{4,4,10,10};
   THnSparseF* fhResponseWRebin = (THnSparseF*)fhResponseW[ipthb]->Rebin(rebinall); 
   
   fhResponseN[ipthb] = (THnSparseF*)fhResponseWRebin->Clone(Form("fhResponseWNormPerJet_pthb%d", ipthb));
   
   // per each pTpar,Mpar bin, the integral of pTrec, Mrec = 1
   Int_t nBinsPPar = fhResponseN[ipthb]->GetAxis(fAxisPtp)->GetNbins();
   Int_t nBinsMPar = fhResponseN[ipthb]->GetAxis(fAxisMp)->GetNbins();
   Int_t nBinsPDet = fhResponseN[ipthb]->GetAxis(fAxisPtd)->GetNbins();
   Int_t nBinsMDet = fhResponseN[ipthb]->GetAxis(fAxisMd)->GetNbins();
   
   fhResponseN[ipthb]->Reset();
   //Printf("Dimensions %d, proj %d, %d", fhResponseWRebin->GetNdimensions(), fAxisMp, fAxisPtp);
   Printf("Starting a loop over %d Pt bins, %d Mass bins at particle level\nPer each of these bins a loop over %d Pt bins and %d M bins at detector level will be run. Might take a while", nBinsPPar, nBinsMPar, nBinsPDet, nBinsMDet);
   
   //particle level projection
   fhptMparW[ipthb] = fhResponseWRebin->Projection(fAxisMp, fAxisPtp);
   fhptMparW[ipthb]->Sumw2();
   fhptMparW[ipthb]->SetName(Form("hptMpar_pthb%d", ipthb));
   
   Int_t count = 0;
   
   //TCanvas *cResDet = new TCanvas("cResDet", "Detector response in each part bin", 1000, 1000);
   //cResDet->Divide(20, 20);
   for(Int_t ibppar = 0; ibppar < nBinsPPar; ibppar++){
      for(Int_t ibmpar = 0; ibmpar < nBinsMPar; ibmpar++){
      	 
      	 if(count % 100 == 0) Printf("Normalizing pTpar , Mpar %d, %d", ibppar, ibmpar);
      	 //project the THnSparse on pTrec, Mrec for that selection of ptpar, Mpar
      	 fhResponseWRebin->GetAxis(fAxisMp)->SetRange(ibmpar, ibmpar);
      	 fhResponseWRebin->GetAxis(fAxisPtp)->SetRange(ibppar, ibppar);
      	 
      	 //integral at particle level
      	 Double_t integralptMparW = fhptMparW[ipthb]->Integral(ibmpar, ibmpar, ibppar, ibppar);
      	 
      	 //detector level projection
      	 fhptMdetW[ipthb] = fhResponseWRebin->Projection(fAxisMd, fAxisPtd); // x = fAxisMd, y = fAxisPtd
      	 fhptMdetW[ipthb]->SetName(Form("fhptMdetWpTP%d_MP%d_pthb%d", ibppar, ibmpar, ipthb));
      	 //integral at detector level
      	 Double_t integralDet = fhptMdetW[ipthb]->Integral();
      	 
      	 //if(integralDet < 1e-13) continue;
      	 
      	 //scale the integral of the projection to one
      	 //fhptMdetW[ipthb]->Scale(1./integralDet);
      	 if(integralptMparW < 1e-13) continue;
      	 //Printf("Integral at particle level = %e", integralptMparW);
      	 //Printf("Integral at detector level = %e", integralDet);
      	 
      	 //loop on pTdet, Mdet
      	 for(Int_t ibpdet = 0; ibpdet < nBinsPDet; ibpdet++){
      	    for(Int_t ibmdet = 0; ibmdet < nBinsMDet; ibmdet++){
      	       
      	       Int_t thisbin[4] = {ibmdet, ibmpar, ibpdet, ibppar}; //same order as the original response
      	       
      	       //use the content of each bin to fill the normalized ThnSparse per each pTdet, Mdet and pTpar,Mpar
      	       Double_t contN = fhResponseWRebin->GetBinContent(thisbin); //x, y
      	       Double_t erroN = fhResponseWRebin->GetBinError(thisbin);
      	       //if(contN>1e-3) Printf("Content %f, Error %f", contN,erroN );
      	       
      	       //fhResponseN[ipthb]->SetBinContent(thisbin, contN/integralDet);
      	       //fhResponseN[ipthb]->SetBinError(thisbin,   TMath::Sqrt(1./integralDet) * erroN);
      	       fhResponseN[ipthb]->SetBinContent(thisbin, contN/integralptMparW);
      	       fhResponseN[ipthb]->SetBinError(thisbin,   TMath::Sqrt(1./integralptMparW) * erroN);
      	    }
      	 }
      	 count++;
      	 //cResDet->cd(count);
      	 //fhptMdetW[ipthb]->Draw("colz");
      }
   }
   for(Int_t ia = 0; ia < 4; ia++){
      fhResponseWRebin->GetAxis(ia)->SetRange(0, -1);
      fhResponseN[ipthb]->GetAxis(ia)->SetRange(0, -1);
   }
   
}

//-------------------------------------------------------------------------------------------------
Int_t WeightClass::CheckBinning(Int_t ipthb, TH1D* h) const {
   /// Compare the particle level binning and range of the current response with that of the input 1D histograms
   /// It returns the rebinninn factor to be applied to the response to match the binning of the input
   
   Double_t rangePPar[2] = {fhResponseW[ipthb]->GetAxis(fAxisPtp)->GetBinLowEdge(1), fhResponseW[ipthb]->GetAxis(fAxisPtp)->GetBinLowEdge(fhResponseW[ipthb]->GetAxis(fAxisPtp)->GetNbins() + 1)};
   Double_t binWidthPPar = fhResponseW[ipthb]->GetAxis(fAxisPtp)->GetBinWidth(5);
   
  Double_t rangeH[2] = {h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinLowEdge(h->GetNbinsX() + 1)};
   Double_t binWidthH = h->GetXaxis()->GetBinWidth(5);
   
   Double_t binWfr = binWidthH/binWidthPPar;
   Printf("Info: Ranges responseW = %.2f, %.2f \tRanges H = %.2f, %.2f  ", rangePPar[0], rangePPar[1], rangeH[0], rangeH[1]);
   if(TMath::Abs(binWfr - 1.) < 1e-6) return 1;
   
   if(binWfr < 1.) return (Int_t) (-1/binWfr);
   else return (Int_t) (binWfr);

}

//-------------------------------------------------------------------------------------------------
void WeightClass::NormalizePerBin(Int_t ipthb){
   
   /// This method loops over each bin of the weighted response and scale the content by the ratio of the number of jets selected, obtained from the response itself and the number found in the MC job
   /// comment To be extended
   /// tricky to treat the bins. They have to be the same number and starting from the same values. use the input from fListNJets to check it
   
   // get the numerator of the normalization factors
   TH1D *hNJetsPerPtPBinThisPtHBin = (TH1D*) fListNJets->FindObject(Form("hNJetsPerPtPBin_PtH%d", ipthb));
    if(!hNJetsPerPtPBinThisPtHBin){
       Printf("List empty, cannot proceed with normalization");
       fListNJets->ls();
       return;
    }
    
   Int_t rebPPt = CheckBinning(ipthb, hNJetsPerPtPBinThisPtHBin); 
   Printf("CheckBinning returns %d", rebPPt);
   Int_t factoM = 2, factorPt = 5;
   if(fFineBinningSlow) {
      //fine binning, takes very long (~7hours)
      factoM = 1;
      factorPt = 5;
   }

   // Define this accordingly to the fListNJets input histograms and the current response
   // Reminder: axis order: Md, Mp, pTd, pTp
   Int_t rebinall[4] = {factoM, factoM, factorPt*rebPPt , factorPt*rebPPt}; //{2,2,2,2}; //{4,4,10,10};
   if(rebPPt < 0) {
      Printf("Bin width of the input histogram is smaller than that of the response. Check that this is working");
      rebinall[2] = rebinall[3] = 1;
      rebPPt = TMath::Abs(rebPPt); 
      hNJetsPerPtPBinThisPtHBin->Rebin(rebPPt);
   }
   THnSparseF* fhResponseWRebin = (THnSparseF*)fhResponseW[ipthb]->Rebin(rebinall); 
   hNJetsPerPtPBinThisPtHBin->Rebin(factorPt);
   
   fhResponseN[ipthb] = (THnSparseF*)fhResponseWRebin->Clone(Form("fhResponseWNormPerBin_pthb%d", ipthb+1));
   
   // number of bins after the rebinning for the loops
   Int_t nBinsPPar = hNJetsPerPtPBinThisPtHBin->GetXaxis()->GetNbins(); //note this!!!
   Int_t nBinsMPar = fhResponseN[ipthb]->GetAxis(fAxisMp)->GetNbins();
   Int_t nBinsPDet = fhResponseN[ipthb]->GetAxis(fAxisPtd)->GetNbins();
   Int_t nBinsMDet = fhResponseN[ipthb]->GetAxis(fAxisMd)->GetNbins();
   
   fhResponseN[ipthb]->Reset();
   //Printf("Dimensions %d, proj %d, %d", fhResponseWRebin->GetNdimensions(), fAxisMp, fAxisPtp);
   Printf("Starting a loop over %d Pt bins, %d Mass bins at particle level\nPer each of these bins a loop over %d Pt bins and %d M bins at detector level will be run. Might take a while", nBinsPPar, nBinsMPar, nBinsPDet, nBinsMDet);
   
   //particle level projection
   fhptMparW[ipthb] = fhResponseWRebin->Projection(fAxisMp, fAxisPtp);
   fhptMparW[ipthb]->Sumw2();
   fhptMparW[ipthb]->SetName(Form("hptMpar_pthb%d", ipthb));
   
   TH1D *hPtparW = fhResponseWRebin->Projection(fAxisPtp);
   hPtparW->SetName(Form("hPtparW_%d", ipthb));
   
   //TCanvas *cResDet = new TCanvas("cResDet", "Detector response in each part bin", 1000, 1000);
   //cResDet->Divide(20, 20);
   for(Int_t ibppar = 0; ibppar < nBinsPPar; ibppar++){ //loops over the bins of input hNJetsPerPtPBinThisPtHBin
      
      Double_t xvalue = hNJetsPerPtPBinThisPtHBin->GetBinCenter(ibppar+1);
      Int_t binRespW = hPtparW->FindBin(xvalue);
      //print for debugging
      //Printf("This is pTP bin %d (input) and %d (resp)", ibppar+1, binRespW);
      
      Double_t factor = 1;
      Double_t njetsPtPsel = hPtparW->GetBinContent(binRespW);
      Double_t njetsPtPtot = hNJetsPerPtPBinThisPtHBin->GetBinContent(ibppar+1);
      
      factor = njetsPtPtot/njetsPtPsel;
      
      //loop over the other bins of the response matrix
      for(Int_t ibmpar = 0; ibmpar < nBinsMPar; ibmpar++){ //M par
      	 //loop on pTdet, Mdet
      	 for(Int_t ibpdet = 0; ibpdet < nBinsPDet; ibpdet++){ // Pt det
      	    
      	    for(Int_t ibmdet = 0; ibmdet < nBinsMDet; ibmdet++){ // M det
      	       
      	       Int_t thisbin[4] = {ibmdet+1, ibmpar+1, ibpdet+1, binRespW}; //same order as the original response
      	       
      	       //use the content of each bin to fill the normalized ThnSparse per each pTdet, Mdet and pTpar,Mpar
      	       Double_t contN = fhResponseWRebin->GetBinContent(thisbin); //x, y
      	       Double_t erroN = fhResponseWRebin->GetBinError(thisbin);
      	       
      	       if(contN < 1e-13) {
      	       	  fhResponseN[ipthb]->SetBinContent(thisbin, 0);
      	       	  continue;
      	       }
      	       
      	       //print for debugging
      	       //Printf("It was %f, fill with %f", contN, factor);
      	       fhResponseN[ipthb]->SetBinContent(thisbin, contN * factor);
      	       fhResponseN[ipthb]->SetBinError(thisbin,   factor * erroN);
      	       
      	       
      	    }
      	 }
      }
   }
   
}

//-------------------------------------------------------------------------------------------------
Bool_t WeightClass::WeightResponse2D(Int_t ipthb){
   
   /// weight response for this pT hard bin
   if(!fh2dResponse[ipthb]){
      Printf("No 2D response to be weighted");
      return 0;
   }
   
   // projections useful for QA
   //if(fAxisPtp == 0) fhptpar[ipthb]   = (TH1*) fh2dResponse[ipthb]->ProjectionX(Form("fhptMpar_pTh%d", ipthb+1));
   //if(fAxisPtp == 1) fhptpar[ipthb]   = (TH1*) fh2dResponse[ipthb]->ProjectionY(Form("fhptMpar_pTh%d", ipthb+1));
   //Double_t integralptpar = fhptpar[ipthb]->Integral();
   
   // definition of the scale according to the input from ReadWeight
   Double_t scale = fxsec/ftrials;
   fhScaleF->Fill(ipthb+1, scale);
   
   fh2dResponseW[ipthb] = (TH2*)fh2dResponse[ipthb]->Clone(Form("fh2dResponseW_pTh%d", ipthb+1));
   fh2dResponseW[ipthb]->Sumw2();
   //fh2dResponseW[ipthb]->Scale(1./integralptMpar);
   fh2dResponseW[ipthb]->Scale(scale);
   
   //if (fDoNormPerInt) NormalizePerIntegral(ipthb);
   //else {
   //   fh2dResponseN[ipthb] = 0;
   //   fhptMparW[ipthb] = fh2dResponseW[ipthb]->Projection(fAxisMd, fAxisPtd); // x = fAxisMd, y = fAxisPtd
   //   fhptMparW[ipthb]->SetName(Form("fhptMparW_pTh%d", ipthb+1));
   //   fhptMdetW[ipthb] = fh2dResponseW[ipthb]->Projection(fAxisMd, fAxisPtd); // x = fAxisMd, y = fAxisPtd
   //   fhptMdetW[ipthb]->SetName(Form("fhptMdetW_pTh%d", ipthb+1));
   //}

   return 1;
}

//-------------------------------------------------------------------------------------------------
Bool_t WeightClass::Result(){
   
   /// sum each weighted pT hard bin up
   Printf("Summing response in pT hard bins %d-%d", fFirstpthB, fNpthB);
   
   for(Int_t ipthb = fFirstpthB; ipthb < fNpthB; ipthb++){
      if(fIs2D) {
      	 if(!fh2dResponseW[ipthb]){
      	    Printf("Weighted response %d not found, return", ipthb);
      	    return 0;
      	 }
      	 if(!fh2DRespFinal) {
      	    fh2DRespFinal = (TH2*)fh2dResponseW[ipthb]->Clone("fh2DResponseFinal");
      	    fh2DRespFinal->Sumw2();
      	 } else {
      	    fh2DRespFinal->Add(fh2dResponseW[ipthb]);
      	 }
      } else {
      	 if(!fhResponseW[ipthb]) {
      	    Printf("Weighted response %d not found, return", ipthb);
      	    return 0;
      	 } else {
      	    if(fDoNormPerInt){ // in this case I have the response in fhResponseN
      	       if(!fhResponseN[ipthb]) {
      	       	  Printf("Weighted response %d not found, return", ipthb);
      	       	  return 0;
      	       }
      	    } else { // fhResponseN will be used in the sum, so I put the content of fhResponseW in it
      	       fhResponseN[ipthb] = fhResponseW[ipthb];
      	       Printf("not normalised, Normalized %p, original %p", fhResponseN[ipthb], fhResponseW[ipthb]);
      	    }
      	 }
      	 if(!fhResponseWFinal) {
      	 	 Printf("Pointer %p", fhResponseN[ipthb]);
      	    fhResponseWFinal = (THnSparseF*)fhResponseN[ipthb]->Clone("fhResponseFinal");
      	    fhResponseWFinal->Sumw2();
      	 } else {
      	    fhResponseWFinal->Add(fhResponseN[ipthb]);
      	 }
      	 
      	 if(fAxisMdEmb != -1 && fAxisPtdEmb != -1){
      	    if(!fhMPt3DW[ipthb]) {
      	       Printf("Weighted MPt3D %d not found, return", ipthb);
      	       return 0;
      	    }
      	    if(!fhMPt3DWFinal) {
      	       fhMPt3DWFinal = (THnSparseF*)fhMPt3DW[ipthb]->Clone("fhMPt3DFinal");
      	       fhMPt3DWFinal->Sumw2();
      	    } else {
      	       fhMPt3DWFinal->Add(fhMPt3DW[ipthb]);
      	    }
      	 }
      }
      
   }
   
   if(fIs2D) Printf("The total integral is %f",fh2DRespFinal->Integral());
   else Printf("The total integral is %f",fhResponseWFinal->ComputeIntegral());
   return 1;
}

//-------------------------------------------------------------------------------------------------
TH1D** WeightClass::All1DProjections(THnSparseF *hnsp, Int_t color, TString nameh) const{
   
   const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed-9, kGreen-10, kBlue - 8};
   if(!hnsp){
      Printf("Input is null, cannot project");
      return 0;
   }
   hnsp->Sumw2();
   const Int_t ndim = hnsp->GetNdimensions();
   TH1D **hproj = new TH1D*[ndim];
   for(Int_t iax = 0; iax < ndim; iax++){
      hproj[iax] = hnsp->Projection(iax);
      if(nameh.IsNull()) hproj[iax]->SetName(Form("%s_%d", hnsp->GetName(), iax));
      else hproj[iax]->SetName(Form("%s_%d", nameh.Data(), iax));
      Printf("Name %s", hproj[iax]->GetName());
      hproj[iax]->SetTitle(Form("%s_%d", nameh.Data(), iax));
      hproj[iax]->SetLineColor(colors[color]);
      hproj[iax]->SetLineWidth(2);
      hproj[iax]->SetMarkerColor(colors[color]);
   }

   return hproj;
}

//-------------------------------------------------------------------------------------------------

TCanvas* WeightClass::CompareResponse(THnSparseF *response2) const{
   
   /// compare the 1D projections of the final response and response2. It requires to run the weighting procedure first
   
   return CompareResponses(fhResponseWFinal, response2);
}

//-------------------------------------------------------------------------------------------------
TCanvas* WeightClass::CompareResponses(THnSparseF *response1, THnSparseF *response2, Int_t color1, Int_t color2, TString leg1, TString leg2) const{
   
   ///Compare the 1D projections of response1 and response2. Does not require any of the methods  of the class
      
   if(!response1) {
      Printf("Reponse 1 not found");
      return 0x0;  
   }
   if(!response2) {
      Printf("Reponse 2 not found");
      return 0x0;  
   }
   
   Printf("Draw comparison %s vs %s ", response1->GetName(), response2->GetName());
   
   TH1D** hProjWFin = 0;
   hProjWFin = All1DProjections(response1, color1, "");
   if(!response2) {
      Printf("Input is null, exit");
      return 0;
   }
   
   TH1D** hProj2 = All1DProjections(response2, color2, "");
   //TH1D** hProj2 = All1DProjections(response2, 0, response2->GetName());
   
   TCanvas *cCompare = new TCanvas(Form("cCompare%s%s", response1->GetName(), response2->GetName()), Form("Comp resps %s and %s", response1->GetName(), response2->GetName()), 1000, 1000);
   cCompare->Divide(2,2);
   TCanvas *cRatios = new TCanvas(Form("cRatios%s%s", response1->GetName(), response2->GetName()), Form("Ratios resps %s and %s", response1->GetName(), response2->GetName()), 1000, 1000);
   cRatios->Divide(2,2);
   
   TLegend *leg = new TLegend(0.1, 0.6, 0.4, 0.8);
   leg->SetFillStyle(0);
   
   for(Int_t ipj = 0; ipj<4; ipj++){
      Printf("Axis %d: %s -> %s,  %s -> %s", ipj, response1->GetName(), response1->GetAxis(ipj)->GetTitle(),  response2->GetName(), response2->GetAxis(ipj)->GetTitle());
      
      TH1* h1 = 0x0;
      TH1* h2 = 0x0;

      UniformTH1FForDivide(hProjWFin[ipj], hProj2[ipj], h1, h2, "TH1D");

      cCompare->cd(ipj+1);
      gPad->SetLogy();
      h1->Draw();
      h2->Draw("sames");
      if(ipj == 0) {
      	 leg->AddEntry(hProjWFin[ipj], leg1, "L");
      	 leg->AddEntry(hProj2[ipj], leg2, "L");
      	 leg->Draw();
      }
      TH1* hRatio = (TH1*)h1->Clone(TString::Format("hRatio%s%s", h1->GetName(), h2->GetName()));
      hRatio->SetTitle(TString::Format("; %s; %s/%s", h1->GetXaxis()->GetTitle(), leg1.Data(), leg2.Data()));
      hRatio->Divide(h2);
      cRatios->cd(ipj+1);
      hRatio->Draw();
   }
   return cCompare;
}

//-------------------------------------------------------------------------------------------------
TCanvas* WeightClass::DrawHistogramTH1(TH1* histogram, Bool_t logy) const {
   
   /// draw one of the output histograms
   Printf("Drawing %s , integral %f", histogram->GetName(), histogram->Integral());
   TCanvas *c = new TCanvas(Form("c%s", histogram->GetName()), Form("%s", histogram->GetName()), 500, 500);
   c->cd();
   if(logy) gPad->SetLogy();
   histogram->DrawClone();
   return c;
}

//-------------------------------------------------------------------------------------------------

TH1* WeightClass::GetProj(TH2* h2, Int_t ax, TString name, Int_t color) const {
   Printf("ax = %d", ax);
   const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed-9, kGreen-10, kBlue - 8};
   
   if(!h2){
      Printf("Intput h2 is null return");
      return 0;
   }
   
   TString namehisto = "";
   TString axisS = "X";
   if(ax == 1) axisS = "Y";
   
   if(name.IsNull()) namehisto = Form("%s_pj%d%s", h2->GetName(), ax, axisS.Data());
   else namehisto = Form("%s_%d%s", name.Data(), ax, axisS.Data());
   
   TH1* hproj = 0x0;
   if(ax == 0) hproj = h2->ProjectionX(namehisto, 0, -1, "E");
   if(ax == 1) hproj = h2->ProjectionY(namehisto, 0, -1, "E");
   
   hproj->SetLineColor(colors[color]);
   hproj->SetLineWidth(2);
   return hproj;
}

//-------------------------------------------------------------------------------------------------

void WeightClass::SaveNumberOfJetPerPtHBin(TString filename, Int_t level) const{
   /// save the number of jets per pT hard bin before weighting. This number has to be calculated in the MC jobs and then used in the embedding jobs for weighting with normalization
   /// level indicated if the uncorrected or the weighted histograms is used.  0 = uncorr, 1 = weighted
   
   TFile *fout = new TFile(filename, "recreate");
   TH1D *hNJetsPerPtPBin[fNpthB];
   
    for(Int_t ipthb = fFirstpthB; ipthb < fNpthB; ipthb++){
          
      Double_t integralptpar = 0;
      hNJetsPerPtPBin[ipthb] = IntegralParPt(ipthb, level, integralptpar);
      if(integralptpar < 0) {
      	 Printf("Warning, integral smaller than 0, skipping");
      	 return;
      }
      TParameter<double> value(Form("IntegralPtH%d", ipthb+1), integralptpar);
      fout->cd();
      value.Write();
      hNJetsPerPtPBin[ipthb]->Write();
   }
   return;
}

//-------------------------------------------------------------------------------------------------

void WeightClass::ReadNumberOfJetPerPtHBin(TString filename) {
   /// read the number of jets per pT hard bin before weighting and save it into a TParameter list. This numbers are used in the embedding jobs for weighting with normalization
   
   TFile *fin = new TFile(filename);
   
   if(!fin->IsOpen()){
      Printf("File %s not found, cannot proceed", filename.Data());
      return;
   }
   
   fListNJets = new TList();
   
   for(Int_t ipthb = fFirstpthB; ipthb < fNpthB; ipthb++){
      
      TParameter<double> *param = (TParameter<double>*)fin->Get(Form("IntegralPtH%d", ipthb+1));
      fListNJets->Add(param);
      TH1D *hNJet = (TH1D*)fin->Get(Form("hNJetsPerPtPBin_PtH%d", ipthb));
      fListNJets->Add(hNJet);
      
   }
   return;
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

//to be removed when library available
TList* WeightClass::ReadFile(TString strIn, TString strLst){
   
   TFile *f = new TFile(strIn.Data());
   if(!f->IsOpen()){
      Printf("File %s not found", strIn.Data());
      return 0x0;
   }
   TList *lst = static_cast<TList*>(f->Get(strLst.Data()));
   if(!lst){
      Printf("Error, list %s not found", strLst.Data());
      f->ls();
      return 0x0;
   }
   return lst;
}

//-------------------------------------------------------------------------------------------------

void WeightClass::UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type) const{
   if(!h1 || !h2){
      Printf("Need valid input!");
      return;
   }
   
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h2->GetNbinsX();
   Double_t binW1 = h1->GetBinWidth(3), binW2 = h2->GetBinWidth(3);
   Double_t low1 = h1->GetBinLowEdge(1), low2 = h2->GetBinLowEdge(1);
   Double_t up1 = h1->GetBinLowEdge(nBins1+1), up2 = h2->GetBinLowEdge(nBins2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBins2, binW2, low2, up2);
   if((nBins1 == nBins2) && ((binW1-binW2) < 1e-6) &&  ((low1-low2) < 1e-6) && ((up1-up2) < 1e-6)) {
      Printf("Histograms can be divided");
      h1bis = (TH1*)h1->Clone(Form("%sb", h1->GetName()));
      h2bis = (TH1*)h2->Clone(Form("%sb", h2->GetName()));
      return;
   }
   Double_t commonlow = TMath::Min(low1, low2), commonup = TMath::Max(up1, up2);
   Double_t commonW = TMath::Max(binW1, binW2);
   Int_t commonNBins = (commonup - commonlow)/commonW;
   Printf("Common) N bins = %d, Width = %f, Range = %f - %f", commonNBins, commonW, commonlow, commonup);
   Double_t minWidth = TMath::Min(binW1, binW2);
   
   Int_t rebin = (Int_t)(commonW/minWidth);
   Printf("%f/%f = %f or %d", commonW, minWidth, commonW/minWidth, rebin);
   
   if(minWidth == binW1) h1->Rebin(rebin);
   else h2->Rebin(rebin);
   
   Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsX(), h2->GetName(), h2->GetNbinsX());
   if(type == "TH1F"){
      h1bis = new TH1F(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      h2bis = new TH1F(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
   } else {
      if(type == "TH1D"){
      	 h1bis = new TH1D(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      	 h2bis = new TH1D(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      } else {
      	 Printf("%s Not defined!", type.Data());
      	 return;
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
      if(fillingbin>0 && fillingbin <= nBins1) {
      	 //Printf("Here, debug %d", nBins1);
      	 fillWith = (Double_t)h1->GetBinContent(fillingbin);
      	 error = (Double_t)h1->GetBinError(fillingbin);
      }
      h1bis ->SetBinContent(i+1, fillWith);
      h1bis ->SetBinError(i+1, error);
      //Printf("Filling %s x = %.3f new bin %d with old bin %d content %e", h1->GetName(), bincentre, i+1, fillingbin, fillWith);
      
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
      //Printf("Filling %s x = %.3f new bin %d with old bin %d content %e", h2->GetName(), bincentre, i+1, fillingbin, fillWith);
   }
   
   
}
