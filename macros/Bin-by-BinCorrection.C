#include <TF1.h>
#include <TList.h>
#include <THnSparse.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TClass.h>
#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>

//global
//binning
const Int_t nptbins = 4;
Double_t ptlims[nptbins+1] = {40., 60., 80., 100., 120};
//const Int_t nptbins = 1;
//Double_t ptlims[nptbins+1] = {5., 120.};
Double_t maxRangeMass[nptbins] = {16., 20., 20., 30.};   
THnSparse* ReadMassResponse(TString filename = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root", TString listname = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged", TString hname = "fhnMassResponse");

TH3F* ReadData(TString infile = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/AnalysisResultsMB.root", TString listname = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv", TString histname = "fh3PtJet1VsMassVsLeadPtTagged_0");
TH1D* ReadCorrection(TFile *file, Int_t bin);
TH1D* ReadParM(TFile *file, Int_t bin);
TH1D* ReadDetM(TFile *file, Int_t bin);
TH1D* ReadFromBinbyBinCorrection(TFile *file, Int_t bin, TString histname = "hFacPoverD_pt");
TH2D* Read2DFromBinbyBinCorrection(TFile *file, TString histname);

void CompareResults(const Int_t ninputs, TString files[], TString hnamebase[], TString legs[], Int_t offset[], Bool_t changeColor = kFALSE, Bool_t writeRatios = kFALSE, Bool_t noUniform = kFALSE);

//Int_t FitFactors(TFile *binbybincorrfile);

//_________________________________________________________________________________________________
//method implementations
void CompareResults(const Int_t ninputs, TString files[], TString hnamebase[], TString legs[], Int_t offset[], Bool_t changeColor, Bool_t writeRatios, Bool_t noUniform){
   // the pt bin number will be appended to the hnamebase
   
   TH1 *histogram[ninputs][nptbins];
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cMass = new TCanvas(Form("cMass"), Form("Mass"), dx, dy);
   cMass->Divide(nx, ny);
   
   TCanvas *cMasR = new TCanvas("cMasR", "Mass Ratios", dx, dy);
   cMasR->Divide(nx, ny);
   
   TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.8);
   leg->SetFillStyle(0);
   
   

   for(Int_t ifile = 0; ifile<ninputs; ifile++){
      
      //need to do it only once per file
      TString hnamesuff = "", hnamepref = "";
      Bool_t isdifferentnamingscheme = kFALSE;
      if(hnamebase[ifile].Contains("Iter")){
      	 isdifferentnamingscheme = kTRUE;
      	 hnamesuff = hnamebase[ifile](6,11);
      	 hnamepref = hnamebase[ifile](0, 6);
      	 Printf("Pref %s Suff %s -> %d", hnamepref.Data(), hnamesuff.Data(), isdifferentnamingscheme);
      }
      cMass->SetName(Form("%s%s", cMass->GetName(), legs[ifile].Data()));
      cMasR->SetName(Form("%s%s", cMasR->GetName(), legs[ifile].Data()));
      //read input file
      TFile *fin = new TFile(files[ifile]);
      if(!fin->IsOpen()) {
      	 Printf("%s not found, skip", files[ifile].Data());
      	 continue;
      }
      Printf("Reading file %s", files[ifile].Data());
      
      //read the histograms from the files
      //for the moment one histogram type at a time
      for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	 histogram[ifile][ipt] = 0x0;
      	 
      	 TString hnamethisloop = "";
      	 if(isdifferentnamingscheme) {
      	    Printf("isdifferentnamingscheme");
      	    hnamethisloop = Form("%s%d%s", hnamepref.Data(),  ipt + offset[ifile], hnamesuff.Data());
      	    
      	     // for the comparison of the unfolded new
      	    
      	 }
      	 else hnamethisloop = Form("%s%d", hnamebase[ifile].Data() ,  ipt + offset[ifile]);
      	 
      	 histogram[ifile][ipt] = (TH1*)fin->Get(hnamethisloop);
      	 if(!histogram[ifile][ipt]) {
      	    Printf("Histogram %s not found, continue", hnamethisloop.Data());
      	    fin->ls();
      	    continue;
      	 } Printf("Read hisrogram %s", hnamethisloop.Data());
      	 
      }
   }

   TFile *fout = 0x0;
   
   if(writeRatios) {
      fout = new TFile(Form("Ratio%sOver%s.root", legs[0].Data(), legs[1].Data()), "recreate");
   }
   
   for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      
      cMass->cd(ipt+1);
      TH1* h1 = 0x0;
      TH1* h2 = 0x0;
      TH1* hR = 0x0;
      
      Int_t color = 0;
      if(changeColor) {
      	 
      	 histogram[0][ipt]->SetLineColor(colors[color]);
      	 histogram[0][ipt]->SetMarkerColor(colors[color]);
      	 
      	 histogram[1][ipt]->SetLineColor(colors[color+1]);
      	 histogram[1][ipt]->SetMarkerColor(colors[color+1]);
      }
      
      if(!noUniform) UniformTH1FForDivide(histogram[0][ipt], histogram[1][ipt], h1, h2, "TH1D", noUniform);
      else {
      	 Printf("The histogram were not manipulated");
      	 h1 = histogram[0][ipt];
      	 h2 = histogram[1][ipt];
      }
      
      if(ipt == 1) {
      	 leg->AddEntry(histogram[0][ipt], legs[0], "LP");
      	 leg->AddEntry(histogram[1][ipt], legs[1], "LP");
      }

      if(h1) {
      	 if(!noUniform) h1->Scale(1./h1->Integral());
      	 h1->GetXaxis()->SetRangeUser(-5, maxRangeMass[ipt]);
      	 h1->GetYaxis()->SetRangeUser(0., 1.5*h1->GetMaximum());
      	 h1->Draw();
      	 hR = (TH1*)h1->Clone(Form("hR%d", ipt));
      	 hR->GetYaxis()->SetRangeUser(-0., 2.);
      }
      if(h2){
      	 if(!noUniform) h2->Scale(1./h2->Integral());
      	 h2->GetXaxis()->SetRangeUser(-5, maxRangeMass[ipt]);
      	 h2->Draw("sames");
      	 if(hR && h2) hR->Divide(h2);
      }
      leg->Draw();

      if(!noUniform){
      	 cMasR->cd(ipt+1);
      	 hR->Draw();
      }
      if(writeRatios) {
      	 h1->Write();
      	 h2->Write();
      	 if(!noUniform) hR->Write();
      }
   }
   if(writeRatios) Printf("The output was saved in a ROOT file");
    
   SaveCv(cMass);
   SaveCv(cMasR);
}
//_________________________________________________________________________________________________
THnSparse* ReadMassResponse(TString filename, TString listname, TString hname){
   
   TList *list = ReadFile(filename, listname);
   if(!list) return 0x0;
   THnSparseF *hspmassresp = dynamic_cast<THnSparseF*>(list->FindObject(hname));
   return hspmassresp;
}

//_________________________________________________________________________________________________

TH3F* ReadData(TString infile, TString listname, TString histname){
   TList *list = ReadFile(infile, listname);
   if(!list) return 0x0;
   TH3F *h3 = dynamic_cast<TH3F*>(list->FindObject(histname));
   return h3;
}

//_________________________________________________________________________________________________

TH1D* ProjectMassPtBin(TH3F* h3, Float_t ptMin, Float_t ptMax){
   
   if(!h3) {
      Printf("3D histogram null");
      return 0x0;
   }
   
   Int_t bins[2] = {h3->GetXaxis()->FindBin(ptMin), h3->GetXaxis()->FindBin(ptMax) -1 };
   
   TH1D* mass = h3->ProjectionY(Form("hMass_pt%.0f-%.0f", ptMin, ptMax), bins[0], bins[1], 0, -1);
   return mass;

}

//_________________________________________________________________________________________________

TH1D* ReadCorrection(TFile *file, Int_t bin){
   
  
   TString histname = "hFacPoverD_pt";
   return ReadFromBinbyBinCorrection(file, bin, histname);
   
}

//_________________________________________________________________________________________________

TH1D* ReadParM(TFile *file, Int_t bin){
   
  
   TString histname = "hMPar_pt";
   return ReadFromBinbyBinCorrection(file, bin, histname);
   
}

//_________________________________________________________________________________________________

TH1D* ReadDetM(TFile *file, Int_t bin){
   
  
   TString histname = "hMDet_pt";
   return ReadFromBinbyBinCorrection(file, bin, histname);
   
}

//______________________________________________________________________________________________________

TH1D* ReadFromBinbyBinCorrection(TFile *file, Int_t bin, TString histname){
   
   if(!file->IsOpen()){
      Printf("No file");
      return 0x0;
   }
   TH1D *hcorr = 0x0;
   if(bin<0) hcorr = dynamic_cast<TH1D*>(file->Get(Form("%s", histname.Data())));
   else hcorr = dynamic_cast<TH1D*>(file->Get(Form("%s%d", histname.Data(), bin)));
   if(!hcorr) Printf("Error!! %s (%d) not found", histname.Data(), bin);
   else Printf("ReadFromBinbyBinCorrection:: return %s", hcorr->GetName());
   return hcorr;
}

//______________________________________________________________________________________________________
TH2D* Read2DFromBinbyBinCorrection(TFile *file, TString histname){
   
   if(!file->IsOpen()){
      Printf("No file");
      return 0x0;
   }
   TH2D *hcorr = 0x0;
   hcorr = dynamic_cast<TH2D*>(file->Get(Form("%s", histname.Data())));
   if(!hcorr) Printf("Error!! %s not found", histname.Data());
   else Printf("ReadFromBinbyBinCorrection:: return %s", hcorr->GetName());
   return hcorr;
}


//______________________________________________________________________________________________________

//main
void WriteMCCorrectionFactors(Int_t type, TString suff = "", Double_t minPtDetCut = -1, Int_t binAx = 0, Bool_t allowneg = kFALSE){
   
   // type 1 = output of task MassResponse
   // type 2 = output of task JetShape
   // type 3 = output of Weighted response
   // the pT cut simulates what happens in the embedding
   // binAx = axis for the pT binning
   //         0 -> what is needed for correction, namely par level for par mass, det or det+fluc for mass det/det plus fluc
   //         1 -> pTpar -> both det/det plus fluc and par are binned in ptpar
   //         2 -> pTdet -> cut at detector level, mass at particle level
   
   // implement different pT cut possibilities: on particle level, on detector level, on det plus fluc level
   
   suff += Form("%d", type);
   if(minPtDetCut>-1) suff += Form("%.0f", minPtDetCut);
   suff += Form("BinT%d", binAx);
   Int_t typecp = type;
   TString inputs[3][3];
   TString ptstring = "det/par";
   if(binAx == 1) ptstring = "par";
   if(binAx == 2) ptstring = "det";
   
   // important! definition of the axes!
   Int_t axisptDet = 2; // pT det level
   Int_t axisptPar = 3; // pT part level
   //Int_t axisptDFl = -1; // pT det plus fluc level   
   Int_t axisPjDet = 0; // M det level
   Int_t axisPjPar = 1; // M part level
   //Int_t axisPjDFl = -1; // M det plus fluc level

   Float_t minmass = 0, minpt = 0, minmasswider = 0;
   if(allowneg) {
      minmass = -5;
      minmasswider = -20;
      minpt = -10;
   }
   //type 1
   inputs[0][0] = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root";
   inputs[1][0] = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged";
   inputs[2][0] = "fhnMassResponse";
   //type 2
   inputs[0][1] = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1055/output/AnalysisResults.root";
   inputs[1][1] = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracksEmb_pT0150_E_scheme_TC";
   inputs[2][1] = "fhnMassResponse_0";
   
   if(type == 11){
      inputs[0][0] = "AnalysisResults.root";
      type = 1;
      
   }
   if(type == 12){
      inputs[0][0] = "AnalysisResultsWeighted.root";
      type = 1;
      
   }
   if(type == 111){
      inputs[0][0] = "AnalysisResults.root";
      inputs[1][0] = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme";
      type = 1;
      
   }
   if(type == 21){
      inputs[0][1] = "AnalysisResults.root";
      type = 2;
      
   }
   if(type == 211){
      inputs[0][1] = "AnalysisResults.root";
      inputs[1][1] = "JetShapeDeriv_JetEmb_AKTChargedR040_PicoTracksEmb_pT0150_E_scheme_TC";
      type = 2;
      suff+="Deriv";
      
   }
   if(type == 22){
      inputs[0][1] = "AnalysisResults.root";
      inputs[1][1] = "JetShapeConst_JetRhosub_AKTChargedR040_PicoTracksTree_pT0150_E_scheme_TC";
      type = 2;
      
   }
   if(type == 2){
      axisptDet = 2; // pT det level
      axisptPar = 3; // pT part level
      axisPjDet = 0; // M det level
      axisPjPar = 1; // M part level
   }
   
   if(type == 3){
      //this is what I have after weighting (coming from JetShape)
      inputs[0][2] = "ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root";
      inputs[1][2] = "";
      inputs[2][2] = "fhResponseFinal";
      
   }
   if(type == 31){
      //this is what I have after weighting (coming from JetMassResponse)
      inputs[0][2] = "ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root";
      inputs[1][2] = "";
      inputs[2][2] = "fhResponseFinal";
      type = 3;
   }
   if(type == 311){
      //this is what I have after weighting (coming from JetMassResponse)
      inputs[0][2] = "ResponseWJetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme.root";
      inputs[1][2] = "";
      inputs[2][2] = "hresponsefinalRhoW";
      type = 3;
   }
   if(type == 32){
      //this is what I have after weighting (coming from Tree)
      inputs[0][2] = "ResponseW.root";
      inputs[1][2] = "";
      inputs[2][2] = "fhResponseFinal";
      type = 3;
   }
   if(type == 33){
      //this is what I have after weighting (coming from JetShape)
      inputs[0][2] = "ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root";
      inputs[1][2] = "";
      inputs[2][2] = "fhResponseFinal";
       type = 3;
   }
   if(type == 331){
      //this is what I have after weighting and rho scaling with method RunWeightForRhoBins of WeightResponse.C (original task JetShape)
      inputs[0][2] = "ResposeJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracksRhoW.root";
      inputs[1][2] = "";
      inputs[2][2] = "hresponsefinalRhoW";
      type = 3;
   }
   if(type == 3321){
      //this is what I have after weighting and rho scaling with method ReAnalyseRespRhoBins of WeightResponse.C (original task JetShape)
      inputs[0][2] = "ResposeJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracksRhoWOut.root";
      inputs[1][2] = "";
      inputs[2][2] = "hRespFacFit1";
      type = 3;
   }
   
   if(type == 3322){
      //this is what I have after weighting and rho scaling with method ReAnalyseRespRhoBins of WeightResponse.C (original task JetShape)
      inputs[0][2] = "ResposeJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracksRhoWOut.root";
      inputs[1][2] = "";
      inputs[2][2] = "hRespFacCalc1";
      type = 3;
   }
   
   if(type == 3) {
      axisptDet = 2; // pT det level (or det plus fluc)
      axisptPar = 3; // pT part level
      //axisptDFl = 5; // pT det plus fluc level 
      axisPjDet = 0; // M det level (or det plus fluc)
      axisPjPar = 1; // M part level
      //axisPjDFl = 4; // M det plus fluc level

   }
   
   Printf("Reading\n %s\n %s\n %s", inputs[0][type-1].Data(), inputs[1][type-1].Data(), inputs[2][type-1].Data());
   
   THnSparse* hspresp = 0x0;
   if(type<2) hspresp = ReadMassResponse(inputs[0][type-1], inputs[1][type-1], inputs[2][type-1]);
   else {
      TFile *fin = new TFile(inputs[0][2]);
      if(!fin->IsOpen()){
      	 Printf("File %s not found", inputs[0][2].Data());
      	 return;
      }
      hspresp = (THnSparseF*)fin->Get(inputs[2][2]);
      fin->ls();
   }
   if(!hspresp){
      Printf("Error, THnSparse not found");
      return;
   }
      
   Int_t axisSelMPar = axisptPar;
   Int_t axisSelMDet = axisptDet;
   Int_t axisPjMPar = axisPjPar;
   Int_t axisPjMDet = axisPjDet;

   if(binAx == 0){
      axisSelMPar = axisptPar;
      axisSelMDet = axisptDet;
   }
   if(binAx == 1){
      axisSelMPar = axisptPar;
      axisSelMDet = axisptPar;
   }
   
   if(binAx == 2){
      axisSelMPar = axisptDet;
      axisSelMDet = axisptDet;
      //axisPjMPar  = axisPjPar;
      //axisPjMDet = axisPjPar;
   }
   
   Printf("########## Check axes: \n pTPar = %d, pTDet = %d, MPar = %d, MDet = %d", axisSelMPar, axisSelMDet, axisPjPar, axisPjDet);
   const Int_t dimresp = hspresp->GetNdimensions();
   Double_t widths[dimresp];
   for(Int_t iax = 0; iax < hspresp->GetNdimensions(); iax++){
      
      widths[iax] = hspresp->GetAxis(iax)->GetBinWidth(1);
      Printf("Ax %d (%s) binW = %f", iax, hspresp->GetAxis(iax)->GetTitle(), widths[iax]);
   }
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   
   TCanvas *cSpectra = new TCanvas("cSpectra", "Mass distributions", dx, dy);
   cSpectra->Divide(nx, ny);
   
   TCanvas *cFactors = new TCanvas("cFactors", "Correction Factors", dx, dy);
   cFactors->Divide(nx, ny);
   
   TCanvas *cpTDistr = new TCanvas("cpTDistr", "Pt distribution", 800, 800);
   
   TCanvas *cMPtPar = new TCanvas("cMPtPar", "Pt vs Mass Particle level", 920, 600);
   TCanvas *cMPtDet = new TCanvas("cMPtDet", "Pt vs Mass Detector level", 920, 600);
   
   TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.8);
   leg->SetFillStyle(0);
      
   TFile *fout = new TFile(Form("BinbyBinCorrection_%dBinT%d.root", typecp, binAx), "recreate");
   
   //apply a pT selection on pTDet 40-100
   Double_t wpT = 1, wM = 1;
   Int_t rebinA[dimresp];
   
   for(Int_t iax = 0; iax < hspresp->GetNdimensions(); iax++) rebinA[iax] = (Int_t)(wpT/widths[iax]);
   
   //the if is temporary, see how it behaves
   if(type != 3) hspresp = hspresp->Rebin(rebinA);
   
   Printf("Projecting %d vs %d and %d vs %d", axisPjMPar, axisSelMPar, axisPjMDet, axisSelMDet);
   //2D distributions
   TH2D* hMPtPar = hspresp->Projection(axisPjMPar, axisSelMPar);
   TH2D* hMPtDet = hspresp->Projection(axisPjMDet, axisSelMDet);
   cMPtPar->cd();
   gPad->SetLogz();
   //hMPtPar->GetZaxis()->SetRangeUser(1e-11, 1e-2);
   hMPtPar->GetZaxis()->SetRangeUser(1e-11, 1e-4);
   hMPtPar->GetXaxis()->SetRangeUser(minpt, 150);
   hMPtPar->GetYaxis()->SetRangeUser(minmasswider, 50);
   hMPtPar->Draw("colz");
   cMPtDet->cd();
   gPad->SetLogz();
   //hMPtDet->GetZaxis()->SetRangeUser(1e-11, 1e-2);
   hMPtDet->GetZaxis()->SetRangeUser(1e-11, 1e-4);
   hMPtDet->GetXaxis()->SetRangeUser(minpt, 150);
   hMPtDet->GetYaxis()->SetRangeUser(minmass, 50);
   hMPtDet->Draw("colz");
   
   Int_t binptAll[2] = {hspresp->GetAxis(axisSelMDet)->FindBin(ptlims[0]), hspresp->GetAxis(axisSelMDet)->FindBin(ptlims[nptbins])-1 };
   hspresp->GetAxis(axisSelMDet)->SetRange(binptAll[0], binptAll[1]);
   
   //determination of the bin corresponding to a minimum pt cut, whe requested
   Int_t minpTDetCutB = 0;
   if(minPtDetCut> -1) minpTDetCutB = hspresp->GetAxis(axisSelMDet)->FindBin(minPtDetCut);
   
   //here no minimum pT cut is applied
   TH1D *hpTD = hspresp->Projection(axisptDet);
   hpTD->SetName(Form("hptD"));
   hpTD->GetXaxis()->SetTitle("#it{p}_{T}");
   hpTD->SetLineColor(kRed);
   hpTD->SetLineWidth(2);
   TH1D *hpTP = hspresp->Projection(axisptPar);
   hpTP->SetName(Form("hptP"));
   hpTP->GetXaxis()->SetTitle("#it{p}_{T}");
   hpTP->SetLineColor(kBlue);
   hpTP->SetLineWidth(2);
   
   hspresp->GetAxis(axisSelMDet)->SetRange(0, -1);

   for(Int_t ipt = 0; ipt<nptbins; ipt++){
      //Use 2D distributions
      //to use the THnSparse methos comment FROM here
      
      //detector level
      Int_t binpt[2] = {hMPtDet->GetXaxis()->FindBin(ptlims[ipt]), hMPtDet->GetXaxis()->FindBin(ptlims[ipt+1])-1 };
      
      hMPtDet->GetXaxis()->SetRange( binpt[0], binpt[1] );
      TH1D *hMDet = hMPtDet->ProjectionY(Form("hMDet_pt%d", ipt), binpt[0], binpt[1]);
      hMDet->SetLineColor(kRed);
      hMDet->SetLineWidth(2);
      hMDet->Sumw2();
      hMDet->Scale(1./hMDet->Integral());
      
      //particle level
      Int_t binptP[2] = {hMPtPar->GetXaxis()->FindBin(ptlims[ipt]), hMPtPar->GetXaxis()->FindBin(ptlims[ipt+1])-1 };
      
      hMPtPar->GetXaxis()->SetRange(binptP[0], binptP[1]);
      
      TH1D *hMPar = hMPtPar->ProjectionY(Form("hMPar_pt%d", ipt), binptP[0], binptP[1]);
      hMPar->SetLineColor(kBlue);
      hMPar->SetLineWidth(2);
      hMPar->Sumw2();
      hMPar->Scale(1./hMPar->Integral());
       
      //pT distributions (note: in mode 1 and 2 they'll be the same plot)
      //det level
      TH1D *hPtDet = hMPtDet->ProjectionX(Form("hPtDet_pt%d", ipt),  0, -1);
      hPtDet->SetLineColor(kRed);
      hPtDet->SetLineWidth(2);
      hPtDet->Sumw2();
      
      //part level
      TH1D *hPtPar = hMPtPar->ProjectionX(Form("hPtPar_pt%d", ipt), 0, -1);
      hPtPar->SetLineColor(kBlue);
      hPtPar->SetLineWidth(2);
      hPtPar->Sumw2();
      
      //to use the THnSparse methos comment TILL here
      
      //use THnSparse
      //to use the TH2D methos comment FROM here
      /*
      //detector level
      Int_t binpt[2] = {hspresp->GetAxis(axisRange)->FindBin(ptlims[ipt]), hspresp->GetAxis(axisRange)->FindBin(ptlims[ipt+1])-1 };
      if(minpTDetCutB < binpt[0]) hspresp->GetAxis(axisRange)->SetRange(binpt[0], binpt[1]);
      else if(minpTDetCutB < binpt[1]) hspresp->GetAxis(axisRange)->SetRange(minpTDetCutB, binpt[1]);
      else {
      	 Printf("Minimum pT cut is larger than this pT bin, skip");
      	 continue;
      }
 
      TH1D *hMDet = hspresp->Projection(axisPjDet);
      hMDet->SetName(Form("hMDet_pt%d", ipt));
      hMDet->SetTitle(Form("Mass %.1f < #it{p}_{T, det} < %.1f GeV/c; #it{M} (GeV)", ptlims[ipt], ptlims[ipt+1]));
      hMDet->SetLineColor(kBlue);
      hMDet->SetLineWidth(2);
      hMDet->Sumw2();
      hMDet->Scale(1./hMDet->Integral());
      
      //to select at pT detector level for both mass at particle and detectol level comment FORM here (to be used only for exercise, not for the real analysis)
      
      //reset the cut in pT dedector level
      hspresp->GetAxis(axisRange)->SetRange(minpTDetCutB, -1);
      
      //particle level
      Int_t binptP[2] = {hspresp->GetAxis(axisptPar)->FindBin(ptlims[ipt]), hspresp->GetAxis(axisptPar)->FindBin(ptlims[ipt+1])-1 };
      hspresp->GetAxis(axisptPar)->SetRange(binptP[0], binptP[1]);
      
      //to select at pT detector level for both mass at particle and detectol level comment TILL here
      
      TH1D *hMPar = hspresp->Projection(axisPjPar);
      hMPar->SetName(Form("hMPar_pt%d", ipt));
      hMPar->SetTitle(Form("Mass %.1f < #it{p}_{T, det} < %.1f GeV/c; #it{M} (GeV)", ptlims[ipt], ptlims[ipt+1]));
      hMPar->SetLineColor(kRed);
      hMPar->SetLineWidth(2);
      hMPar->Sumw2();
      hMPar->Scale(1./hMPar->Integral());
      */
      //to use the TH2D methos comment TILL here (to be used only for exercise, not for the real analysis)
      
      hMDet->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      hMPar->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      
      cSpectra->cd(ipt+1);
      hMDet->Draw("E");
      hMPar->Draw("Esames");
      hspresp->GetAxis(axisSelMPar)->SetRange(0, -1);
      
      TH1D *hFacPoverD = dynamic_cast<TH1D*>(hMPar->Clone(Form("hFacPoverD_pt%d", ipt)));
      hFacPoverD->SetTitle(Form("Mass Particle/Detector %.1f < #it{p}_{T, %s} < %.1f GeV/c; #it{M} (GeV)", ptlims[ipt], ptstring.Data(), ptlims[ipt+1]));
      hFacPoverD->SetLineColor(kBlack);
      hFacPoverD->SetLineWidth(2);
      hFacPoverD->Divide(hMDet);
      
      hFacPoverD->GetXaxis()->SetRangeUser(minmass,  maxRangeMass[ipt]);
      
      cFactors->cd(ipt+1);
      gPad->SetLogy();
      hFacPoverD->Draw("E");
      
      if(ipt == 0){
      	 leg->AddEntry(hMDet, "Detector", "L");
      	 leg->AddEntry(hMPar, "Particle", "L");
      }
      cSpectra->cd(ipt+1);
      leg->Draw();
      
      fout->cd();
      hMDet->Write();
      hMPar->Write();
      hFacPoverD->Write();
      hPtDet->Write();
      hPtPar->Write();
      
      hMPtDet->GetXaxis()->SetRange(0, -1);
      hMPtPar->GetXaxis()->SetRange(0, -1);
   }
   
   fout->cd();
   hpTD->Write();
   hpTP->Write();
   hMPtPar->Write();
   hMPtDet->Write();
   
   cpTDistr->cd();
   gPad->SetLogy();
   hpTD->Draw();
   hpTP->Draw("sames");
   leg->Draw();
   
   SaveCv(cSpectra, suff);
   SaveCv(cFactors, suff);
   SaveCv(cpTDistr, suff);
   SaveCv(cMPtPar, suff);
   SaveCv(cMPtDet, suff);
   
}

//_________________________________________________________________________________________________

void CompareFactors(Int_t mode = 0, TString fac1name = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det/BinbyBinCorrection_1.root", TString fac2name = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/BinbyBinCorrection_2.root", TString leg1 = "Detector", TString leg2 = "Det+Back", Bool_t allowneg = kFALSE){
   // mode 0 uses the BinbyBinCorrection_*.root files
   // mode 1 uses the BinbyBinCorrection_*.root for the detector level and the FullResponseOutput.root for the det+fluc
   
   TFile *file1 = new TFile(fac1name);
   if(!file1){
      Printf("File %s not found",fac1name.Data());
      return;
   }
   TFile *file2 = new TFile(fac2name);
   
   if(!file2){
      Printf("File %s not found",fac2name.Data());
      return;
   }
   
   Float_t minmass = 0, minpt = 0;
   if(allowneg) {
      minmass = -5;
      minpt = -10;
   }
   
   TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.8);
   leg->SetFillStyle(0);
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   
   TCanvas *cCorr = new TCanvas("cCorr", Form("Correction factors (%s - %s)", leg1.Data(), leg2.Data()), dx, dy);
   cCorr->Divide(nx, ny);
   TCanvas *cMassP = new TCanvas("cMassP", Form("Mass particle level (%s - %s)", leg1.Data(), leg2.Data()), dx, dy);
   cMassP->Divide(nx, ny);
   TCanvas *cMassD = new TCanvas("cMassD", Form("Mass detector level (%s - %s)", leg1.Data(), leg2.Data()), dx, dy);
   cMassD->Divide(nx, ny);
   TCanvas *cPtPBins = new TCanvas("cPtPBins", Form("Pt particle level (%s - %s)", leg1.Data(), leg2.Data()), dx, dy);
   cPtPBins->Divide(nx, ny);
   TCanvas *cPtDBins = new TCanvas("cPtDBins", Form("Pt detector level (%s - %s)", leg1.Data(), leg2.Data()), dx, dy);
   cPtDBins->Divide(nx, ny);
   
   TCanvas *cpTP = new TCanvas("cpTP", Form("Pt particle level (%s - %s)", leg1.Data(), leg2.Data()), 800, 800);
   
   TCanvas *cpTD = new TCanvas("cpTD", Form("Pt detector level (%s - %s)", leg1.Data(), leg2.Data()), 800, 800);
   
   TCanvas *cRCor = new TCanvas("cRCor", Form("Correction factor Ratio (%s/%s)", leg1.Data(), leg2.Data()), dx, dy);
   cRCor->Divide(nx, ny);
   
   TCanvas *cRecalCorr = new TCanvas("cRecalCorr", Form("Correction factors recalc (%s/%s)", leg1.Data(), leg2.Data()), dx, dy);
   cRecalCorr->Divide(nx, ny);
   TCanvas *cRecalRCor = new TCanvas("cRecalRCor", Form("Correction factor recalc Ratio (%s/%s)", leg1.Data(), leg2.Data()), dx, dy);
   cRecalRCor->Divide(nx, ny);
   
   TString hname = "", nameP = "", nameD = "";
   
   //TString hname2DPar = "";
   //TH2D *hMPtPar = Read2DFromBinbyBinCorrection(file1);
   
   for(Int_t i = 0; i< nptbins; i++){
      
      TH1D* hFac1T= 0x0;
      TH1D* hFac2T = 0x0;
      
      hFac1T = ReadCorrection(file1, i);
      if(!hFac1T){
      	 Printf("Factor 1 pT %d not found", i);
      	 continue;
      }
      
      if(mode == 0){
      	 hFac2T = ReadCorrection(file2, i);
      	 if(!hFac2T){
      	    Printf("Factor 2 pT %d not found", i);
      	    continue;
      	 }
      }
      
      if(mode == 1){
      	 hname = Form("hFracPoverD_pt%d", i);
      	 nameP = Form("hPtD%d_1", i);
      	 
      	 hFac2T = (TH1D*)file2->Get(hname);
      }
      
      TH1* hFac1 = 0x0;
      TH1* hFac2 = 0x0;
      UniformTH1FForDivide(hFac1T, hFac2T, hFac1, hFac2, "TH1D");
      hFac1->SetLineColor(colors[0]);
      hFac2->SetLineColor(colors[1]);
      //hFac1->Scale(1./hFac1->Integral());
      //hFac2->Scale(1./hFac2->Integral());
      if(i == 0){
      	 leg->AddEntry(hFac1, Form("%s", leg1.Data()), "L");
	 leg->AddEntry(hFac2, Form("%s", leg2.Data()), "L");
      }
      hFac1->GetXaxis()->SetRangeUser(minmass, maxRangeMass[i]);
      hFac2->GetXaxis()->SetRangeUser(minmass, maxRangeMass[i]);
      cCorr->cd(i+1);
      gPad->SetLogy();
      hFac1->Draw();
      hFac2->Draw("sames");
      leg->Draw();
      
      TH1* hFacR = (TH1*)hFac1->Clone(Form("hFac1_pt%d", i));
      hFacR->SetTitle(Form("Ratio; %s; %s/%s ", hFac1->GetXaxis()->GetTitle(), leg1.Data(), leg2.Data()));
      hFacR->Divide(hFac2);
      cRCor->cd(i+1);
      hFacR->Draw();
      
      TH1D* hMP1T = 0x0;
      TH1D* hMP2T = 0x0;
      
      hname = "hMPar_pt";
      hMP1T = ReadFromBinbyBinCorrection(file1, i, hname);
      if(!hMP1T){
      	 Printf("Mass particle level %d not found", i);
      	 continue;
      }
      if(mode == 0){
      	 hMP2T = ReadFromBinbyBinCorrection(file2, i, hname);
      	 if(!hMP2T){
      	    Printf("Mass particle level %d not found", i);
      	    continue;
      	 }
      }
      if(mode == 1){
      	 nameP = Form("hpTD%d_1", i);      	 
      	 hMP2T = (TH1D*)file2->Get(nameP);
      	 if(!hMP2T){
      	    Printf("Mass particle level %s not found", nameP.Data());
      	    continue;
      	 }
      }
      
      TH1* hMP1 = 0x0;
      TH1* hMP2 = 0x0;
      UniformTH1FForDivide(hMP1T, hMP2T, hMP1, hMP2, "TH1D");
      hMP1->SetLineColor(colors[0]);
      hMP2->SetLineColor(colors[1]);
      hMP1T->SetLineColor(colors[0]);
      hMP2T->SetLineColor(colors[1]);
      hMP1->GetXaxis()->SetRangeUser(minmass, maxRangeMass[i]);
      hMP2->GetXaxis()->SetRangeUser(minmass, maxRangeMass[i]);
      cMassP->cd(i+1);
      hMP1->Scale(1./hMP1->Integral());
      hMP2->Scale(1./hMP2->Integral());
      hMP1->Draw();
      hMP2->Draw("sames");
      leg->Draw();
      
      TH1D* hMD1T = 0x0;
      TH1D* hMD2T = 0x0;
      
      hname = "hMDet_pt";
      hMD1T = ReadFromBinbyBinCorrection(file1, i, hname);
      if(!hMD1T){
      	 Printf("Mass particle level %d not found", i);
      	 continue;
      }
      if(mode == 0){
      	 hMD2T = ReadFromBinbyBinCorrection(file2, i, hname);
      	 if(!hMD2T){
      	    Printf("Mass particle level %d not found", i);
      	    continue;
      	 }
      }
      if(mode == 1){
      	 nameD = Form("hpTD%d_0", i);      	 
      	 hMD2T = (TH1D*)file2->Get(nameD);
      }
      
      TH1* hMD1 = 0x0;
      TH1* hMD2 = 0x0;
      UniformTH1FForDivide(hMD1T, hMD2T, hMD1, hMD2, "TH1D");
      hMD1->SetLineColor(colors[0]);
      hMD2->SetLineColor(colors[1]);
      //hMP1T->SetLineColor(colors[0]);
      //hMP2T->SetLineColor(colors[1]);
      hMD1->GetXaxis()->SetRangeUser(minmass, maxRangeMass[i]);
      hMD2->GetXaxis()->SetRangeUser(minmass, maxRangeMass[i]);
      cMassD->cd(i+1);
      hMD1->Scale(1./hMD1->Integral());
      hMD2->Scale(1./hMD2->Integral());
      hMD1->Draw();
      hMD2->Draw("sames");
      leg->Draw();
      
      //recalculate factors with the rebinned mass spectra
      TH1* hRecFac1 = (TH1*)hMP1->Clone("hRecFac1");
      TH1* hRecFac2 = (TH1*)hMP2->Clone("hRecFac2");
      hRecFac1->Divide(hMD1);
      hRecFac2->Divide(hMD2);
      
      cRecalCorr->cd(i+1);
      gPad->SetLogy();
      hRecFac1->Draw();
      hRecFac2->Draw("sames");
      leg->Draw();
      
      TH1* hRecFacR = (TH1*)hRecFac1->Clone(Form("hRecFac1_pt%d", i));
      hRecFacR->SetTitle(Form("Recalculated Ratio; %s; %s/%s ", hRecFac1->GetXaxis()->GetTitle(), leg1.Data(), leg2.Data()));
      hRecFacR->Divide(hFac2);
      cRecalRCor->cd(i+1);
      hRecFacR->Draw();
      
      //pT distribution
      TH1D* hpTBinsD1T = 0x0;
      TH1D* hpTBinsP1T = 0x0;
      
      hname = "hPtDet_pt";
      hpTBinsD1T = ReadFromBinbyBinCorrection(file1, i, hname);
      hname = "hPtPar_pt";
      hpTBinsP1T = ReadFromBinbyBinCorrection(file1, i, hname);
      
      TH1D* hpTBinsD2T = 0x0;
      TH1D* hpTBinsP2T = 0x0;

      hname = "hPtDet_pt";
      hpTBinsD2T = ReadFromBinbyBinCorrection(file2, i, hname);
      hname = "hPtPar_pt";                        
      hpTBinsP2T = ReadFromBinbyBinCorrection(file2, i, hname);
      if(hpTBinsD1T && hpTBinsP1T && hpTBinsD2T && hpTBinsP2T){
      	 hpTBinsD1T->SetLineColor(colors[0]);
      	 hpTBinsD2T->SetLineColor(colors[1]);
      	 
      	 hpTBinsP1T->SetLineColor(colors[0]);
      	 hpTBinsP2T->SetLineColor(colors[1]);
      	 
      	 cPtPBins->cd(i+1);
      	 gPad->SetLogy();
      	 hpTBinsP1T->Draw("");
      	 hpTBinsP2T->Draw("sames");
      	 leg->Draw();
      	 
      	 cPtDBins->cd(i+1);
      	 gPad->SetLogy();
      	 hpTBinsD1T->Draw("");
      	 hpTBinsD2T->Draw("sames");
      	 leg->Draw();
      }
   }
   
   TH1D* hpTD1T = 0x0;
   TH1D* hpTD2T = 0x0;
   
   hname = "hptD";
   hpTD1T = ReadFromBinbyBinCorrection(file1, -1, hname);
   if(!hpTD1T){
      Printf("pT detector level %s in file %s not found", hname.Data(), fac1name.Data());
      return;
   }
   if(mode == 0){
      hpTD2T = ReadFromBinbyBinCorrection(file2, -1, hname);
      if(!hpTD2T){
      	 Printf("pT detector level %s in file %s not found", hname.Data(), fac2name.Data());
      	 return;
      }
   }
   if(mode == 1){
      nameD = Form("hRespWProj_2");
      
      hpTD2T = (TH1D*)file2->Get(nameD);
   }
   TH1* hpTD1 = 0x0;
   TH1* hpTD2 = 0x0;
   UniformTH1FForDivide(hpTD1T, hpTD2T, hpTD1, hpTD2, "TH1D");
   hpTD1->SetLineColor(colors[0]);
   hpTD2->SetLineColor(colors[1]);
   hpTD1->GetXaxis()->SetRangeUser(0, 100);
   hpTD2->GetXaxis()->SetRangeUser(0, 100);
   cpTD->cd();
   hpTD1->Scale(1./hpTD1->Integral());
   hpTD2->Scale(1./hpTD2->Integral());
   hpTD1->Draw();
   hpTD2->Draw("sames");
   leg->Draw();
   
   TH1D* hpTP1T = 0x0;
   TH1D* hpTP2T = 0x0;
   
   hname = "hptP";
   hpTP1T = ReadFromBinbyBinCorrection(file1, -1, hname);
   if(!hpTP1T){
      Printf("pT particle level %s in file %s not found", hname.Data(), fac1name.Data());
      return;
   }
   if(mode == 0){
      hpTP2T = ReadFromBinbyBinCorrection(file2, -1, hname);
      if(!hpTP2T){
      	 Printf("pT particle level %s in file %s not found", hname.Data(), fac2name.Data());
      	 return;
      }
   }
   if(mode == 1){
      
      nameP = Form("hRespWProj_3");
      
      hpTP2T = (TH1D*)file2->Get(nameP);
      
   }
   TH1* hpTP1 = 0x0;
   TH1* hpTP2 = 0x0;
   UniformTH1FForDivide(hpTP1T, hpTP2T, hpTP1, hpTP2, "TH1D");
   hpTP1->SetLineColor(colors[0]);
   hpTP2->SetLineColor(colors[1]);
   hpTP1->GetXaxis()->SetRangeUser(0, 100);
   hpTP2->GetXaxis()->SetRangeUser(0, 100);
   cpTP->cd();
   //temporary
   //Int_t bincut = hpTP1->GetXaxis()->FindBin(20.);
   //Int_t binmax = hpTP1->GetXaxis()->FindBin(200.);
   //Double_t norm = hpTP1->Integral(bincut, binmax);
   hpTP1->Scale(1./hpTP1->Integral());
   hpTP2->Scale(1./hpTP2->Integral());
   hpTP1->Draw();
   hpTP2->Draw("sames");
   leg->Draw();
   
   SaveCv(cCorr);
   SaveCv(cMassP);
   SaveCv(cMassD);
   SaveCv(cpTP);
   SaveCv(cpTD);
   SaveCv(cRCor);
   SaveCv(cRecalCorr);
   SaveCv(cRecalRCor);
}



//_________________________________________________________________________________________________

void CorrectBinbyBin(Int_t mode = 3, TString filebinbybinname = "BinbyBinCorrection_1.root", TString infileMBEJE80 = "/data/Work/jets/JetMass/pPbJetMassAnalysis/EJE/output806-807/resultMBbelow80EJEabove/MassOutput.root", TString parametrization = "",  Int_t bkgSub = 0, Int_t rbmass = 1, Bool_t allowneg = kFALSE, Int_t colorCorr = 3){

   /// mode = 1 derivative method
   /// mode = 2 constituent method
   /// mode = 3 read directly the sum of the MB and EJE output, background sub type to be set by hand, see below -> __ recommended__
   /// mode = 4 closure test, "data" is detector level distribution
   ///
   /// parametrization : if (!parametrization.IsNull()) use the parametrization of the response found there
   ///  bkgSub: when mode == 3 need to choose the background subtraction method to be selected : = 0; // const sub, = 1;       // no sub,  2;       // deriv sub
   
   
   //read data
   TH3F *h3DataMB = 0x0;
   TH3F *h3DataJT = 0x0;
   // to be set MANUALLY!!!!
   Int_t firstPtBinData = 1;
   if(mode == 4) firstPtBinData = 0;
   
   Printf("GOING TO START FROM BIN %d IN DATA", firstPtBinData);
   
   Float_t minmass = 0, minpt = 0;
   if(allowneg || bkgSub == 2) {
      minmass = -5;
      minpt = -10;
   }

   Int_t nhist = 1;
   if(mode == 4) nhist = nptbins; 
   if(mode == 3) nhist = 5;
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nhist, nx, ny, dx, dy);
   TString bkgSubString = "";
   const Int_t nhistInFile = nhist;
   TH1D *hmassRawRead[nhistInFile];
   
   
   if(mode < 3){
      TString infileMB = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/AnalysisResultsMB.root";
      TString listMB;
      if(mode == 1) listMB = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDeriv";
      if(mode == 2) listMB = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeConstSub_TC";
      h3DataMB = ReadData(infileMB, listMB);
      if(!h3DataMB) {
      	 Printf("File MB not found");
      	 return;
      }
      
      TString infileJT = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/AnalysisResultsEJE.root";
      TString listJT;
      if(mode == 1) {
      	 listJT = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCDerivJ1";
      	 bkgSubString = "Deriv";
      }
      if(mode == 2) {
      	 listJT = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeConstSub_TCConstJ1";
      	 bkgSubString = "Const";
      }
      
      h3DataJT = ReadData(infileJT, listJT);
      if(!h3DataJT) {
      	 Printf("File triggered not found");
      	 return;
      }
   } else {
      TString hnamebase = "";
      TFile *fin = 0x0;
      if(mode == 3){
      	 //here can still specify which background subtraction I want
      	 if(bkgSub == 2) bkgSubString = "Deriv";
      	 if(bkgSub == 0) bkgSubString = "Const";
      	 if(bkgSub == 1) bkgSubString = "NoSub";
      	 
      	 fin = new TFile(infileMBEJE80);
      	 if(!fin->IsOpen()){
      	    Printf("Error! File %s not found", infileMBEJE80.Data());
      	    return;
      	 }
      	 
      	 
      	 hnamebase = "hMassTagged";
      	 hnamebase+=bkgSub;
      	 hnamebase+="pTj";
      	 
      	 
      }
      
      if(mode == 4){
      	 fin = new TFile(filebinbybinname);
      	 hnamebase = "hMDet_pt";
      	 if(!fin->IsOpen()){
      	    Printf("Error! File %s not found", filebinbybinname.Data());
      	    return;
      	 }
      	 
      }
      TCanvas *cMassOrig = new TCanvas(Form("cMassOrig%d%s", mode, bkgSubString.Data()), Form("Mass before manipulation (%s)", bkgSubString.Data()), dx, dy);
      cMassOrig->Divide(nx, ny);
   
      for(Int_t ipt = 0; ipt < nhistInFile; ipt++){
      	 hmassRawRead[ipt] = (TH1D*)fin->Get(Form("%s%d", hnamebase.Data(), ipt));
      	 if(!hmassRawRead[ipt]) {
      	    Printf("Error at bin %d, mass not found", ipt);
      	    continue;
      	 }
      	 if(ipt == nhistInFile -1 ) Printf("Histo %s mean = %f", hmassRawRead[ipt]->GetName(), hmassRawRead[ipt]->GetMean());
      	 cMassOrig->cd(ipt+1);
      	 hmassRawRead[ipt]->Draw();
      }
   }
  
   
   //read correction factors
   TFile *fbinbybinCorr = new TFile(filebinbybinname);
   if(!fbinbybinCorr->IsOpen()){
      Printf("File for corrections %s not found", filebinbybinname.Data());
      return;
   }
   
   //read parametrization, if any
   TFile *fparam = 0x0;
   if (!parametrization.IsNull()){
      fparam = new TFile(parametrization);
      if(!fparam->IsOpen()) {
      	 Printf("File %s not found", parametrization.Data());
      	 return;
      }
   }
   
   
   TLegend *leg = new TLegend(0.5, 0.6, 0.85, 0.9, "Scaled to 1");
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   TLegend *legR=  new TLegend(0.1, 0.7, 0.45, 0.9);
   legR->SetFillStyle(0);
   legR->SetBorderSize(0);
   
   TString rebst = "";
   if(rbmass>1) {
      Printf("REBINNING");
      rebst = "rb";
      rebst += rbmass;
      
   }

   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cMass = new TCanvas(Form("cMassT%d%s%s", mode, bkgSubString.Data(), rebst.Data()), Form("Mass (%s)", bkgSubString.Data()), dx, dy);
   cMass->Divide(nx, ny);
   TCanvas *cMasR = new TCanvas(Form("cMasRT%d%s%s", mode, bkgSubString.Data(), rebst.Data()), Form("Mass Ratios (%s)", bkgSubString.Data()), dx, dy);
   cMasR->Divide(nx, ny);
   
   TCanvas *cMassMC = new TCanvas(Form("cMassMC%d%s%s", mode, bkgSubString.Data(), rebst.Data()), Form("Mass MC (%s)", bkgSubString.Data()), dx, dy);
   cMassMC->Divide(nx, ny);
   
   TCanvas *cCorF = new TCanvas(Form("cCorF%d%s%s", mode, bkgSubString.Data(), rebst.Data()), Form("Correction factors (%s)", bkgSubString.Data()), dx, dy);
   cCorF->Divide(nx, ny);
   
   TCanvas *cMasRelUnc = new TCanvas(Form("cMasRelUnc%d%s%s", mode, bkgSubString.Data(), rebst.Data()), Form("Mass Relative uncertainty (%s)", bkgSubString.Data()), dx, dy);
   cMasRelUnc->Divide(nx, ny);
   
   //mass
   TH1 *hmassRaw[nptbins];
   TH1 *hbbbCorf[nptbins];
   TH1 *hmassCor[nptbins];
   TH1 *hmassCorParam[nptbins];
   TH1 *hmassPar[nptbins];
   TH1 *hmassDet[nptbins];
   
   // mass relative uncertainty
   TH1F* hmassRawRelUnc[nptbins];
   TH1F* hmassCorRelUnc[nptbins];
   
   //ratios
   TH1 *hmassUC[nptbins];
   TH1 *hmassCo[nptbins];
   TH1 *hmassDetDiv[nptbins];
   TH1 *hmassParDiv[nptbins];
   
   TPaveText** pvpt = GetPavePtBins(nptbins, ptlims);
   TFile *fout = new TFile(Form("BinbyBinCorrectedMass_%d%s%s.root",mode, bkgSubString.Data(), rebst.Data()), "recreate");
   
   for(Int_t ipt = 0; ipt < nptbins; ipt++){
      
      Printf("\n----------------------------------------------------------------------------------------------------------------");
      Printf("THIS IS ipt = %d. Looking at pT bin %d for data and %d for correction", ipt, ipt+firstPtBinData, ipt);
      
      hmassRaw[ipt] = 0x0;
      hmassPar[ipt] = 0x0;
      hmassDet[ipt] = 0x0;
      TH1D *htmpRaw = 0x0;
      TH1D *htmpPar = 0x0;
      TH1D *htmpDet = 0x0;
      TH1D *htmpCorf = 0x0;
     
      
      if(mode < 3) {
      	 //Data mass distribution
      	 if(ptlims[ipt]<80) { //MB
      	    htmpRaw = ProjectMassPtBin(h3DataMB, ptlims[ipt], ptlims[ipt+1]);
      	 } else {
      	    
      	    htmpRaw = ProjectMassPtBin(h3DataJT, ptlims[ipt], ptlims[ipt+1]);
      	    htmpRaw->Add(ProjectMassPtBin(h3DataMB, ptlims[ipt], ptlims[ipt+1]));
      	 }
      	 if(!htmpRaw) continue;
      } else {
      	 htmpRaw = hmassRawRead[ipt+firstPtBinData];
      }
      
      htmpRaw->Rebin(rbmass);
      
      if(ipt == nptbins -1 ) Printf("Loop on factors %d Histo %s, mean = %f", ipt, htmpRaw->GetName(), htmpRaw->GetMean());
      
      //correction factors
      //htmpCorf = ReadCorrection(fbinbybinCorr, ipt);
      
      //if(!htmpCorf) continue;
      
      TH1 *hbbbMRawbD = 0; //not used, it's the same as hmassRaw[ipt]
      
      //MC particle level mass spectrum
      htmpPar = ReadParM(fbinbybinCorr, ipt);
      UniformTH1FForDivide(htmpPar, htmpRaw, hmassPar[ipt], hmassRaw[ipt], "TH1D");
      
      //MC detector level mass spectrum
      htmpDet = ReadDetM(fbinbybinCorr, ipt);
      UniformTH1FForDivide(htmpDet, htmpRaw, hmassDet[ipt], hbbbMRawbD, "TH1D");
      hmassRaw[ipt]->SetName(Form("%sbinFin", hmassRaw[ipt]->GetName()));
       
      //build the correction factors
      hbbbCorf[ipt] = (TH1D*)hmassPar[ipt]->Clone();
      hbbbCorf[ipt]->SetName(Form("hCorrFactor_pT%d", ipt));
      hbbbCorf[ipt]->SetTitle(Form("Correction Factor bin %d; #it{M} (GeV/#it{c}); Par/Det", ipt));
      
      Bool_t candivide = hbbbCorf[ipt]->Divide(hmassDet[ipt]);
      if(!candivide) {
      	 Printf("Could not divide, factors not available. Exit");
      	 continue;
      }
     
      if(ipt == nptbins -1 ) Printf("After Uniform %d Histo %s, mean = %f", ipt, hmassRaw[ipt]->GetName(), hmassRaw[ipt]->GetMean());

      if(ipt == nptbins -1 ) Printf("Detlevel Before %s mean = %f, After uniform %s mean = %f", htmpDet->GetName(), htmpDet->GetMean(), hmassDet[ipt]->GetName(), hmassDet[ipt]->GetMean());
      hmassRaw[ipt]->SetMarkerStyle(20);
      hmassRaw[ipt]->SetMarkerColor(colors[0]);
      hmassRaw[ipt]->SetLineWidth(2);
      hmassRaw[ipt]->SetLineColor(colors[0]);
      hmassRaw[ipt]->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      Double_t areaRaw = hmassRaw[ipt]->Integral();
      
      hmassCor[ipt] = (TH1D*)hmassRaw[ipt]->Clone(Form("hMassCorr_pt%d", ipt));
      //Printf("COMPARE::: lower bin %.2f vs %.2f; width %.2f vs %.2f", hmassCor[ipt]->GetBinLowEdge(1),  hbbbCorf[ipt]->GetBinLowEdge(1), hmassCor[ipt]->GetBinWidth(1), hbbbCorf[ipt]->GetBinWidth(1));
      Bool_t cancorrect = hmassCor[ipt]->Multiply(hbbbCorf[ipt]);
      if(!cancorrect){
      	 Printf("Cannot apply correction factors! Exit");
      	 continue;
      }
      hmassCor[ipt]->SetMarkerStyle(24);
      hmassCor[ipt]->SetMarkerColor(colors[colorCorr]);
      hmassCor[ipt]->SetLineWidth(2);
      hmassCor[ipt]->SetLineColor(colors[colorCorr]);
      hmassCor[ipt]->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);

      //define histograms to be filled with relative uncertainties
      hmassRawRelUnc[ipt] = (TH1F*)hmassRaw[ipt]->Clone(Form("hmassRawRelUnc_pt%d", ipt));
      hmassRawRelUnc[ipt]->Reset(); //reset bin content and errors
      hmassRawRelUnc[ipt]->SetMarkerStyle(20);
      hmassRawRelUnc[ipt]->SetMarkerColor(hmassRaw[ipt]->GetMarkerColor());
      hmassCorRelUnc[ipt] = (TH1F*)hmassRawRelUnc[ipt]->Clone(Form("hmassCorRelUnc_pt%d", ipt));
      hmassCorRelUnc[ipt]->SetMarkerColor(hmassCor[ipt]->GetMarkerColor());
      hmassCorRelUnc[ipt]->SetMarkerStyle(24);
            
      cCorF->cd(ipt+1);
      hbbbCorf[ipt]->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      hbbbCorf[ipt]->Draw();
      //htmpCorf->Draw("sames");
      if(fparam){
      	 TString funcname = Form("fpol5_pT%d", ipt);
      	 TF1 *fCorFParam = (TF1*)fparam->Get(funcname);
      	 hmassCorParam[ipt] = (TH1D*)hmassRaw[ipt]->Clone(Form("hmassCorParam_pT%d", ipt));
      	 hmassCorParam[ipt]->SetMarkerColor(kMagenta);
      	 hmassCorParam[ipt]->SetMarkerStyle(24);
      	 hmassCorParam[ipt]->Multiply(fCorFParam);
      	 
      }
      //ratios
      Int_t thisiszero =  UniformTH1FForDivide(hmassRaw[ipt], hmassDet[ipt], hmassUC[ipt], hmassDetDiv[ipt], "TH1D"); //should simply do the clones of hmassRaw[ipt] and hmassDet[ipt]
      if(thisiszero != 0) Printf("Strange -- det level!!");
      hmassUC[ipt]->SetName(Form("hmassUC%d", ipt));
      hmassUC[ipt]->SetTitle(Form("Ratio Data/MC (%.0f < #it{p}_{T} < %.0f); #it{M}; Raw/Det", ptlims[ipt], ptlims[ipt+1]));
      hmassUC[ipt]->SetMarkerStyle(20);
      hmassUC[ipt]->SetMarkerColor(colors[0]);
      hmassUC[ipt]->SetLineWidth(2);
      hmassUC[ipt]->SetLineColor(colors[0]);
      hmassUC[ipt]->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      
      thisiszero = UniformTH1FForDivide(hmassCor[ipt], hmassPar[ipt], hmassCo[ipt], hmassParDiv[ipt], "TH1D"); //should simply do the clones of hmassCor[ipt] and hmassPar[ipt]
      if(thisiszero != 0) Printf("Strange -- par level!!");
      //hmassParDiv[ipt]= (TH1*)hmassPar[ipt]->Clone(Form("hmassParDiv%d", ipt));
      
      hmassCo[ipt] =  (TH1*)hmassCor[ipt]->Clone();
      hmassCo[ipt]->SetName(Form("hmassCo%d", ipt));
      hmassCo[ipt]->SetTitle(Form("Ratio Data/MC (%.0f < #it{p}_{T} < %.0f); #it{M}; Corr/Par", ptlims[ipt], ptlims[ipt+1]));
      hmassCo[ipt]->SetMarkerStyle(24);
      hmassCo[ipt]->SetMarkerColor(colors[3]);
      hmassCo[ipt]->SetLineWidth(2);
      hmassCo[ipt]->SetLineColor(colors[3]);
      hmassCo[ipt]->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      
      hmassUC[ipt]->Divide(hmassDetDiv[ipt]);
      hmassCo[ipt]->Divide(hmassParDiv[ipt]);
      //hmassParDiv[ipt]->Divide(hmassCo[ipt]);  //debugging
      if(ipt == 0){
      	 legR->AddEntry(hmassUC[ipt], "Uncorrected", "PL");
      	 legR->AddEntry(hmassCo[ipt], "Corrected", "PL");
      }
      
      	 
      hmassPar[ipt]->SetMarkerStyle(30);
      hmassPar[ipt]->SetMarkerColor(colors[2]);
      hmassPar[ipt]->SetLineColor(colors[2]);
      //if(mode != 4) hmassPar[ipt]->Scale(areaRaw/hmassPar[ipt]->Integral());
      
      hmassDet[ipt]->SetMarkerStyle(30);
      hmassDet[ipt]->SetMarkerColor(colors[1]);
      hmassDet[ipt]->SetLineColor(colors[1]);
      //if(mode != 4) hmassDet[ipt]->Scale(areaRaw/hmassDet[ipt]->Integral());
      
      if(ipt == 0){
      	 leg->AddEntry(hmassCor[ipt], "Corrected", "LP");
	 leg->AddEntry(hmassRaw[ipt], "Raw", "LP");
	 leg->AddEntry(hmassPar[ipt], "(MC) Particle", "LP");
	 leg->AddEntry(hmassDet[ipt], "(MC) Detector", "LP");
      }
      cMass->cd(ipt+1);
      hmassCor[ipt]->Draw();
      hmassRaw[ipt]->Draw("sames");
      if(fparam) hmassCorParam[ipt]->Draw("sames");
      
      pvpt[ipt]->DrawClone();
      legR->Draw();

      cMassMC->cd(ipt+1);
      TH1*hmassParCp =  hmassPar[ipt]->DrawNormalized(); //
      hmassParCp->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      hmassParCp->GetYaxis()->SetRangeUser(0, hmassParCp->GetMaximum()*1.5);
      
      hmassDet[ipt]->DrawNormalized("sames");
      hmassRaw[ipt]->DrawNormalized("sames");
      hmassCor[ipt]->DrawNormalized("sames");
      leg->Draw();
      
      cMasR->cd(ipt+1);
      hmassUC[ipt]->Draw();
      //hmassParDiv[ipt]->Draw(); //debugging
      hmassCo[ipt]->Draw("sames");
      //hbbbCorf[ipt]->Draw("sames"); //debugging
      legR->Draw();
      
      fout->cd();
      hmassRaw[ipt]->Write();
      hmassCor[ipt]->Write();
      hmassUC[ipt] ->Write();
      hmassCo[ipt] ->Write();
      hbbbCorf[ipt] ->Write();
      
      // loop on the mass entries to calculate relative uncertainties
      for(Int_t ibin = 0; ibin < hmassRaw[ipt]->GetNbinsX(); ibin++){
      	 if(hmassRaw[ipt]->GetBinContent(ibin+1) > 1e-10){
      	    hmassRawRelUnc[ipt]->SetBinContent(ibin+1, hmassRaw[ipt]->GetBinError(ibin+1)/hmassRaw[ipt]->GetBinContent(ibin+1));
      	    if(hmassCor[ipt]->GetBinContent(ibin+1)  > 1e-10) hmassCorRelUnc[ipt]->SetBinContent(ibin+1, hmassCor[ipt]->GetBinError(ibin+1)/hmassCor[ipt]->GetBinContent(ibin+1));
      	 }
      }
      
      cMasRelUnc->cd(ipt+1);
      hmassRawRelUnc[ipt]->Draw("P");
      hmassCorRelUnc[ipt]->Draw("Psames");
   }
   
   SaveCv(cMass);
   SaveCv(cMassMC);
   SaveCv(cMasR);
   SaveCv(cMasRelUnc);
   SaveCv(cCorF);
   
   
}

//_________________________________________________________________________________________________

void CompareDataMC(TString finMC = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/binbybin/BinbyBinCorrection_1.root", TString findata = "BinbyBinCorrectedMass_2.root", Bool_t allowneg = kFALSE){
   
   TFile *fMC = new TFile(finMC);
   TFile *fData = new TFile(findata);

   Float_t minmass = 0, minpt = 0;
   if(allowneg) {
      minmass = -5;
      minpt = -10;
   }

   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cMass = new TCanvas(Form("cMass"), Form("Mass (%s - %s)", finMC.Data(), findata.Data()), dx, dy);
   cMass->Divide(nx, ny);
   TCanvas *cMasR = new TCanvas(Form("cMasR"), Form("Mass Ratio Data/Part (%s - %s)", finMC.Data(), findata.Data()), dx, dy);
   cMasR->Divide(nx, ny);
   
   TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.8);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   
   for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      TH1 *hPar = 0x0;
      TH1 *hDet = 0x0;
      TH1 *hCor = 0x0;
      TH1 *hRatio = 0x0;
      TH1D *htmpPar = ReadFromBinbyBinCorrection(fMC, ipt, "hMPar_pt");
      TH1D *htmpDet = ReadFromBinbyBinCorrection(fMC, ipt, "hMDet_pt");
      TH1D *htmpCor = ReadFromBinbyBinCorrection(fData, ipt, "hMassCorr_pt");
      cMass->cd(ipt+1);
      if(htmpCor) {
      	 if(ipt == 0) leg->AddEntry(htmpCor, "Corrected", "LP");
      	 htmpCor->Scale(1./htmpCor->Integral());
      	 htmpCor->Draw();
      }
      if(htmpPar) {
      	 UniformTH1FForDivide(htmpPar, htmpCor, hPar, hRatio, "TH1D");
      	 if(ipt == 0) leg->AddEntry(hPar, "Particle", "LP");
      	 hPar->Scale(1./hPar->Integral());
      	 hPar->Draw("sames");
      	 
      	 hRatio->Divide(hPar);
      }
      if(htmpDet) {
      	 UniformTH1FForDivide(htmpDet, htmpCor, hDet, hCor, "TH1D");
      	 if(ipt == 0) leg->AddEntry(hDet, "Detector", "LP");
      	 hDet->Scale(1./hDet->Integral());
      	 hDet->Draw("sames");
      }
      
      leg->Draw();
      
      hRatio->GetXaxis()->SetRangeUser(minmass, maxRangeMass[ipt]);
      cMasR->cd(ipt+1);
      hRatio->Draw();
      
   }
   
   SaveCv(cMass);
   SaveCv(cMasR);
}

//_____________________________________________________________________________

void CompareDataMCRec(TString finMC = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/BinByBinCorrDetPtParPt/BinbyBinCorrection_111BinT0.root", TString findata = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/InputDataEJEMB/dataRerunMacroMBEJETrigComb1/MassOutput.root", Int_t bkgsubtype = 2, Int_t offsetOnData = 1, TString mcspecleg = "", Bool_t changeBinning = kFALSE){ //, TString* leg, Bool_t allowneg = kFALSE
	
	/// bkgsubtype: 0 = Const sub, 1 = no sub, 2 = Deriv Sub
	/// the response has to be chosen with the corresponding subtraction method
   
   TFile *fMC = new TFile(finMC);
   TFile *fData = new TFile(findata);
   
   if(!fMC->IsOpen()){
   	   Printf("File MC not found (%s)", finMC.Data());
   	   return;
   }
   if(!fData->IsOpen()){
   	   Printf("File Data not found (%s)", findata.Data());
   	   return;
   }
   
   TString hnameData = Form("hMassTagged%dpTj", bkgsubtype);
   TString hnameMC   = "hMDet_pt";
   TString bkgSubName = "NoBkg";
   if(bkgsubtype == 0) bkgSubName = "ConstSub";
   if(bkgsubtype == 2) bkgSubName = "DerivSub";
   
   Int_t nptbins = 4;
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   
   TCanvas *cMassDataMCRec = new TCanvas(Form("cMassData%sMCRec%s", bkgSubName.Data(), mcspecleg.Data()), Form("Mass Data%s and MCRec %s", bkgSubName.Data(), mcspecleg.Data()), 900, 800);
   cMassDataMCRec->Divide(nx, ny);
   
   TLegend *leg = new TLegend(0.5, 0.4, 0.9, 0.6);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   
   
   for(Int_t ipt = 0; ipt<nptbins; ipt++){
   
   	   TH1D* hMC = (TH1D*)fMC->Get(Form("%s%d", hnameMC.Data(), ipt));
   	   if(!hMC){
   	   	   Printf("MC bin %s%d not found", hnameMC.Data(), ipt);
   	   	   fMC->ls();
   	   	   continue;
   	   }
   	   hMC->Scale(1./hMC->Integral("width"));
   	   TH1D* hDa = (TH1D*)fData->Get(Form("%s%d", hnameData.Data(), ipt+offsetOnData));
   	   if(!hDa){
   	   	   Printf("MC bin %s%d not found", hnameData.Data(), ipt);
   	   	   fData->ls();
   	   	   continue;
   	   }
   	   hDa->Scale(1./hDa->Integral("width"));
   	   
   	   if(ipt == 0){
   	   	   leg->AddEntry(hMC, Form("MCRec%s", mcspecleg.Data()), "LP");
   	   	   leg->AddEntry(hDa, Form("Data%s", bkgSubName.Data()), "LP");
   	   }
   	   cMassDataMCRec->cd(ipt+1);
   	   hMC->Draw();
   	   hDa->Draw("Psames");
   	   leg->Draw();
   
   }
   
   SaveCv(cMassDataMCRec);
   
   
}
//_________________________________________________________________________________________________

void CompareResultFixedInput(){
   
   const Int_t ninputs = 2;
   
   TString files[ninputs] = //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_2.root", "/data/cernbox/JetMassAnalysis/Data/pPb/150501/TranformedIntoTH1.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det/BinbyBinCorrectedMass_2.root", "/data/cernbox/JetMassAnalysis/Data/pPb/150501/TranformedIntoTH1.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_2.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det/BinbyBinCorrectedMass_2.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPtParPt/Det/BinbyBinCorrectedMass_2.root", "/data/cernbox/JetMassAnalysis/Data/pPb/150501/TranformedIntoTH1.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPtParPt/DetPlusFluc/BinbyBinCorrectedMass_2.root", "/data/cernbox/JetMassAnalysis/Data/pPb/150501/TranformedIntoTH1.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPtParPt/Det576/BinbyBinCorrectedMass_2.root", "/data/cernbox/JetMassAnalysis/Data/pPb/150501/TranformedIntoTH1.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_3Derivrb2.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; //
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/outputpThardBins/analysis/BinByBinCorrDetPtParPt/BinbyBinCorrection_31.root", "/data/cernbox/JetMassAnalysis/Data/pPb/150501/TranformedIntoTH1.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/systematics/ovlExclu/BinbyBinCorrectedMass_3.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160330pth2_PPtC20/BinbyBinCorrectedMass_3.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Const.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160330pth2_PPtC10/BinbyBinCorrectedMass_3.root"};
   //files[0] = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Deriv.root";
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Const.root", "/data/Work/jets/JetMass/PbPbResults/TranformedIntoTH1.root"};
   //{"/data/Work/jets/JetMass/PbPbResults/TranformedIntoTH1.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Derivrb2.root"};
   //{"/data/Work/jets/JetMass/PbPbResults/DistJetMassAllSyst_Area.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/BinByBinCorrDetPtParPt/BinbyBinCorrection_111BinT0.root"};//
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/CmpPbPb/RatioPbPb2DUnfOverBbyB-B+F-DerivSub.root", "/data/Work/jets/JetMass/PbPbResults/CmpPYTHIA/RatioPYTHIA576OverPbPb2DUnf.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/CmpConstDeriv/BinbyBinCorrectedMass_3Constrb2.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/CmpConstDeriv/BinbyBinCorrectedMass_3Derivrb2.root"};
   /// compare bin-by-bin corrected Const and Deriv rebin 2
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/CmpConstDeriv/CmpRb1/BinbyBinCorrectedMass_3Const.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/CmpConstDeriv/CmpRb1/BinbyBinCorrectedMass_3Deriv.root"}; /// compare bin-by-bin corrected Const and Deriv 
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MBEJE/00resp_pt20_ptT10/MassUnfSumMBEJE.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/CmpConstDeriv/CmpRb1/BinbyBinCorrectedMass_3Const.root"}; /// Const sub compare bin-by-bin corrected and unfolded
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MBEJE/00resp_pt20_ptT10/MassUnfSumMBEJE.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/CmpConstDeriv/CmpRb1/BinbyBinCorrectedMass_3Deriv.root"}; /// Deriv sub compare bin-by-bin corrected and unfolded
   
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Derivrb2.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_3Derivrb2.root"}; /// Compare corr det+fluc-derivSub vs corr det; Raw DerivSub rebin 2
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Deriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_3Deriv.root"}; /// Compare corr det+fluc-derivSub vs corr det; Raw DerivSub
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Derivrb2.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; //
   //{"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/BinByBinCorrParPt/CmpPbPbPYTHIA/RatioPYTHIA276TeVOverPYTHIA1134.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/CmpPbPb/RatioPbPb2DUnfOverBbyB-B+F-DerivSub.root"}; /// comparison ratios PbPb/pPb vs Pyhtia PbPb/ pythia pPb
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MB/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; /// compare my unfolding MB with Marta's results pPb
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_3Const.root"}; /// compare unf and bin-by-bin for detector only (using Marta results for unfolding to have all the stat)
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MB/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; /// compare my unfolding det MB with Marta's results pPb
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/MBbelow80EJEabove/UnfoldedDistributionsPrior0MB80EJEabove.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; /// compare my unfolding det MB+EJE with Marta's results pPb
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/MBbelow80EJEabove/UseMatrixMB/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; /// compare my unfolding det MB+EJE using the response calculared in /data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB with Marta's results pPb -> this was a cross-checks, gives the same as the one above
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/SumUnfolded/MassUnfSumMBEJE.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; /// compare my unfolding det MB+EJE (separate and summed) with Marta's results pPb
   
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/SumUnfolded/MassUnfSumMBEJE.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/MBbelow80EJEabove/UnfoldedDistributionsPrior0MB80EJEabove.root"}; /// compare my unfolding det MB+EJE (separate and summed) with MB+EJEscaled raw unfolded
   
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/EJE/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; /// compare my unfolding det+fluc EJE with Marta's results pPb
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/UnfoldedDistributionsPrior0Test.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160419pth2_PPtC10/BinbyBinCorrectedMass_3Deriv.root"}; /// compare my unfolding det+fluc MB with bin-by-bin det+fluc Deriv
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/UnfoldedDistributionsPrior0ConstMB.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160330pth2_PPtC10/BinbyBinCorrectedMass_3Const.root"}; /// compare my unfolding det+fluc MB with bin-by-bin det+fluc Const
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160330pth2_PPtC10/BinbyBinCorrectedMass_3Const.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_3.root"}; /// bin bin by bin corr det fluc const over det only
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160330pth2_PPtC10/CmpWithUnf2D/DetFlucMB/RatioNewUnfDetFlcCnstMBOverBbyB-B+F-ConstSub.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/CmpWithBinbyBin/Ratio2DunfOverBbyB1134.root"}; /// Compare ratio unfolding/bin-by bin for the resp det+ fluc (const) and det only (unfoldig det only is from Marta, all stat)
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160330pth2_PPtC10/BinbyBinCorrectedMass_3Constrb2.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_3Constrb2.root"}; /// comparison bin-by-bin fluc+det constsub vs bin-by-bin det only
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/UnfoldedDistributionsPrior0ConstMB.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"}; ///unfolding det+fluc (const) vs unfolding det only (Marta) 
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/CmpDet/RatioNewUnfDetFlcCnstMBOver2Dunf.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160330/pTH2/pTCut10/analysis/WeightNormPerBin/DetFluc/FineBinning/BinByBinCorrectionDetPtParPt/CmpDetOnly/RatioBbyB-B+F-ConstSubOverBbyB1134.root"}; /// ratio det+fluc (const) over det only, comparison ratio for bin-by-bin and unfolding
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgDeriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgConst.root"};
   // data Deriv sub unfolded with Deriv sub compared to unfolded with Const sub (MB only)
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgConst.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/UnfoldDataConstSubWithDeriv/UnfoldedDistributionsPrior0.root"};
   // data Const sub unfolded with Const sub compared to unfolded with Deriv sub (MB only)   
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgDeriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/OneRhoBinNoWeight/UnfoldedDistributionsPrior0.root"};
   // compare standard unfolding deriv with one rho bin Deriv
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgConst.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/OneRhoBinNoWeight/UnfoldedDistributionsPrior0.root"};
   // compare standard unfolding const with one rho bin Deriv
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/OneRhoBinNoWeight/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root"};
   // compare unfolding Marta pPb  with one rho bin Deriv
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgConst.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/UnfoldedDistributionsPrior0.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/updatedRooUnfold/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgConst.root"};
   //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetFlucNoBkgSub/PtDetPtPar/BinbyBinCorrectedMass_3NoSub.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/Det1134/BinbyBinCorrectedMass_3NoSub.root"};
   {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MBEJE/00resp_pt20_80or70_110_ptT10_140or50_140_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/ConstSub/00resp_pt20_80or70_110_ptT10_140or50_140_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root"};
    
   TString hnamebase[ninputs] = 
   //{"hMassCorr_pt", "hUnfMass_PtBin"};       // this is for the comparison my pPb corrected data vs Marta's unfolding
   //{"hMassCorr_pt", "hMassCorr_pt"};          // this is for the comparison among different versions of my pPb corrected data
   //{"hMPar_pt", "hUnfMass_PtBin"};             // this is for the comparison of my  PYTHIA vs Marta's unfolding
   //{"hMPyt_PtBin", "hMPar_pt"};         // this is for the comparison of Marta's PYTHIA with mine
   //{"hUnfMass_PtBin" , "hMassCorr_pt"};                 // this is for the comparison Marta's corrected data with pPb corrected data
   //{"hR", "hR"};                                              // this is the ratio gotten from this same method
   //{"hMUnf__Iter3", "hUnfMass_PtBin"};          // unfolding partial with Marta results
   //{"hMUnf__Iter3", "hMassCorr_pt"};              // unfolding partial with my bin-by-bin correction
   //{"hUnfMpj_Itr3_ptb", "hUnfMass_PtBin"};          // unfolding sum trig and Marta's results
   //{"hUnfMpj_Itr3_ptb", "hMUnf__Iter3"};          // unfolding sum MB and trig and alltogether
   {"hMUnf__Iter3", "hMUnf__Iter3"};
   
   TString legs[ninputs] = 
   //{"BbyB-B+F", "2Dunf"};
   //{"2Dunf", "BbyB1134"};
   //{"576Par", "2Dunf"};
   //{"1134", "576"};
   //{"BbyB-B+F-noOvl", "BbyB-B+F"};
   //{"BbyB-B+F-ConstSub", "BbyB-B+F-DerivSub"};
   //{"BbyB-B+F-ConstSub", "PbPb2DUnf"};
   //legs[0].Prepend("RawDerivSub");
   //{"PYTHIA276TeV", "PYTHIA1134"};
   //{"PbPb2DUnf", "BbyB-B+F-DerivSub"};
   //{"Pyhia276O502", "PbPbOpPb"};
   //{"BbyB-B+F-DerivSub", "BbyB1134"};
   //{"BbyB-B+F-ConstSub", "BbyB1134"};
   //{"NewUnfMB", "2Dunf"};
   //{"NewUnfMBEJE", "2Dunf"};
   //{"NewUnfMBEJErespMB", "2Dunf"};
   //{"NewUnfMBUnfEJESum", "2Dunf"};
   //{"NewUnfMBUnfEJEStick", "2Dunf"};
   //{"NewUnfMBUnfEJESum", "NewUnfMBEJE"};
   //{"NewUnfDetFlcMB", "2Dunf"};
   //{"NewUnfDetFlcEJE", "2Dunf"};
   //{"NewUnfDetFlcDrivMB", "BbyB-B+F-DerivSub"};
   //{"UnfDetFlcDriv", "BbyB-B+F-Deriv"};
   //{"NewUnfDetFlcCnstMB", "BbyB-B+F-ConstSub"};
   //{"UnfDetFlcConst", "BbyB-B+F-Const"};
   //{"NewUnfDetFlcCnstMB", "2Dunf"};
   //{"UnfOBbyB-FConst", "UnfOBbyBDetOnly"};
   //{"DetFlConstODet-unf", "DetFlConstDetOnly-BbyB-F"};
   //{"UnfDetFlDeriv", "DataDerivResConst"};
   //{"UnfDetFlDeriv", "UnfDetFlConst"};
   //{"UnfDetFlConst", "DataConstResDeriv"};
   //{"UnfDetFlDeriv","UnfDetFlDerivRho15_30"};
   //{"UnfDetFlConst","UnfDetFlDerivRho15_30"};
   //{"UnfDetFlDerivRho15_30", "2Dunf"};
   //{"UnfDetFlConst", "UnfNoBkgSub"};
   //{"UnfOldConst", "UnfNewConts"};
   //{"BbBDetFlNoSub", "BbBDetNoSub"};
   {"UnfDetFlNoSub", "UnfDetNoSub"};
   Int_t firstMatchingBin[ninputs] = 
   //{1, 0}; // 2D unfolding has different numbering of bins (check this again and again)
   {0, 0};
   //{1, 1}; // 2D unfolding and Marta's PYTHIA has different numbering of bins (check this again and again)
   //{0, 1};
   
   Bool_t changeColor = kFALSE;
   if(hnamebase[0] == hnamebase[1]) changeColor = kTRUE;
   //changeColor = kTRUE;
   Bool_t writeout = kTRUE;
   Bool_t nouniform = kFALSE;
   if(hnamebase[0] =="hR") nouniform = kTRUE;
   CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform);
}


void TransformCustom( TString namein = "/data/cernbox/JetMassAnalysis/Data/pPb/150501/DistJetMassAllSyst_Area.root", TString nameout = "TranformedIntoTH1.root"){
   
   TFile *fin = new TFile(namein);
   if(!fin->IsOpen()){
      Printf("File %s not found, exit", namein.Data());
      return;
   }
   TList *list = fin->GetListOfKeys();
   Int_t nkeys = list->GetEntries();
   TFile *fout = new TFile(nameout, "recreate");
   Printf("%d entries in file", nkeys);
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nkeys, nx, ny, dx, dy);
   TCanvas *cMass = new TCanvas(Form("cMass"), Form("Mass"), dx, dy);
   cMass->Divide(nx, ny);
   Int_t cc = 0;
   for(Int_t i = 0; i< nkeys; i++){
      TString nameh = list->At(i)->GetName();
      TClass* objtype=fin->Get(nameh)->IsA();
      TString tpname=objtype->GetName();
      Printf("Read %s %s", tpname.Data(), nameh.Data());      
      if(tpname == "TGraphErrors"){
      	 
      	 TGraphErrors *graph = (TGraphErrors*)fin->Get(nameh);
      	 Printf("(%p)", graph);
      	 TString newname = nameh(2, nameh.Length()-1);
      	 newname.Prepend("h");
      	 //newname.Append(nameh());
      	 Int_t nBins = graph->GetN();
      	 cMass->cd(cc+1);
      	 graph->Draw("AP");
      	 
      	 //TH1F *hrg = graph->GetHistogram();
      	 Double_t minx, maxx, y, erminx, ermaxx, er;
      	 graph->GetPoint(0, minx, y);
      	 graph->GetPoint(nBins-1, maxx, y);
      	 erminx = graph->GetErrorX(0);
      	 ermaxx = graph->GetErrorX(nBins-1);
      	 Printf("Max %f + %f", maxx, ermaxx);
      	 TH1F *h = new TH1F(Form("%s", newname.Data()), Form("%s;%s;%s", graph->GetTitle(),  graph->GetXaxis()->GetTitle(),  graph->GetYaxis()->GetTitle()), nBins, minx - erminx, maxx + ermaxx);
      	 h->SetName(newname);
      	 h->SetMarkerStyle(20);
      	 h->SetMarkerColor(kOrange+7);
      	 h->SetLineColor(kOrange+7);
      	 h->SetLineWidth(2);
      	 Printf("Range %s (%f, %f), Nbins %d", h->GetName(), h->GetBinLowEdge(1), h->GetBinLowEdge(nBins + 1), nBins);
      	 cc++;
      	 for(Int_t ibin = 0; ibin<nBins; ibin++){
      	    Double_t px, py, ex, ey;
      	    graph->GetPoint(ibin, px, py);
      	    ex = graph->GetErrorX(ibin);
      	    ey = graph->GetErrorY(ibin);
      	    h->SetBinContent(ibin+1, py);
      	    h->SetBinError(ibin+1, ey);
      	    Printf("Point %d , %f, %f", ibin, px, py);      	    
      	 }
      	 
      	 fout->cd();
      	 h->Write();
      
      }
   }
}

Int_t RunFitFactors(TString filebinbybinname){
   
   /// fits the factor with a polinomial function to be used for correction
   /// examine the distance between the fit function and the histogram to assing a systematic
   /// ? examine different orders 5-6

   TFile *binbybincorrfile = new TFile(filebinbybinname);
//   FitFactors(fin);
//   
//   return;
//}
//
//Int_t FitFactors(TFile *binbybincorrfile){
   if(!binbybincorrfile->IsOpen()){
      Printf("File not found");
      return -1;
   }
   TString facname = "hPtDet_pt";
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cFacFit = new TCanvas(Form("cFacFit"), Form("Factors fitted"), dx, dy);
   cFacFit->Divide(nx, ny);
  TCanvas *cDeltaDataFit = new TCanvas(Form("cDeltaDataFit"), Form("Differce between the factors  and the fit function"), dx, dy);
  cDeltaDataFit->Divide(nx, ny);
  
  const Int_t nfunctype = 2;
  TF1 ***fpol = new TF1**[nfunctype];
  TString funcname[nfunctype] = {"fpol4", "fpol5"};
  TString funcformula[nfunctype] = { "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x",  "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x"};
   
  TGraph ***gDiff = new TGraph**[nfunctype];
  
  TLegend *leg = new TLegend(0.2, 0.7, 0.45, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
   
   Int_t success = 0;
   TFile *fileout = new TFile(Form("Parametrization%s", filebinbybinname.Data()), "recreate");
   
   for(Int_t ifunc = 0; ifunc<nfunctype; ifunc++){
      
      fpol[ifunc] = new TF1*[nptbins];
      gDiff[ifunc] = new TGraph*[nptbins];
      
      for(Int_t ipt = 0; ipt<nptbins; ipt++){
      	 
      	 fpol[ifunc][ipt] = new TF1(Form("%s_pT%d", funcname[ifunc].Data(), ipt), funcformula[ifunc].Data(), 0, maxRangeMass[ipt]);
      	 fpol[ifunc][ipt]->SetParameters(1., 0.5 , 0.1 , 1.e-3, 1.e-5);
      	 //fpol->SetParameters(0.796462,
      	 //   0.370401,
      	 //   -0.189753,
      	 //   0.0358117,
      	 //   -0.00249496,
      	 //   7.92823e-05);
      	 fpol[ifunc][ipt]->SetLineColor(colors[ifunc+1]);
      	 if(ipt == 0) leg->AddEntry(fpol[ifunc][ipt], funcname[ifunc], "L");
      	 TH1D* hFac = ReadCorrection(binbybincorrfile, ipt);
      	 if(!hFac){
      	    Printf("Histogram %d not found", ipt);
      	    continue;
      	 }
      	 
      	 Int_t nbins = hFac->GetNbinsX();
      	 gDiff[ifunc][ipt] = new TGraph(nbins);
      	 gDiff[ifunc][ipt]->SetMarkerColor(colors[ifunc+1]);
      	 gDiff[ifunc][ipt]->SetMarkerStyle(20);
      	 gDiff[ifunc][ipt]->SetName(Form("gDiff%s_pT%d", funcname[ifunc].Data(), ipt));
      	 gDiff[ifunc][ipt]->SetTitle("Diff (histogram - func)/histogram; #it{M} (GeV/#it{c}); (H - pol)/H");
      	 Int_t fitstatus = hFac->Fit(fpol[ifunc][ipt]->GetName(), "L0R+");
      	 //TFitResultPtr *fitres = hFac->Fit(fpol->GetName(), "L0SR+");
      	 if(fitstatus) continue;
      	 
      	 success++;
      	 for(Int_t ib = 0; ib < nbins; ib++){
      	    Double_t xvalue = hFac->GetBinCenter(ib+1);
      	    Double_t yvalue = hFac->GetBinContent(ib+1);
      	    if(TMath::Abs(yvalue) > 1.e-5){
      	       Double_t diff = (yvalue - fpol[ifunc][ipt]->Eval(xvalue))/yvalue;
      	       gDiff[ifunc][ipt]->SetPoint(ib, xvalue, diff);
      	    }
      	 }
      	 
      	 cFacFit->cd(ipt+1);
      	 if(ifunc == 0) 
      	    hFac->DrawClone();
      	 else leg->Draw();
      	 fpol[ifunc][ipt]->Draw("sames");
      	 
      	 cDeltaDataFit->cd(ipt+1);
      	 gDiff[ifunc][ipt]->GetXaxis()->SetRangeUser(0, maxRangeMass[ipt]);
      	  if(ifunc == 0)  gDiff[ifunc][ipt]->Draw("AP");
      	  else {
      	     gDiff[ifunc][ipt]->Draw("P");
      	     leg->Draw();
      	  }
      	 
      	 fileout->cd();
      	 fpol[ifunc][ipt]->Write();
      	 gDiff[ifunc][ipt]->Write();
      }
     
   }
    SaveCv(cFacFit);
    SaveCv(cDeltaDataFit);
   
   return success;
   
}
