#include <TF1.h>
#include <TProfile.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>
#include </data/macros/LoadALICEFigures.C>
#include <Riostream.h>

#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand



void LoopOverTree(TTree *tree, TString names[], TString suffix, TString cmpEmbFilename = "", Bool_t applyptcut = kFALSE);
void LoopOverTreeWithDoubles(TTree *tree, TString names[], TString suffix);
TH1D* DrawSampleTask(TString strIn =  "AnalysisResultsH.root",  TString strLst = "EmcalJetSampleCheckDistr_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TPC_histos", TString hname = "fHistJetsPtArea_0");
TTree* ReadTree (TString strIn, TString nameTree);
void ScaleFactors(Int_t neventsLHC10b = 2e7, TString filename = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1054/outputpTHardBins/mergeRuns/plainMerge/out.root", TString treename = "fTreeJet");
TH1F* GetHistoEmb(TString hname, TList *list);
TH1D* GetHistoEmb(TString hname, TList *list, Int_t axisPj, Int_t* axisSel = 0x0, Float_t* cutL = 0x0, Float_t* cutH = 0x0);

void DrawTVectorInTree(TString strIn =  "/data/Work/jets/JetMass/DetectorCorrections/TaskInputEmbeddingJetVector/test/OutputTree/AnalysisResults.root",  TString strLst = "PrepareInputForEmbedding", Bool_t btlvec = kFALSE){
   
   TList* list = ReadFile(strIn, strLst);
   if(!list) return;
   
   TString nameArea = "fTreeJet";
   //TString nameCons = "fTreeJetC";
   TTree *treeArea = (TTree*) list->FindObject(nameArea);
   //TTree *treeCons = (TTree*) list->FindObject(nameCons);
   
   //if(!treeArea || !treeCons) {
   //   Printf("Tree %s not found %p, %p", nameCons.Data(), treeCons);
   //   return;
   //}
   if(!treeArea) {
      Printf("Tree %s not found %p", nameArea.Data(), treeArea);
      list->ls();
      return;
   }
   
   TFile *fout = new TFile(Form("Tree%s.root", strLst.Data()), "recreate");
   fout->cd();
   treeArea->Write();
   //treeCons->Write();
   if(btlvec){
      TString branches[2] = {"fJetDet", "fJetPart"};
      LoopOverTree(treeArea, branches, "");
   } else {
      TString branches[8] = {"fPt", "fEta", "fPhi", "fMass", "fPtP", "fEtaP", "fPhiP", "fMassP"};
      LoopOverTreeWithDoubles(treeArea, branches, "");
   }
}

//________________________________________________________________________________________________

void DrawTVectorInTree2(TString strIn = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1054/outputpTHardBins/mergeRuns/plainMerge/Randomized/outR.root", TString nameArea = "fTreeJetR", TString suffixoutput = "", TString comparison = "", Bool_t applyptcut = kFALSE){
   TTree* treeArea = ReadTree (strIn, nameArea);
   /* 
   TFile *fout = new TFile(Form("Tree%s.root", suffixoutput.Data()), "recreate");
   fout->cd();
   treeArea->Write();
   //treeCons->Write();
   */
   TString branches[2] = {"fJetDet.", "fJetPart."};
   LoopOverTree(treeArea, branches, "", comparison, applyptcut);
   
}

//________________________________________________________________________________________________

TTree* ReadTree (TString strIn, TString nameTree){
   TFile *fin = new TFile (strIn);
   if(!fin->IsOpen()){
      Printf("%s not found", strIn.Data());
      return 0x0;
   }
   TTree *treeArea = (TTree*) fin->Get(nameTree);
   if(!treeArea) {
      Printf("Tree %s not found %p", nameTree.Data(), treeArea);
      fin->ls();
      return 0x0;
   }
   
   return treeArea;
}

//________________________________________________________________________________________________

TH1D* DrawSampleTask(TString strIn,  TString strLst, TString hname){
   TList *list = ReadFile(strIn, strLst);
   if(!list) {
      Printf("List %s not found in %s", strLst.Data(), strIn.Data());
      list->ls();
      return 0x0;
   }
   
   
   TH2F *hpTArea = (TH2F*)list->FindObject(hname);
   TH1D *hpT = hpTArea->ProjectionX(); 
   TCanvas *cpT = new TCanvas("cpT", "pT");
   cpT->cd();
   gPad->SetLogy();
   hpT->Draw();

   return hpT;
}


void LoopOverTree(TTree *tree, TString names[], TString suffix, TString cmpEmbFilename, Bool_t applyptcut){

   
   //Warning!! Assume that the tree contains TLorentzVectors
   
   if(!tree){
      Printf("Tree null");
      return;
   
   }
   
   TH1F* hpTDEmb = 0x0;
   TH1F* hpTPEmb = 0x0;
   TH1F* hMDEmb  = 0x0;
   TH1F* hMPEmb  = 0x0;
   TH1F* hdpTEmb = 0x0;
   TH1F* hdMEmb  = 0x0;
   if(!cmpEmbFilename.IsNull()){
      TList *listEmb = ReadFile(cmpEmbFilename, "SingleTrackEmbedding");
      if(listEmb){
      	 hpTDEmb = dynamic_cast<TH1F*>(listEmb->FindObject("fhpTEmb"));
      	 hpTPEmb = dynamic_cast<TH1F*>(listEmb->FindObject("fhpTPart"));
      	 hMDEmb  = dynamic_cast<TH1F*>(listEmb->FindObject("fhMEmb"));
      	 hMPEmb  = dynamic_cast<TH1F*>(listEmb->FindObject("fhMPart"));
      	 hdpTEmb = dynamic_cast<TH1F*>(listEmb->FindObject("fhDeltapT"));
      	 hdMEmb  = dynamic_cast<TH1F*>(listEmb->FindObject("fhDeltaM"));
      	 if(!hpTPEmb) {
      	    THnSparse *hsp = dynamic_cast<THnSparse*>(listEmb->FindObject("fhPartJet"));
      	    if(hsp) {
      	       //hsp->GetAxis()->SetRangeUser();
      	       hdMEmb  = (TH1F*)hsp->Projection(0);
      	       hdpTEmb = (TH1F*)hsp->Projection(1);
      	       hMPEmb  = (TH1F*)hsp->Projection(2);
      	       hpTPEmb = (TH1F*)hsp->Projection(3);
      	       
      	       hdMEmb  ->SetName("hdMEmb");
      	       hdpTEmb ->SetName("hdpTEmb");
      	       hMPEmb  ->SetName("hMPEmb");
      	       hpTPEmb ->SetName("hpTPEmb");
      	    }
      	 }
      }
   }
   
   TH1F *hpTR = new TH1F(Form("hpTR_%s", suffix.Data()), Form("hpTR_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 140, -20, 120);
   hpTR->SetLineWidth(2);
   hpTR->SetLineColor(colors[0]);
   hpTR->Sumw2();
   TH1F *hpTP = new TH1F(Form("hpTP_%s", suffix.Data()), Form("hpTP_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 140, -20, 120);
   hpTP->SetLineColor(colors[1]);
   hpTP->Sumw2();
   
   TH1F *hMR = new TH1F(Form("hMR_%s", suffix.Data()), Form("hMR_%s;#it{M} (GeV)", suffix.Data()), 1400, -20, 100);
   hMR->SetLineWidth(2);
   hMR->SetLineColor(colors[0]);
   hMR->Sumw2();
   TH1F *hMP = new TH1F(Form("hMP_%s", suffix.Data()), Form("hMP_%s;#it{M} (GeV)", suffix.Data()), 140, -20, 100);
   hMP->SetLineWidth(2);
   hMP->SetLineColor(colors[1]);
   hMP->Sumw2();
   
   TH2F *hResponsepT = new TH2F("hResponsepT", "Response pT; #it{p}_{T,gen};#it{p}_{T,rec};", 60, -10, 100, 60, -10, 100);
   TH2F *hResponseM  = new TH2F("hResponseM", "Response M; #it{M}_{gen};#it{M}_{rec};", 60, -10, 20, 60, -10, 20);
   
   TH1F *hdM = new TH1F(Form("hdM_%s", suffix.Data()), Form("hdM_%s;#it{#Delta M} (GeV)", suffix.Data()), 60, -30, 30);
   hdM->Sumw2();
   TH1F *hdpT = new TH1F(Form("hdpT_%s", suffix.Data()), Form("hdpT_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -30, 30);
   hdpT->Sumw2();
   TLegend *legDelta = new TLegend(0.5, 0.6, 0.8, 0.8);
   legDelta->SetFillStyle(0);
   legDelta->AddEntry(hdM, "Tree", "PL");
   legDelta->AddEntry(hdMEmb, "Emb Task", "PL");
   
   //bin-by-bin correction 
   Int_t nptbins = 4;
   Double_t ptlims[nptbins+1] = {40, 60, 80, 100, 120};
   TProfile *hMgenOverMrec[nptbins];
   TH1F     *hMRpTRange[nptbins];
   TH1F     *hMPpTRange[nptbins];
   for(Int_t ipt = 0; ipt<nptbins; ipt++){
      hMgenOverMrec[ipt] = new TProfile(Form("hMgenOverMrec_%d", ipt), Form("#it{M}_{rec}/#it{M}_{gen} bin %.0f<#it{p}_{T}<%.0f GeV/#it{c};#it{M}_{rec}; Entries", ptlims[ipt], ptlims[ipt+1]), 20, 0., 20./*, "s"*/);
      
      hMRpTRange[ipt] = new TH1F(Form("hMRpTRange_%d", ipt), Form("#it{M}_{rec} bin %.0f<#it{p}_{T}<%.0f GeV/#it{c};#it{M}_{rec}; Entries", ptlims[ipt], ptlims[ipt+1]), 60, -20, 100);
      hMRpTRange[ipt]->SetLineColor(colors[0]);
      hMRpTRange[ipt]->GetXaxis()->SetRangeUser(-5, 20);
      
      hMPpTRange[ipt] = new TH1F(Form("hMPpTRange_%d", ipt), Form("#it{M}_{rec} bin %.0f<#it{p}_{T}<%.0f GeV/#it{c};#it{M}_{rec}; Entries", ptlims[ipt], ptlims[ipt+1]), 60, -20, 100);
      hMPpTRange[ipt]->SetLineColor(colors[1]);
      hMPpTRange[ipt]->GetXaxis()->SetRangeUser(-5, 20);
   }

   TLegend *legMass = new TLegend(0.5, 0.6, 0.8, 0.8);
   legMass->SetFillStyle(0);
   legMass->AddEntry(hMRpTRange[0], "Detector", "L");
   legMass->AddEntry(hMPpTRange[0], "Particle", "L");
   
   Int_t nentries = tree->GetEntries();
   const Int_t nbranches = tree->GetNbranches();
   Printf("Number of branches %d", nbranches);
   
   TLorentzVector *vecsD = 0x0; // 0 reco
   TLorentzVector *vecsP = 0x0; // 1 particle
   
   
   TBranch *bJD = 0x0;
   TBranch *bJP = 0x0;
   //for(Int_t k = 0; k< nbranches; k++)
   //for(Int_t k = 0; k< 1; k++){
   //bJ[k] = 0;
   Printf("Set Branch Address %p, %s,  %p %p",tree, names[0].Data(), vecsD, bJD);
   Int_t checkBranch = tree->SetBranchAddress(Form("%s", names[0].Data()),&vecsD, &bJD);
   if(checkBranch< 0) {
      Printf("Wrong");
      
   }
   tree->SetBranchAddress(Form("%s", names[1].Data()), &vecsP, &bJP);
   
   //}
   
   Printf("Entries =  %d", nentries);
   if(!bJD){
      Printf("No branch %s found", names[0].Data());
      return;
   }   
   for(Int_t i = 0; i<nentries; i++){
   //for(Int_t i = 0; i<50000; i++){
      //Printf("GetEntry %d", i);
      //Long64_t tentry = tree->LoadTree(i);
      //cout<<"tentry "<< tentry<<" and bJ[0] = "<< bJ[0]<<endl;
      //Printf("tentry = %lld", tentry);
      Int_t ient = bJD->GetEntry(i);
      //Printf("i entry = %d", ient);
      Double_t ptgen, ptrec, mrec, mgen;
      ptrec = vecsD->Pt();
      mrec = vecsD->M();
      if(applyptcut && ptrec >= ptlims[0] && ptrec <= ptlims[nptbins]) continue;
      //if(ptrec < 10) continue;
      hpTR->Fill(ptrec);
      hMR->Fill(mrec);
      ient = bJP->GetEntry(i);
      //if(vecsD) Printf("Read: pT %f, eta %f, phi %f, m  %f", vecsD->Pt(), vecsD->Eta(), vecsD->Phi(), vecsD->M());

      //Printf("i entry = %d", ient);
      ptgen = vecsP->Pt();
      mgen = vecsP->M();
      
      hpTP->Fill(ptgen);
      hMP->Fill(mgen);
      
      hResponsepT->Fill(ptgen, ptrec);
      hResponseM->Fill(mgen, mrec);
      
      hdpT->Fill(ptrec - ptgen);
      hdM->Fill(mrec - mgen);
      
      Int_t ptbin = 0;
      ptbin = FindBinInArray(ptgen, ptlims, nptbins);
      if(ptbin>-1) {
      	 hMgenOverMrec[ptbin]->Fill(mrec, mrec/mgen);
      	 hMRpTRange[ptbin]->Fill(mrec);
      	 hMPpTRange[ptbin]->Fill(mgen);
      }
   }

   TH1D* hpThisto = 0x0;
   //TH1D* hpThisto = DrawSampleTask("AnalysisResults.root", "EmcalJetSampleCheckDistr_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TPC_histos");
   if(0){
      hpThisto = DrawSampleTask("AnalysisResults.root", "PrepareInputForEmbedding", "fhFractionSharedpT");
      hpThisto->Scale(1./hpThisto->Integral("width"));
      hpThisto->SetLineColor(colors[5]);
   }
   TCanvas *cpT = new TCanvas(Form("cpT%s", suffix.Data()), Form("cpT%s", suffix.Data()), 600, 600);
   cpT->cd();
   gPad->SetLogy();
   hpTR->Scale(1./hpTR->Integral("width"));
   hpTP->Scale(1./hpTP->Integral("width"));
   hpTR->Draw("E");
   hpTP->Draw("Esames");
   if(hpThisto) hpThisto->Draw("Esames");
   if(hpTDEmb) {
      hpTDEmb->Scale(hpTR->Integral()/hpTDEmb->Integral());
      hpTDEmb->SetMarkerColor(colors[0]);
      hpTDEmb->SetMarkerStyle(25);
      hpTDEmb->Draw("Esames");
   }
   if(hpTPEmb) {
      hpTPEmb->Scale(hpTP->Integral()/hpTPEmb->Integral());
      hpTPEmb->SetMarkerColor(colors[1]);
      hpTPEmb->SetMarkerStyle(21);
      hpTPEmb->Draw("Esames");
   }
   TCanvas *cM = new TCanvas(Form("cM%s", suffix.Data()), Form("cM%s", suffix.Data()), 600, 600);
   cM->cd();
   gPad->SetLogy();
   hMR->Draw("E");
   hMP->Draw("Esames");
   if(hMDEmb) {
      hMDEmb->Scale(hMR->Integral()/hMDEmb->Integral());
      hMDEmb->SetMarkerStyle(25);
      hMDEmb->SetMarkerColor(colors[0]);
      hMDEmb->Draw("Esames");
   }
   if(hMPEmb) {
      hMPEmb->Scale(hMP->Integral()/hMPEmb->Integral());
      hMPEmb->SetMarkerStyle(21);
      hMPEmb->SetMarkerColor(colors[1]);
      hMPEmb->Draw("Esames");
   }
   TCanvas *cpTResp = new TCanvas("cpTResp", "Response pt");
   cpTResp->cd();
   gPad->SetLogz();
   hResponsepT->Draw("colz");
   
   TCanvas *cMResp = new TCanvas("cMResp", "Response M");
   cMResp->cd();
   gPad->SetLogz();
   hResponseM->Draw("colz");
   
   TCanvas *cdpT = new TCanvas(Form("cdpT%s", suffix.Data()), Form("cdpT%s", suffix.Data()), 600, 600);
   cdpT->cd();
   hdpT->Scale(1./hdpT->Integral());
   hdpT->Draw("E");
   if(hdpTEmb) {
      hdpTEmb->Scale(1./hdpTEmb->Integral());
      hdpTEmb->SetMarkerStyle(29);
      hdpTEmb->Draw("Esames");
   }
   
   TCanvas *cdM = new TCanvas(Form("cdM%s", suffix.Data()), Form("cdM%s", suffix.Data()), 600, 600);
   cdM->cd();
   hdM->Scale(1./hdM->Integral());
   hdM->Draw("E");
   if(hdMEmb) {
      hdMEmb->Scale(1./hdMEmb->Integral());
      hdMEmb->SetMarkerStyle(29);
      hdMEmb->Draw("Esames");
      legDelta->Draw();
   }
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cFractionM = new TCanvas("cFractionM", "Fraction Mgen/Mrec", dx, dy);
   cFractionM->Divide(nx, ny);
   
   TCanvas *cMpTRanges = new TCanvas("cMpTRanges", "Mgen, Mrec in pt ranges", dx, dy);
   cMpTRanges->Divide(nx, ny);
   
   TFile *fout = new TFile("HistosFromTree.root", "recreate");
   
   for(Int_t ipt = 0; ipt<nptbins; ipt++){
      cFractionM->cd(ipt+1);
      hMgenOverMrec[ipt]->Draw("E");
      fout->cd();
      hMgenOverMrec[ipt]->Write();
      cMpTRanges->cd(ipt+1);
      hMRpTRange[ipt]->Draw("E");
      hMPpTRange[ipt]->Draw("Esames");
      legMass->Draw();
   }  

   fout->cd();
   hpTP          ->Write();
   hpTR          ->Write();
   hMR           ->Write();
   hMP           ->Write();
   hResponsepT   ->Write();
   hResponseM    ->Write();
   hdpT          ->Write();
   hdM           ->Write();
   
   SaveCv(cpT);
   SaveCv(cM);
   SaveCv(cpTResp);
   SaveCv(cMResp);
   SaveCv(cFractionM);
   SaveCv(cdpT);
   SaveCv(cdM);
}

void LoopOverTreeWithDoubles(TTree *tree, TString names[], TString suffix){

   //Warning!! Assume that the tree contains Double_t
   
   if(!tree){
      Printf("Tree null");
      return;
   
   }

   TH1F *hpTR = new TH1F(Form("hpTR_%s", suffix.Data()), Form("hpTR_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hpTR->SetLineWidth(2);
   hpTR->SetLineColor(colors[0]);
   TH1F *hpTP = new TH1F(Form("hpTP_%s", suffix.Data()), Form("hpTP_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   
   
   TH1F *hMR = new TH1F(Form("hMR_%s", suffix.Data()), Form("hMR_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hMR->SetLineWidth(2);
   hMR->SetLineColor(colors[0]);
   TH1F *hMP = new TH1F(Form("hMP_%s", suffix.Data()), Form("hMP_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hMP->SetLineWidth(2);
   hMP->SetLineColor(colors[1]);
   
   TH2F *hResponsepT = new TH2F("hResponsepT", "Response pT; #it{p}_{T,gen};#it{p}_{T,rec};", 60, -10, 100, 60, -10, 100);
   TH2F *hResponseM  = new TH2F("hResponseM", "Response M; #it{M}_{gen};#it{M}_{rec};", 60, -10, 20, 60, -10, 20);
   
   //bin-by-bin correction 
   Int_t nptbins = 4;
   Double_t ptlims[nptbins+1] = {20, 40, 60, 80, 100};
   TProfile *hMgenOverMrec[nptbins];
   for(Int_t ipt = 0; ipt<nptbins; ipt++){
      hMgenOverMrec[ipt] = new TProfile(Form("hMgenOverMrec_%d", ipt), Form("#it{M}_{rec}/#it{M}_{gen} bin %.0f<#it{p}_{T}<%.0f GeV/#it{c};#it{M}_{rec}; Entries", ptlims[ipt], ptlims[ipt+1]), 20, 0., 20./*, "s"*/);
   }

   
   Int_t nentries = tree->GetEntries();
   const Int_t nbranches = tree->GetNbranches();
   Printf("Number of branches %d", nbranches);

   Double_t pt, eta, phi, m;
   Double_t ptp, etap, phip, mp;   
   tree->SetBranchAddress(names[0],&pt);
   tree->SetBranchAddress(names[1],&eta);
   tree->SetBranchAddress(names[2],&phi);
   tree->SetBranchAddress(names[3],&m);
   tree->SetBranchAddress(names[4],&ptp);
   tree->SetBranchAddress(names[5],&etap);
   tree->SetBranchAddress(names[6],&phip);
   tree->SetBranchAddress(names[7],&mp);
   
   //}
   
   Printf("Entries =  %d", nentries);

   for(Int_t i = 0; i<nentries; i++){
      Int_t ient = tree->GetEntry(i);
      //Printf("i entry = %d", ient);
      
      hpTR->Fill(pt);
      hMR->Fill(m);

      
      hpTP->Fill(ptp);
      hMP->Fill(mp);
      
      hResponsepT->Fill(ptp, pt);
      hResponseM->Fill(mp, m);
      
      Int_t ptbin = 0;
      ptbin = FindBinInArray(ptp, ptlims, nptbins);
      if(ptbin>-1) hMgenOverMrec[ptbin]->Fill(m, m/mp);
   }
   //TH1D* hpThisto = DrawSampleTask("AnalysisResults.root", "EmcalJetSampleCheckDistr_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TPC_histos");
   TH1D* hpThisto = DrawSampleTask("AnalysisResults.root", "PrepareInputForEmbedding", "fhFractionSharedpT");
   hpThisto->Scale(1./hpThisto->Integral("width"));
   hpThisto->SetLineColor(colors[5]);
   
   TCanvas *cpT = new TCanvas(Form("cpT%s", suffix.Data()), Form("cpT%s", suffix.Data()), 600, 600);
   cpT->cd();
   hpTR->Scale(1./hpTR->Integral("width"));
   hpTP->Scale(1./hpTP->Integral("width"));
   hpTR->Draw();
   hpTP->Draw("sames");
   if(hpThisto) hpThisto->Draw("sames");
   
   TCanvas *cM = new TCanvas(Form("cM%s", suffix.Data()), Form("cM%s", suffix.Data()), 600, 600);
   cM->cd();
   hMR->Draw();
   hMP->Draw("sames");
   
   TCanvas *cpTResp = new TCanvas("cpTResp", "Response pt");
   cpTResp->cd();
   hResponsepT->Draw("colz");
   
   TCanvas *cMResp = new TCanvas("cMResp", "Response M");
   cMResp->cd();
   hResponseM->Draw("colz");
   

   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cFractionM = new TCanvas("cFractionM", "Fraction Mgen/Mrec", dx, dy);
   cFractionM->Divide(nx, ny);
   for(Int_t ipt = 0; ipt<nptbins; ipt++){
      cFractionM->cd(ipt+1);
      hMgenOverMrec[ipt]->Draw();
   }
}

//____________________________________________________________________
void TestTree(){
   
   TLorentzVector fJetDet(0,0,0,0); //TLorentzVector();
   TTree *fTreeJets = new TTree("fTreeJet", "fTreeJet");
   fTreeJets->Branch("fJetDet.",  &fJetDet);
   //fTreeJets->Branch("fJetPart",&fJetPart);
   
   TRandom3 *rndm = new TRandom3(1234);
   TFile *fout = new TFile("output.root", "update");
   Int_t nevents = 1000;
   for(Int_t i = 0; i<nevents; i++){
      
      Double_t pt = rndm->Uniform(20, 100);
      Double_t eta = rndm->Uniform(-1, 1);
      Double_t phi = rndm->Uniform(0, 2*TMath::Pi());
      Double_t m = rndm->Uniform(1,10);
      fJetDet.SetPtEtaPhiM(pt, eta, phi, m);
      
      if(i>0){
      	 TTree *readTree = dynamic_cast<TTree*>(fout->Get("fTreeJet"));
      	 readTree->SetBranchAddress("fJetDet.",  &fJetDet);
      	 readTree->Fill();
      	 readTree->Write();
      } else {
      	 
      	 fTreeJets->Fill();
      	 fTreeJets->Write();
      }
   }
   
}

//____________________________________________________________________

void ScaleFactors(Int_t neventsLHC10b, TString filename, TString treename){
   
   TTree* tree = ReadTree(filename, treename);
   TString branches[2] = {"fJetDet.", "fJetPart."};
   
   TLorentzVector *vecD = 0x0;
   tree->SetBranchAddress(branches[0], &vecD);
   
   Int_t nptbins = 9;
   Float_t ptlims[nptbins] = {5.,15., 25., 40., 55.,60., 75., 90., 105.}; // these are the low limitf of the pT bins
   
   
   Float_t nembedperBin = neventsLHC10b/(Float_t)nptbins;
   
   Float_t entriesTreePerBin[nptbins];
   Float_t factorPerBin[nptbins];
   
   Printf("pT ranges considered");
   for(Int_t i = 0; i<nptbins; i++ ) {
      if(i< nptbins-1) Printf(" * %.1f < #it{p}_{T} < %.1f GeV/c", ptlims[i], ptlims[i+1]);
      else Printf(" * #it{p}_{T} > %.1f GeV/c", ptlims[i]);
   }
   Printf("I want %f embedded particles per pT range", nembedperBin);

   tree->Draw(Form("%s.Pt()", branches[0].Data()));
   TH1F *hpT = (TH1F*)gPad->GetPrimitive("htemp");
   hpT->SetName("hpT");
   
   Printf("Factors");
   for(Int_t i = 0; i < nptbins; i++){
      Int_t bin[2];
      bin[0] = hpT->FindBin(ptlims[i]);
      bin[1] = (i< nptbins-1) ? (hpT->FindBin(ptlims[i+1]) - 1) : hpT->GetNbinsX();
      entriesTreePerBin[i] = hpT->Integral(bin[0], bin[1]);
      factorPerBin[i] = nembedperBin/entriesTreePerBin[i];
      Printf("factor[%d] = %.3e;", i, factorPerBin[i]);
   }

}


//_____________________________________________________________________________

void SelectPtTree(Float_t ptmin = 10, Float_t ptmax = 200, Int_t whichpt = 1, TString filename = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1054/outputpTHardBins/mergeRuns/plainMerge/out.root", TString treename = "fTreeJet"){
   TStopwatch watch;
   watch.Start();
   TTree* tree = ReadTree(filename, treename);
   if(!tree){
      Printf("Tree not found, exit");
      return;
   }
   //tree->LoadBaskets(13ul * 1024 * 1024 * 1024);
   
   TLorentzVector *vecD = 0x0;
   TLorentzVector *vecP = 0x0;
   TString branches[2] = {"fJetDet.", "fJetPart."};
   tree->SetBranchAddress(branches[0], &vecD);
   tree->SetBranchAddress(branches[1], &vecP);
   
   TString types = "Det";
   if(whichpt == 2) types = "Par";
   
   TFile *fOut = TFile::Open(Form("outPt%s%.0f-%.0f.root", types.Data(), ptmin, ptmax), "recreate");
   
   TTree *outTree = new TTree(Form("fTreeJetPt%s%.0f%.0f", types.Data(), ptmin, ptmax), Form("fTreeJetPt%.0f-%.0f", ptmin, ptmax));
   TLorentzVector *vecDN = 0x0;
   TLorentzVector *vecPN = 0x0;
   outTree->Branch(branches[0], &vecDN);
   outTree->Branch(branches[1], &vecPN);
   
   const Int_t nEntries = tree->GetEntries();
    
   TRandom3 *rnd = new TRandom3(123);
   for(Int_t i = 0; i < nEntries; i++){
      //Printf("Entry %d", i);
      if (i % 100000 == 0) Printf("Entry %d", i);
      Int_t entry = tree->GetEntry(i);
      Float_t pt = 0;
      switch (whichpt){
      case 1:
      	 pt = vecD->Pt();
      	 break;
      	 
      case 2:
      	 pt = vecP->Pt();
      	 break;
      }
      
      if(pt>ptmin && pt<ptmax){
      	 vecDN = dynamic_cast<TLorentzVector*>(vecD->Clone());
      	 vecPN = dynamic_cast<TLorentzVector*>(vecP->Clone());
      	 outTree->Fill();
      }
      
   }
   Printf("Done, Saving now...");
   fOut->WriteTObject(outTree);
   fOut->Close();
   Printf("Done.");
   
   watch.Stop();
   watch.Print();
}

//_____________________________________________________________________________

TH1F* GetHistoEmb(TString hname, TList *list){
   if(!list){
      Printf("List is null");
      return 0x0;
   }
   TH1F *h = dynamic_cast<TH1F*>(list->FindObject(hname));
   if(!h) {
      Printf("%s was not found", hname.Data());
   }
   return h;
}

TH1D* GetHistoEmb(TString hname, TList *list, Int_t axisPj, Int_t* axisSel, Float_t* cutL, Float_t* cutH){
   if(!list){
      Printf("List is null");
      return 0x0;
   }
   THnSparseF *hnsp = dynamic_cast<THnSparseF*>(list->FindObject(hname));
   if(!hnsp) {
      Printf("%s was not found", hname.Data());
      return 0;
   }
   Int_t ndim = hnsp->GetNdimensions();
   if(axisSel && cutL && cutH){
      
      for(Int_t i = 0; i<ndim; i++){
      	 if(i == axisPj) continue;
      	 Int_t range[2] = {hnsp->GetAxis(i)->FindBin(cutL[i]), hnsp->GetAxis(i)->FindBin(cutL[i]) - 1};
      	 hnsp->GetAxis(i)->SetRange(range[0], range[1]);
      	 
      }
   }
   if(axisPj<0 || axisPj>=ndim){ 
      Printf("Axis %d not available", axisPj);
      return 0;
   }
   TH1D *h = hnsp->Projection(axisPj);
   Printf("Projection %s", h->GetName());
   return h;
}

//_____________________________________________________________________________

void WeightTreeDistributions(TString inputfile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1067/output/AnalysisResults.root", TString treename = "fTreeJet", Bool_t useweights = kTRUE, TString suffix = "", Float_t pTrecCut = 10., Bool_t drawEmbTask = kFALSE, TString fileEmb = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1068/analysis/AnalysisResults.root", TString listEmbN = "SingleTrackEmbedding"){
   
   // drawEmbTask =  kTRUE overlaps the pT and mass distribution from the embedding task
   // the weighting done here is wrong I think
   
   TFile *fout = new TFile("HistosFromTree.root", "recreate");
   
   TTree* tree = ReadTree (inputfile, treename);
   if(!tree) return;
   const Int_t nbranches = 5;
   TString branches[nbranches] = {"fJetDet.", "fJetPart.", "fXsection", "fNTrials", "fPthardB"};
   TLorentzVector *vecD = 0x0;
   TLorentzVector *vecP = 0x0;
   Float_t xsec = 0, trials = 0 ;
   Int_t pthardbin = 0;
   
   tree->SetBranchAddress(branches[0], &vecD);
   tree->SetBranchAddress(branches[1], &vecP);
   tree->SetBranchAddress(branches[2], &xsec);
   tree->SetBranchAddress(branches[3], &trials);
   tree->SetBranchAddress(branches[4], &pthardbin);
   
   TH1F *hpTR = new TH1F(Form("hpTR_%s", suffix.Data()), Form("hpTR_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hpTR->SetLineWidth(2);
   hpTR->SetLineColor(colors[0]);
   hpTR->Sumw2();
   TH1F *hpTRUW = new TH1F(Form("hpTRUW_%s", suffix.Data()), Form("hpTRUW_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hpTRUW->SetLineWidth(2);
   hpTRUW->SetLineColor(colors[0]);
   hpTRUW->SetMarkerStyle(30);
   hpTRUW->SetMarkerColor(colors[0]);
   hpTRUW->Sumw2();
   
   TH1F *hpTP = new TH1F(Form("hpTP_%s", suffix.Data()), Form("hpTP_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hpTP->SetLineColor(colors[1]);
   hpTP->SetLineWidth(2);
   hpTP->Sumw2();
   
   TH1F *hMR = new TH1F(Form("hMR_%s", suffix.Data()), Form("hMR_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hMR->SetLineWidth(2);
   hMR->SetLineColor(colors[0]);
   hMR->Sumw2();
   TH1F *hMP = new TH1F(Form("hMP_%s", suffix.Data()), Form("hMP_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), 60, -20, 100);
   hMP->SetLineWidth(2);
   hMP->SetLineColor(colors[1]);
   hMP->Sumw2();
   
   TH1F *hxsec = new TH1F("hxsec", "Cross section; x-sec;", 1e3, 1e-8, 2e-6);
   TH1F *htria = new TH1F("htria", "Trials; ntrials;", 1e3, 1, 3e3);
   TH1F *hNpthB= new TH1F("hNpthB", "Entries per pthard bin; pT-hard Bin", 10, 0.5, 10.5 );
   TH1F *hweig = new TH1F("hweig", "Weight per jet; weight;", 1e3, 1e-10, 1e-8);
   
   // thnsparse with the axes like MassResponse Task
   const Int_t naxes = 4;
   Double_t minRange1[naxes] = {0., 0., 0., 0.};
   Double_t maxRange1[naxes] = {50., 50., 200., 200.};
   Int_t nbins1[naxes] = {100, 100, 200, 200};
   TString titles = "Response; #it{M}_{det};#it{M}_{par}; #it{p}_{T, det}; #it{p}_{T, par}"; 
   THnSparseF *hMassRespTree = new THnSparseF("hMassRespTree", titles, naxes, nbins1, minRange1, maxRange1);
   hMassRespTree->Sumw2();
   
   // thnsparse with the axes like JetShape Task
   Double_t minRange2[naxes] = {-20., -20., -50., -50.};
   Double_t maxRange2[naxes] = {80., 80., 150., 150.};
   Int_t nbins2[naxes] = {100, 100, 200, 200};
   THnSparseF *hShapeRespTree = new THnSparseF("hShapeRespTree", titles, naxes, nbins2, minRange2, maxRange2);
   hShapeRespTree->Sumw2();
   
   TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.8);
   leg->SetFillStyle(0);
   leg->AddEntry(hMR, "Detector", "L");
   leg->AddEntry(hMP, "Particle", "L");
   
   TString sw = "width";
   
   const Int_t nEntries = tree->GetEntries();
   
   for(Int_t i = 0; i < nEntries; i++){
      Int_t entry = tree->GetEntry(i);
      hNpthB->Fill(pthardbin);
   }
   
   TCanvas *cpthardB = new TCanvas("cpthardB", "Pt hard bins entries", 400, 400);
   cpthardB->cd();
   hNpthB->Draw();
   
   for(Int_t i = 0; i < nEntries; i++){
      
      Int_t entry = tree->GetEntry(i);
      if(vecD->Pt() < pTrecCut) continue;
      Double_t values[naxes] = {vecD->M(), vecP->M(), vecD->Pt(), vecP->Pt()};
      if(useweights){
      	 
      	 Double_t pthBEntries = hNpthB->GetBinContent(hNpthB->FindBin(pthardbin));
      	 if (i % 100000 == 0) Printf("Entry %d: This is pt hard %d, entries %f", i, pthardbin, pthBEntries);
      	 Double_t ww = 0;
      	 if(trials>0 && pthBEntries > 0) ww = xsec/trials/pthBEntries;
      	 hxsec->Fill(xsec);
      	 htria->Fill(trials);
      	 hweig->Fill(ww);
      	 
      	 hpTR->Fill(vecD->Pt(), ww);
      	 hpTRUW->Fill(vecD->Pt());
      	 hpTP->Fill(vecP->Pt(), ww);
      	 
      	 hMR->Fill(vecD->M(), ww);
      	 hMP->Fill(vecP->M(), ww);
      	 
      	 hMassRespTree-> Fill(values, ww);
      	 hShapeRespTree->Fill(values, ww);
      	 
      } else {
      	 hpTR->Fill(vecD->Pt());
      	 hpTRUW->Fill(vecD->Pt());
      	 hpTP->Fill(vecP->Pt());
      	 
      	 hMR->Fill(vecD->M());
      	 hMP->Fill(vecP->M());
      	 
      	 hMassRespTree-> Fill(values);
      	 hShapeRespTree->Fill(values);
      }
   }
   
   fout->cd();
   hMassRespTree-> Write();
   hShapeRespTree->Write();
   
   TCanvas *cpT = new TCanvas(Form("cpT%s", suffix.Data()), Form("cpT%s", suffix.Data()), 600, 600);
   cpT->cd();
   gPad->SetLogy();
   Double_t refIntpT = hpTR->Integral(sw);
   hpTRUW->Scale(refIntpT/hpTRUW->Integral(sw));
   hpTP->Scale(refIntpT/hpTP->Integral(sw));
   hpTR->Draw("E");
   hpTP->Draw("Esames");
   hpTRUW->Draw("Esames");
   leg->Draw();
   
   TCanvas *cM = new TCanvas(Form("cM%s", suffix.Data()), Form("cM%s", suffix.Data()), 600, 600);
   cM->cd();
   gPad->SetLogy();
   hMR->Draw("E");
   hMP->Draw("Esames");
   leg->Draw();

   TCanvas *cweights = new TCanvas("cweights", "xsec and trials", 1200, 400);
   cweights->Divide(3,1);
   cweights->cd(1);
   hxsec->Draw();
   cweights->cd(2);
   htria->Draw();
   cweights->cd(3);
   hweig->Draw();
   
   if(drawEmbTask){
      TList *listEmb = ReadFile(fileEmb, listEmbN);
      
      TString hName = "fhpTEmb";                
      TH1F* hpTDEmb = GetHistoEmb(hName, listEmb);
      if(hpTDEmb){
      	 hpTDEmb->SetMarkerColor(colors[0]);
      	 hpTDEmb->SetMarkerStyle(20);
      	 hpTDEmb->Scale(refIntpT/hpTDEmb->Integral(sw));
      	 cpT->cd();
      	 hpTDEmb->Draw("Esames");
      }
      
      hName = "fhMEmb";            
      TH1F* hMDEmb = GetHistoEmb(hName, listEmb);
      if(hMDEmb){
      	 hMDEmb->SetMarkerColor(colors[0]);
      	 hMDEmb->SetMarkerStyle(20);
      	 hMDEmb->Scale(refIntpT/hMDEmb->Integral(sw));
      	 cM->cd();
      	 hMDEmb->Draw("Esames");
      }
      
      hName = "fhPartJet";
      //Int_t axesSel[3] = {0, 1, 2};
      Int_t projax = 3; //pT
      TH1D *hpTPEmb = GetHistoEmb(hName, listEmb, projax);
      if(hpTPEmb){
      	 hpTPEmb->SetMarkerColor(colors[1]);
      	 hpTPEmb->SetMarkerStyle(20);
      	 hpTPEmb->Scale(refIntpT/hpTPEmb->Integral(sw));
      	 cpT->cd();
      	 hpTPEmb->Draw("Esames");
      	 projax = 2; //M
      	 TH1D *hMPEmb = GetHistoEmb(hName, listEmb, projax);
      	 hMPEmb->SetMarkerColor(colors[1]);
      	 hMPEmb->SetMarkerStyle(20);
      	 hMPEmb->Scale(refIntpT/hMPEmb->Integral(sw));
      	 cM->cd();
      	 hMPEmb->Draw("Esames");
      }
      
      hName = "fhnMassResponse";
      projax = 3; //pT part
      TH1D *hpTP576 = GetHistoEmb(hName, listEmb, projax);
      if(hpTP576){
      	 hpTP576->SetMarkerColor(colors[1]);
      	 hpTP576->SetMarkerStyle(20);
      	 hpTP576->Scale(refIntpT/hpTP576->Integral(sw));
      	 cpT->cd();
      	 hpTP576->Draw("Esames");
      	 projax = 1; //M
      	 TH1D *hMP576 = GetHistoEmb(hName, listEmb, projax);
      	 hMP576->SetMarkerColor(colors[1]);
      	 hMP576->SetMarkerStyle(20);
      	 hMP576->Scale(refIntpT/hMP576->Integral(sw));
      	 cM->cd();
      	 hMP576->Draw("Esames");
      	 
      }
      
      
      
      
   }
}

void DrawEntriesUsed(TString fileEmb, TString listEmb, TString listResp, TString fileTreeIn, Double_t embPtCut = 0, Bool_t savetree = kFALSE, Int_t pTHB = 1, Bool_t fillTreeNtimes = kFALSE){
  
   //cross check the expected distribution of the embedded tracks
   /// fileEmb = path and filename of the embedding job
   /// listEmb = list name in the embedding output -> e.g. "SingleTrackEmbedding9"
   /// listResp= list name of the response matrix  -> e.g. "JetShapeConst_JetRhosub_AKTChargedR040_PicoTracksPtH8"
   /// fileTreeIn = path to the randomized Tree used as input to the embedding
   /// embPtCut = pT particle level cut applied to the embedding task. If it's the same as the tree production, it can be 0 = default
   /// * It was particle level, changed to pT particle level, because this has changed in the embedding (20160325)
   /// savetree = kTRUE, if you want to save the tree with only the used entries
   /// need to specify the pT hard bin number, for directory (must be already created) and file name
   /// fillTreeNtimes: if kTRUE, the histograms filled with the tree entries are filled the number of times the entry was used. If kFALSE, the histograms from the tree entries are normalised to the integral of the response projections
   
   TString suffix = "";
   TList* list = ReadFile(fileEmb, listEmb);
   if(!list) return;
   
   TH1F *entriesIn = (TH1F*)list->FindObject("fhTreeEntriesUsed");
   if(!entriesIn) {
      Printf("File not found");
      return;
   }
   
   TTree* tree = ReadTree (fileTreeIn, "fTreeJetR");
   if(!tree) {
      Printf("Tree not found");
      return;
   }
   
   TLorentzVector *vecD = 0x0;
   TLorentzVector *vecP = 0x0;
   TString branches[2] = {"fJetDet.", "fJetPart."};
   
   TFile* fout = 0x0;
   if(savetree) //fout = new TFile(Form("%2d/TreeEntriesUsed%d.root", pTHB, pTHB), "recreate"); //need to make the directory first
   fout = new TFile(Form("TreeEntriesUsed%d.root", pTHB), "recreate"); 

   TTree* treeNEntries = new TTree("fTreeNEntries", "Tree with only entries used");
   treeNEntries->Branch(branches[0], &vecD);
   treeNEntries->Branch(branches[1], &vecP);
   
   tree->SetBranchAddress(branches[0], &vecD);
   tree->SetBranchAddress(branches[1], &vecP);
   
   TList* listRs = ReadFile(fileEmb, listResp);
   if(!listRs) return;
   
   // histograms from the embedding list
   TH1F *hpTEmbR = (TH1F*)list->FindObject("fhpTEmb");
   TH1F *hMEmbR = (TH1F*)list->FindObject("fhMEmb");
   

   hpTEmbR->SetMarkerStyle(21);
   hpTEmbR->SetMarkerColor(colors[5]);
   hMEmbR->SetMarkerStyle(21);
   hMEmbR->SetMarkerColor(colors[5]);
   
   THnSparseF *hspEmb = (THnSparseF*)list->FindObject("fhPartJet");
   TH1F *hMEmbP  = (TH1F*)hspEmb->Projection(2);
   TH1F *hpTEmbP = (TH1F*)hspEmb->Projection(3);
   
   
   hpTEmbP->SetMarkerStyle(21);
   hpTEmbP->SetMarkerColor(colors[5]);
   hMEmbP->SetMarkerStyle(21);
   hMEmbP->SetMarkerColor(colors[5]);

   TH1  *hMEmbPForDiv  =  0x0;
   TH1  *hpTEmbPForDiv =  0x0;
   TH1  *hMEmbRForDiv  =  0x0;
   TH1  *hpTEmbRForDiv =  0x0;
   
   // histograms from the response matrix
   THnSparseF *hnsp = (THnSparseF*)listRs->FindObject("fhnMassResponse_0");
   //selection on pT emb
   Int_t axMSub   = 0;
   Int_t axMTru   = 1;
   Int_t axPtSub  = 2;
   Int_t axPtTru  = 3;
   Int_t axMEmb   = 4;
   Int_t axPtEmb  = 5;
   
   //old
   //Int_t ptB = hnsp->GetAxis(axPtEmb)->FindBin(embPtCut), maxAx = hnsp->GetAxis(axPtEmb)->GetNbins();
   //hnsp->GetAxis(axPtEmb)->SetRange(ptB, maxAx);
   Int_t ptB = hnsp->GetAxis(axPtTru)->FindBin(embPtCut), maxAx = hnsp->GetAxis(axPtTru)->GetNbins();
   hnsp->GetAxis(axPtTru)->SetRange(ptB, maxAx);

   TH1D *horigRespMP  = hnsp->Projection(axMTru);                
   TH1D *horigRespPtP = hnsp->Projection(axPtTru);
   TH1D *horigRespMEmb = hnsp->Projection(axMEmb);
   TH1D *horigRespPtEmb= hnsp->Projection(axPtEmb);
   TH1D *horigRespMDF  = hnsp->Projection(axMSub);
   TH1D *horigRespPtDF = hnsp->Projection(axPtSub);
   TH1  *hRespMP  =  0x0;
   TH1  *hRespPtP =  0x0;
   TH1  *hRespMEmb = 0x0;
   TH1  *hRespPtEmb= 0x0;
   TH1  *hRespMDF  = 0x0;
   TH1  *hRespPtDF = 0x0;
   
   //rebinning / re-ranging
   UniformTH1FForDivide(horigRespMP, hMEmbP, hRespMP, hMEmbPForDiv);
   UniformTH1FForDivide(horigRespPtP, hpTEmbP, hRespPtP, hpTEmbPForDiv);
   UniformTH1FForDivide(horigRespMEmb, hMEmbR, hRespMEmb, hMEmbRForDiv);
   UniformTH1FForDivide(horigRespPtEmb, hpTEmbR, hRespPtEmb, hpTEmbRForDiv);
   hMEmbRForDiv = 0x0; hpTEmbRForDiv = 0x0;
   UniformTH1FForDivide(horigRespMDF, hMEmbR, hRespMDF, hMEmbRForDiv);
   UniformTH1FForDivide(horigRespPtDF, hpTEmbR, hRespPtDF, hpTEmbRForDiv);

   //scale to the same integral
   
   //hpTEmbPForDiv->Scale(1./hpTEmbPForDiv->Integral());
   //hMEmbPForDiv-> Scale(1./hMEmbPForDiv ->Integral());
   //hpTEmbRForDiv->Scale(1./hpTEmbRForDiv->Integral());
   //hMEmbRForDiv-> Scale(1./hMEmbRForDiv ->Integral());
   //
   //hRespMP    ->Scale(1./hRespMP   ->Integral());
   //hRespPtP   ->Scale(1./hRespPtP  ->Integral());
   //hRespMEmb  ->Scale(1./hRespMEmb ->Integral());
   //hRespPtEmb ->Scale(1./hRespPtEmb->Integral());
   //hRespMDF   ->Scale(1./hRespMDF  ->Integral());
   //hRespPtDF  ->Scale(1./hRespPtDF ->Integral());
   Double_t integMP  = hRespMP   ->Integral();
   Double_t integPtP = hRespPtP  ->Integral();
   Double_t integMR  = hRespMEmb ->Integral();
   Double_t integPtR = hRespPtEmb->Integral();
   
   hRespMP    ->SetMarkerStyle(24);
   hRespPtP   ->SetMarkerStyle(24);
   hRespMEmb  ->SetMarkerStyle(24);
   hRespPtEmb ->SetMarkerStyle(24);
   hRespMDF   ->SetMarkerStyle(24);
   hRespPtDF  ->SetMarkerStyle(24);
   
   hRespMP    ->SetMarkerColor(kBlue);
   hRespPtP   ->SetMarkerColor(kBlue);
   hRespMEmb  ->SetMarkerColor(kBlue);
   hRespPtEmb ->SetMarkerColor(kBlue);
   hRespMDF   ->SetMarkerColor(kRed);
   hRespPtDF  ->SetMarkerColor(kRed);

   Int_t nbinsMP = hRespMP->GetNbinsX(), nbinspTP = hRespPtP->GetNbinsX();
   Double_t minMP = hRespMP->GetBinLowEdge(1), minpTP = hRespPtP->GetBinLowEdge(1), maxMP = hRespMP->GetBinLowEdge(nbinsMP+1), maxpTP = hRespPtP->GetBinLowEdge(nbinspTP+1);
   Int_t nbinsMR = hRespMEmb->GetNbinsX(), nbinspTR = hRespPtEmb->GetNbinsX();
   Double_t minMR = hRespMEmb->GetBinLowEdge(1), minpTR = hRespPtEmb->GetBinLowEdge(1), maxMR = hRespMEmb->GetBinLowEdge(nbinsMR+1), maxpTR = hRespPtEmb->GetBinLowEdge(nbinspTR+1);
   
   // histograms to be filled with the Tree   
   TH1F *hpTR = new TH1F(Form("hpTR_%s", suffix.Data()), Form("hpTR_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), nbinspTR, minpTR, maxpTR);
   hpTR->SetMarkerStyle(20);
   //hpTR->SetLineWidth(2);
   hpTR->SetLineColor(colors[0]);
   hpTR->Sumw2();
  
   TH1F *hpTP = new TH1F(Form("hpTP_%s", suffix.Data()), Form("hpTP_%s;#it{p}_{T} (GeV/#it{c})", suffix.Data()), nbinspTP, minpTP, maxpTP);
   hpTP->SetMarkerStyle(20);
   hpTP->SetLineColor(colors[0]);
   hpTP->Sumw2();
   
   TH1F *hMR = new TH1F(Form("hMR_%s", suffix.Data()), Form("hMR_%s;#it{M} (GeV)", suffix.Data()), nbinsMR, minMR, maxMR);
   //hMR->SetLineWidth(2);
   hMR->SetMarkerStyle(20);
   hMR->SetLineColor(colors[0]);
   hMR->Sumw2();
   
   TH1F *hMP = new TH1F(Form("hMP_%s", suffix.Data()), Form("hMP_%s;#it{M} (GeV)", suffix.Data()), nbinsMP, minMP, maxMP);
   //hMP->SetLineWidth(2);
   hMP->SetMarkerStyle(20);
   hMP->SetLineColor(colors[0]);
   hMP->Sumw2();
   
   TLegend *leg = new TLegend(0.4, 0.3, 0.7, 0.6);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   
   for(Int_t ib = 0; ib<entriesIn->GetNbinsX(); ib++){
      
      Int_t nentry = (Int_t)entriesIn->GetBinCenter(ib+1);
      Float_t nTimes = entriesIn->GetBinContent(ib+1);
      if(!fillTreeNtimes) nTimes = 1;
      tree->GetEntry(nentry);
      //old
      //if(vecD->Pt() < embPtCut) continue;
      if(vecP->Pt() < embPtCut) continue;
      
      for(Int_t in = 0; in<nTimes; in++){
      	 hMP-> Fill(vecP->M());
      	 hMR-> Fill(vecD->M());
      	 hpTP->Fill(vecP->Pt());
      	 hpTR->Fill(vecD->Pt());
      }
      
      treeNEntries->Fill();
   }

   //hMP-> Scale(1./hMP ->Integral("width"));
   //hMR-> Scale(1./hMR ->Integral("width"));
   //hpTP->Scale(1./hpTP->Integral("width"));
   //hpTR->Scale(1./hpTR->Integral("width"));
   if(!fillTreeNtimes){
      hMP-> Scale(integMP /hMP ->Integral());
      hMR-> Scale(integMR /hMR ->Integral());
      hpTP->Scale(integPtP/hpTP->Integral());
      hpTR->Scale(integPtR/hpTR->Integral());
   }
   leg->AddEntry(hMP, Form("From Tree %s", fillTreeNtimes ? "(#times N embed)" : "(normalised resp)"), "PL");
   leg->AddEntry(hMEmbPForDiv, "Embedded", "PL");
   leg->AddEntry(hRespMP, "Embedded in Response", "PL");
   leg->AddEntry(hRespMDF, "Bkg+Fluc in Response", "PL");
   
   TCanvas* cEqual = new TCanvas(Form("cEqual%d", pTHB), Form("Are they equal? (pTH bin %d) ", pTHB), 1000, 1000);
   cEqual->Divide(2,2);
   cEqual->cd(1);
   gPad->SetLogy();
   hMP-> Draw();
   hMEmbPForDiv->Draw("sames"); 
   hRespMP->Draw("Psames"); 
   leg->Draw();
   
   
   cEqual->cd(2);
   gPad->SetLogy();
   hMR-> Draw();
   hMEmbRForDiv->Draw("sames");
   hRespMEmb->Draw("Psames");
   hRespMDF->Draw("Psames");
   
   cEqual->cd(3);
   gPad->SetLogy();
   hpTP->Draw();
   hpTEmbPForDiv->Draw("sames");
   hRespPtP->Draw("Psames");
   
   cEqual->cd(4);
   gPad->SetLogy();
   hpTR->Draw();
   hpTEmbRForDiv->Draw("sames");
   hRespPtEmb->Draw("Psames");
   hRespPtDF->Draw("Psames");
   
   
   SaveCv(cEqual);
   
   if(savetree) {
      fout->cd();
      treeNEntries->Write();
   }
}
