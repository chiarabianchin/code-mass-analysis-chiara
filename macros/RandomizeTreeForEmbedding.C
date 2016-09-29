#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TRandom3.h>

#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

TTree* ReadTree (TString strIn, TString nameTree);

// random generator function:
int myrandom (int i) { return std::rand()%i;}

// definition
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

void RandomizeTreeForEmbedding(TString filename = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1064/output/AnalysisResults.root", TString treename = "fTreeJet", Int_t version = 1, Int_t pthardB = 1, Double_t minpTDetCut = -1, Double_t minpTParCut = 15, Int_t maxentriesbin = 150000, TString treeFileName = "outR"){
   // input file history 
   // "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1054/outputpTHardBins/mergeRuns/plainMerge/out.root"
   //version
   // 0 = 2 branches "fJetDet.", "fJetPart."
   // 1 = 4 branches "fJetDet.", "fJetPart.", "fXsection", "fNTrials"
   // 2 = 5 branches "fJetDet.", "fJetPart.", "fXsection", "fNTrials" "fNCountsPthardB"
   
   TStopwatch watch;
   watch.Start();
   if(minpTDetCut > -1) Printf("WARNING!!!!!!!!! A pt cut at detector level will be applied, are you sure it's what you want?");
   Printf("pT part cut = %f", minpTParCut);

   TTree* tree = ReadTree(filename, treename);
   if(!tree){
      Printf("Tree not found, exit");
      return;
   }
   tree->LoadBaskets(13ul * 1024 * 1024 * 1024);
   
   //define vector of ttree entries
   std::vector<int> myvector;
   std::srand ( unsigned ( std::time(0) ) );
   
   const Int_t nEntries = tree->GetEntries();
   // set entry values
   for (int i=0; i < nEntries; i++)
      myvector.push_back(i);

   Printf("Push_back vector done");
   // random swap
   std::random_shuffle ( myvector.begin(), myvector.end(), myrandom);
   
   TLorentzVector *vecD = 0x0;
   TLorentzVector *vecP = 0x0;
   Float_t xsec = 0; 
   Float_t ntrials = 0;
   Int_t pthardbin = 0;
   
   Int_t nbranches = 2;
   if(version == 1) nbranches = 4;
   if(version == 2) nbranches = 5;
   const Int_t nbranchesc = nbranches; 
   TString branches[nbranchesc];
   if(version == 0 ){
      branches[0] = "fJetDet.";
      branches[1] = "fJetPart.";
   }
   if(version >= 1 ){
      branches[0] = "fJetDet.";
      branches[1] = "fJetPart.";
      branches[2] = "fXsection";
      branches[3] = "fNTrials";
     
   }
   if(version == 2)  branches[4] = "fPthardB";
   
   tree->SetBranchAddress(branches[0], &vecD);
   tree->SetBranchAddress(branches[1], &vecP);
   if(version == 1) {
      tree->SetBranchAddress(branches[2], &xsec);
      tree->SetBranchAddress(branches[3], &ntrials);
   }
   if(version == 2)  tree->SetBranchAddress(branches[4], &pthardbin);
   
   TFile *fOut = TFile::Open(Form("%s%d.root", treeFileName.Data(), pthardB), "recreate");

   TTree *outTree = new TTree("fTreeJetR", "fTreeJetR");
   Float_t fNCountsPthardB;
   outTree->Branch(branches[0], &vecD);
   outTree->Branch(branches[1], &vecP);
   if(version >= 1) {
      outTree->Branch(branches[2], &xsec);
      outTree->Branch(branches[3], &ntrials);
      
   }
   if(version == 2){
      outTree->Branch("fNCountsPthardB", &fNCountsPthardB); // NOTE: if this name is changed it affects the settings of the task AliAnalysisTaskJetShapeBase!!!
   }
   //histogram used to fill the new branch
   TH1F *hNpthB= new TH1F("hNpthB", "Entries per pthard bin; pT-hard Bin", 10, 0.5, 10.5 );
   if(version == 2){   
      for(Int_t i = 0; i < nEntries; i++){
      	 Int_t entry = tree->GetEntry(i);
      	 hNpthB->Fill(pthardbin);
      }
   }
   
   Int_t numberfilled = 0;
   for(Int_t i = 0; i < nEntries; i++){
      //Printf("Entry %d", i);
      
      if (i % 100000 == 0) Printf("Entry %d", i);
      Int_t entry = tree->GetEntry(myvector[i]);
      if(version == 2) fNCountsPthardB = hNpthB->GetBinContent(hNpthB->FindBin(pthardbin));
      if(entry) {
      	 if(minpTDetCut > -1 && vecD->Pt()< minpTDetCut) continue;
      	 if(minpTParCut > -1 && vecP->Pt()< minpTParCut) continue;
      	 if(numberfilled < maxentriesbin) {
      	    outTree->Fill();
      	    numberfilled++;
      	 }
      	 else break;
      }
   }
   Printf("Done, Saving now...");
   fOut->WriteTObject(outTree);
   fOut->Close();
   Printf("Done.");
   
   watch.Stop();
   watch.Print();
}

