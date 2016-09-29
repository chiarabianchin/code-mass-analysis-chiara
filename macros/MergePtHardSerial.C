void MergePtHardSerial(TString list = "output.list", TString line = "", TString weightList = "JetTagger_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_Jet_AKTChargedR040_PicoTracks_pT4000_E_scheme_TC", Bool_t mergel=kFALSE, TString linel="JetByJetCorrectionOutput") {
   
   gROOT->LoadMacro("$gitJetMass/PtHardUtil/PtHardBinUtilities.cxx+");
   ifstream read(list.Data());
   if (read.is_open())
   {
      while ( getline (read,list.Data()) )
      {
      	 Printf("line: %s",line.Data());
      	 if(weightList.IsNull()) weightList = line;
      	 Printf("weightList: %s",weightList.Data());
      	 
      	 PtHardBinUtilities *ptHardUtil = new PtHardBinUtilities();
      	 
      	 ptHardUtil->SetNameHistosDir("");
      	 if(mergel) {
      	    ptHardUtil->SetNameHistosDir(line.Data());
      	    ptHardUtil->SetNameHistosList(linel.Data()); //for double linked lists
      	 }
      	 else ptHardUtil->SetNameHistosList(line.Data());
      	 ptHardUtil->SetNameWeightList(weightList.Data());
      	 //ptHardUtil->SetNameWeightList(Form("%s/%s",weightList.Data(),weightList.Data())); //for double linked lists
      	 
      	 ptHardUtil->InitFileList(list);
      	 
      	 //ptHardUtil->SetUseDefaultWeights12a15e(kTRUE);
      	 //ptHardUtil->SetUseDefaultWeights13b4(kTRUE);
      	 ptHardUtil->SetMinStatInBin(3,1);
      	 
      	 ptHardUtil->DoWeighting();
      	 
      	 Printf("Weighting done, store output in root file");
      	 
      	 TFile *fout = new TFile(Form("%s.root",line.Data()),"RECREATE");
      	 TList *listOut = ptHardUtil->GetOutputList();
      	 listOut->Print();
      	 fout->WriteTObject(listOut,listOut->GetName());
      	 //listOut->Write();
      	 fout->Write();
      	 fout->Close();
      	 
      	 ptHardUtil->CleanMemory();
      }
   }
}
