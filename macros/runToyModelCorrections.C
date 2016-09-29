void runToyModelCorrections(TString rootFileDetEff = "AnalysisResultsWeighted.root", TString rootFileDeltaPt = "/data/wrk/postdoc/kt/Analysis/Data/grid/LegoMerge/MergeLHC13bLHC13c/AnalysisResultsLHC13bLHC13c.root") {

  gROOT->LoadMacro("/data/wrk/postdoc/kt/Analysis/PlotMacros/ToyModelCorrections.cxx++");

  ToyModelCorrections *toy = new ToyModelCorrections();
  toy->SetFileName(rootFileDetEff);
  toy->LoadFile();
  toy->LoadDiJetResponse();

  toy->SetFileNameDeltaPt(rootFileDeltaPt);

  toy->SetNToyEvents(1e7);

  toy->DoToy();
  toy->WriteResponseToy();

  toy->PlotDiJetResponseKt();

}


//"/data/wrk/postdoc/kt/Analysis/Data/DeltaPtChris/Smear2Output.root"
