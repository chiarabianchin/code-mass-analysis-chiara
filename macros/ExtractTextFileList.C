#include <Riostream.h>
#include <TList.h>
#include <TFile.h>
#include <TString.h>

void ExtractTextFileList(TString sampleTFile = "01/AnalysisResults.root", TString outputTextFile = "listOfListNamesToMerge.txt", TString outputTextFile2 = "listOfFilesToMerge.txt"){
   
   TFile *fin = new TFile(sampleTFile);
   if(!fin->IsOpen()){
      Printf("File not found");
      return;
   }
   
   TList *listkeys = fin->GetListOfKeys();
   Int_t nkeys = listkeys->GetEntries();
   
   ofstream output;
   output.open(outputTextFile.Data());
   ofstream output2;
   output2.open(outputTextFile2.Data());
   
   for(Int_t i = 0; i < nkeys; i++){
      TString name = listkeys->At(i)->GetName();
      Printf("%s", name.Data());
      output<<name.Data()<<"\n";
      output2<<name.Data()<<".root\n";
   }
   output.close();
   output2.close();
   Printf("File %s written", outputTextFile.Data());
   Printf("File %s written", outputTextFile2.Data());
}

