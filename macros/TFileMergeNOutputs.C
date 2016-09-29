#include <TFileMerger.h>
#include <TFile.h>
#include <TString.h>
#include <fstream>

void TFileMergeNOutputs(Int_t nfiles, TString *inputs, TString outname);


void TFileMergeNOutputs(Int_t nfiles, const char* inputfilelist, TString outname){
	ifstream inputfile(inputfilelist);
	TString inputs[nfiles];
	string bypass;
	if(!inputfile.is_open()){
		Printf("Error");
		return;
	}
	for(Int_t in = 0; in<nfiles;in++) {
		Printf("%d", in);
		getline(inputfile, bypass);
		inputs[in]=bypass;
		in++;
	}
	inputfile.close();
	TFileMergeNOutputs(nfiles, inputs, outname);
	
}

//
void TFileMergeNOutputs(Int_t nfiles, TString *inputs, TString outname){
	TFileMerger merger;
	merger.OutputFile(outname);
	for(Int_t i = 0; i<nfiles; i++){
		merger.AddFile(inputs[i]);
	}
	
	merger.Merge();

}
