#include "/data/Work/MyCodeJetMass/macros/DrawTreesInOutputTask.C"
void RunDrawEntriesUsed(TString fileEmb, TString listEmb, TString listResp, TString fileTreeIn, Double_t embPtCut = 0, Bool_t savetree = kFALSE, Int_t pTHB = 1){
   DrawEntriesUsed(fileEmb, listEmb, listResp, fileTreeIn, embPtCut, savetree, pTHB );
}
