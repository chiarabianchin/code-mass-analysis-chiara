void SaveCv(TCanvas* c,TString suffix="", Int_t format=2);   
TList* ReadFile(TString strIn, TString strLst);
TList* ReadFile(TString strIn, TString strDir, TString strLst);
TObject* ReadObjInFile(TString strIn, TString strLst, TString objname);
TH1* CompareDistributions(TVirtualPad* pad, TH1* h1, TH1* h2);
void CalculatePads(Int_t n, Int_t&nx, Int_t&ny, Int_t&dx, Int_t&dy, Int_t perrow = 2, Int_t stdd = 400);
void UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type = "TH1F");
Int_t FindBinInArray(Double_t q, Double_t array[], Int_t dim, Bool_t verbose = kFALSE);
TH1F* ConvertTH1DinF(TH1D *h1d);
TH2F* TransformAxisRanges(TH2F *h2orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew, Int_t ynbinsnew, Float_t yminnew, Float_t ymaxnew);
TH1F* TransformAxisRanges(TH1F *h1orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew);
void DefineNewAxes(TH2 *hOrig, TH2 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Float_t deltaDownY, Float_t deltaUpY, Int_t& newNbinsX, Int_t& newNbinsY, Float_t& minX, Float_t& maxX, Float_t& minY, Float_t& maxY );
void DefineNewAxes(TH1 *hOrig, TH1 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Int_t& newNbinsX, Float_t& minX, Float_t& maxX );

void PrintHisto2(TH2* h1);

