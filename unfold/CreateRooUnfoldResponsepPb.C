void CreateRooUnfoldResponse(TString strIn = "B.root", TString strL = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", TString tag = "DerivPart", Int_t binWidthPt = 5., Int_t skipBins = 0, Double_t pt_min = -40., Double_t m_min = -20., Int_t colType = -1) {

  //The response matrix will determine in which binning the data will be evaluated
  //RooUnfold needs to be installed

  //colType
  //-1: use default settings (MB data sets)
  // 1: triggered p-Pb data set ptMinMeas=60
  // 2: triggered p-Pb data set ptMinMeas=40

  //RooUnfold library: set path to RooUnfold as environment variable
  gSystem->Load("$ROOUNFOLD/libRooUnfold.so");

  //strIn: root file containing response matrix
  //strL: name of list from which to take response matrix. Default name is fhnMassResponse*
  //tag: suffix for root file which will be created
  //binWidthPt: default 5 GeV. Can be changed to other value but needs to be possible with binning of original response in THnSparse
  //skipBins: exclude first pt bin of response
  //pt_min: minimum jet pT on detector-level axis
  //m_min: minimum jet mass on detector-level axis

  TFile *f = new TFile(strIn.Data());
  f->ls();
  TList *lst = static_cast<TList*>f->Get(strL.Data());
  lst->Print();

  Int_t iMSub   = 0;
  Int_t iMTrue  = 1;
  Int_t iPtSub  = 2;
  Int_t iPtTrue = 3;

  TH2D *fh2Smear = dynamic_cast<TH2D*>(fhnSparseReduced->Projection(0,2,"E")); //projection on pTrec (x), Mrec (y)
  TH2D *fh2Prior = dynamic_cast<TH2D*>(fhnSparseReduced->Projection(1,3,"E")); //projection on pTpar (x), Mpar (y)

  Int_t nBinPt[2] = {36,30};//meas(5),true
  Double_t ptmin[2] = {pt_min,0.};
  Double_t ptmax[2] = {140.,150.};
  Int_t nBinM[2] = {30,20};//meas(2),true(2)
  Double_t mmin[2] = {m_min,0.};
  Double_t mmax[2] = {40.,40.};

  //change binning for p-Pb triggered data set
  if(colType==1) {
    ptmin[0] = 60.; //meas
    ptmin[1] = 40.; //gen
  }
  else if(colType==2) {
    ptmin[0] = 40.; //meas
    ptmin[1] = 40.; //gen
  }

  Double_t binWidthM = 2.;
  for(Int_t i = 0; i<2; i++) {
    nBinPt[i] = TMath::Nint((ptmax[i]-ptmin[i])/binWidthPt);
    nBinM[i] = TMath::Nint((mmax[i]-mmin[i])/binWidthM);
  }

  if(skipBins>0) {
    nBinPt[1]--;
    ptmin[1]+=binWidthPt;
  }

  //dimensions of measured axis
  TH2D *fh2RespDimM = new TH2D("fh2RespDimM","fh2RespDimM",nBinPt[0],ptmin[0],ptmax[0],nBinM[0],mmin[0],mmax[0]);
  //dimensions of true axis
  TH2D *fh2RespDimT = new TH2D("fh2RespDimT","fh2RespDimT",nBinPt[1],ptmin[1],ptmax[1],nBinM[1],mmin[1],mmax[1]);
  //feed-out of response
  TH2D *fh2Miss     = new TH2D("fh2Miss","fh2Miss",nBinPt[1],ptmin[1],ptmax[1],nBinM[1],mmin[1],mmax[1]);
  cout << "fh2Smear->GetEntries() " << fh2Smear->GetEntries() << endl;
  
  //fill detector-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimM->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimM->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimM->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Smear->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Smear->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimM->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimM->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimM->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Smear->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Smear->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Smear->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      fh2RespDimM->SetBinContent(ix,iy,con);
      fh2RespDimM->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimM");

  //fill particle-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimT->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimT->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimT->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Prior->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Prior->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimT->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimT->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimT->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Prior->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Prior->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Prior->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      fh2RespDimT->SetBinContent(ix,iy,con);
      fh2RespDimT->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimT");

  //response object for RooUnfold
  // the two histograms are like hSmear ( = rec) and hPrior ( = part), but limited at the axis ranges wanted 
  RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp","resp");
  fResponse->Setup(fh2RespDimM,fh2RespDimT);

  //Fill RooUnfoldResponse object
  Int_t* coord = new Int_t[nDim];
  Int_t nbin = fhnSparseReduced->GetNbins();
  for(Int_t bin=0; bin<nbin; bin++) {
    Double_t w = fhnSparseReduced->GetBinContent(bin,coord);
    Double_t msub = fhnSparseReduced->GetAxis(0)->GetBinCenter(coord[0]);
    Double_t mtru = fhnSparseReduced->GetAxis(1)->GetBinCenter(coord[1]); 
    Double_t ptsub = fhnSparseReduced->GetAxis(2)->GetBinCenter(coord[2]); 
    Double_t pttru = fhnSparseReduced->GetAxis(3)->GetBinCenter(coord[3]);
    if(msub>=mmin[0] && msub<mmax[0]
       && mtru>=mmin[1] && mtru<mmax[1]
       && ptsub>=ptmin[0] && ptsub<ptmax[0]
       && pttru>=ptmin[1] && pttru<ptmax[1]
      )
      fResponse->Fill(ptsub,msub,pttru,mtru,w); // filling the response in the requested ranges
    else {
       //outside the ranges fill the missed events in the true matrix
      fResponse->Miss(pttru,mtru,w); 
      fh2Miss->Fill(pttru,mtru,w);
    }
  }
  
  delete [] coord;

  //Write response + 2D histos to file
  TFile *fout = new TFile(Form("response%s.root",tag.Data()),"RECREATE");
  hn->Write("fhn");
  fhnSparseReduced->Write("fhnReduced");
  fResponse->Write("resp");
  fh2Smear->Write("fh2Smear");
  fh2Prior->Write("fh2Prior");
  fh2RespDimM->Write();
  fh2RespDimT->Write();
  fh2Miss->Write();
  fout->Write();
  fout->Close();
}
