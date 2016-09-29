#include <TF1.h>
#include <TList.h>
#include <THnSparse.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>

TH1F* ConvolutionH1XH2(TH1F *h1, TH1F *h2, Bool_t drawlog = kFALSE);
TH1F* Convolution1DRandom(TH1F *h1, TH1F *h2, Int_t ndraw, Bool_t posx = kTRUE);
TH2F* Convolution2DRandom(TH2F *h1, TH2F *h2, Int_t ndraw, Bool_t posx = kTRUE, Bool_t posy = kTRUE); //true with constituent method, fluctuations only can't go negative. mass true also with detector I think 
void SetToSameBins(TH2* h1, TH2* h2, TH2*& h1bis, TH2*& h2bis);
TH2F* Convolution2D(TH2F *h1, TH2F *h2);

TH2D* SmearOfdxdy(TH2D* h2, Double_t dx, Double_t dy, Double_t f, Int_t uid);
TH2D* LoopOnFluctuations(TH2D *hdMdpT, TH2D *hMpTpart, Int_t iid = 0);

void PrintTHnSparse(THnSparse *hsp, TString appendString = "");
TH2D* MakeGaussiandMdpT(Int_t mode = 1);
void PrintCompareTH2(TH2* h);

//________________________________________________________________________________________________

TH1F* ConvolutionH1XH2(TH1F *h1, TH1F *h2, Bool_t drawlog){
   
   if(!h1 || !h2){
      Printf("Need valid input!");
      return 0x0;
   }
   
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h2->GetNbinsX();
   Double_t binW1 = h1->GetBinWidth(3), binW2 = h2->GetBinWidth(3);
   Double_t low1 = h1->GetBinLowEdge(1), low2 = h2->GetBinLowEdge(1);
   Double_t up1 = h1->GetBinLowEdge(nBins1+1), up2 = h2->GetBinLowEdge(nBins2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBins2, binW2, low2, up2);
   
   Double_t commonlow = TMath::Min(low1, low2), commonup = TMath::Max(up1, up2);
   Double_t commonW = TMath::Max(binW1, binW2);
   Int_t commonNBins = (commonup - commonlow)/commonW;
   Printf("Convol) N bins = %d, Width = %f, Range = %f - %f", commonNBins, commonW, commonlow, commonup);
   Double_t minWidth = TMath::Min(binW1, binW2);
   
   Int_t rebin = (Int_t)commonW/minWidth;
   Printf("%f/%f = %f or %d", commonW, minWidth, commonW/minWidth, rebin);
   
   if(minWidth == binW1) h1->Rebin(rebin);
   else h2->Rebin(rebin);
   
   Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsX(), h2->GetName(), h2->GetNbinsX());
   TH1F *h1bis = new TH1F(Form("%sb", h1->GetName()), Form("%s bis", h1->GetTitle()), commonNBins, commonlow, commonup);
   TH1F *h2bis = new TH1F(Form("%sb", h2->GetName()), Form("%s bis", h2->GetTitle()), commonNBins, commonlow, commonup);
   //refill the original histograms using the same binning
   for(Int_t i = 0; i < commonNBins; i++){
      Int_t fillingbin = h1->FindBin(h1bis ->GetBinCenter(i+1));
      Double_t fillWith = 0, error = 0;
      if(fillingbin>0 && fillingbin <= nBins1) {
      	 fillWith = h1->GetBinContent(fillingbin);
      	 error = h1->GetBinError(fillingbin);
      }
      h1bis ->SetBinContent(i+1, fillWith);
      h1bis ->SetBinError(i+1, error);
      Printf("Filling new bin %d with old bin %d content %f", i+1, fillingbin, fillWith);
      fillingbin = h2->FindBin(h2bis ->GetBinCenter(i+1));
      fillWith = 0; error = 0;
      if(fillingbin>0 && fillingbin <= nBins2) {
      	 fillWith = h2->GetBinContent(fillingbin);
      	 error = h2->GetBinError(fillingbin);
      }
      h2bis ->SetBinContent(i+1, fillWith);
      h2bis ->SetBinError(i+1, error);
   }
   
   h2bis->Scale(1./h2bis->Integral());
      
   TH1F *h1Modif[commonNBins];
   TH1F *h1ModifTot=0x0;
   TCanvas *cConvol = new TCanvas(Form("cConvol%s", h2->GetName()), Form("Convolute with %s", h2->GetName()));
   cConvol->cd();
   h2->Draw("hist");
   h2bis->Draw("histsames");
   TCanvas *cModifbins = new TCanvas(Form("cModifbins%s", h1->GetName()), Form("Modified %s", h1->GetName()));
   TLegend *leg = new TLegend(0.5,0.65, 0.9, 0.8);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   
   
   //Printf("Start loop on %d x %d bins", nBins2/DeltaBin, nBins1);

   for(Int_t jb2 = 0; jb2 < commonNBins ; jb2++ ){
      //this gives the probability of being in the jb2 bin (can it be also GetBinContent?)
      Double_t fract = h2bis->Integral(jb2, jb2, "width");
      Printf("Multiply yield by %f", fract);
      Double_t xValue = h2bis->GetBinCenter(jb2);
      h1Modif[jb2] = new TH1F(Form("%sModif-Bin2%d", h1->GetName(), jb2), "; Mass; Entries", commonNBins, commonlow, commonup);
      h1Modif[jb2]->SetLineColor(1);//(colors[jb2]);
      h1Modif[jb2]->SetMarkerColor(1);//(colors[jb2]);
      
      for(Int_t ib1 = 0; ib1<commonNBins; ib1++){
      	 //smear the X of the histogram to be changed with the corresponding X of the smearing histogram
      	 Double_t newX = h1bis->GetBinCenter(ib1+1) + xValue;
      	 Int_t newB = h1bis->FindBin(newX);
      	 //Fill the newX with the original content multiplied by the probability of the smearing
      	 Double_t fillWith =  h1bis->GetBinContent(ib1+1)  * fract;
      	 h1Modif[jb2]->SetBinContent(newB,fillWith);
      	 h1Modif[jb2]->SetBinError(newB, h1bis->GetBinError(ib1+1) * fract);
      	 
      }
      //sum up all the new X smeared
      if(!h1ModifTot) {
      	 h1ModifTot = (TH1F*)h1Modif[jb2]->Clone(Form("%sModif", h1->GetName()));
      	 h1ModifTot->SetLineColor(colors[2]);
      	 h1ModifTot->SetMarkerStyle(24);
      	 h1ModifTot->SetMarkerColor(colors[2]);
      }
      else h1ModifTot->Add(h1Modif[jb2]);
      
      cModifbins->cd();
      if(jb2 == 0) h1Modif[jb2]->GetYaxis()->SetRangeUser(0,1e-2);
      h1Modif[jb2]->Draw("sames");
      
   }
   
   leg->AddEntry(h1, Form("Original #mu = %.2f, RMS = %.2f", h1->GetMean(), h1->GetRMS()), "P");
   leg->AddEntry(h1ModifTot, Form("Convoluted #mu = %.2f, RMS = %.2f", h1ModifTot->GetMean(), h1ModifTot->GetRMS()), "P");
   TCanvas *cFinal = new TCanvas(Form("cFinal%s", h1->GetName()), Form("Final %s", h1->GetName()));
   cFinal->cd();
   if(drawlog) gPad->SetLogy();
   //h1bis->Draw("sames");
   h1->Draw("sames");
   h1ModifTot->Draw("sames");
   leg->Draw();
   
   //SaveCv();
   SaveCv(cConvol);
   SaveCv(cFinal);
   return h1ModifTot;
}

//________________________________________________________________________________________________
TH1F* Convolution1DRandom(TH1F *h1, TH1F *h2, Int_t ndraw, Bool_t posx){
   if(!h1 || !h2){
      Printf("Need valid input!");
      return 0x0;
   }
   
   Int_t nbins = h1->GetNbinsX();
   TH1F *hNew = new TH1F(Form("%sModif", h1->GetName()), Form("%s;%s;%s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), nbins , h1->GetBinLowEdge(1), h1->GetBinLowEdge(nbins+1));
   hNew->Sumw2();
   for(Int_t ibin = 0 ; ibin<nbins; ibin++){
      Float_t xorig = h1->GetBinCenter(ibin);
      Float_t content = h1->GetBinContent(ibin);
      //Printf("Debug -> bin %d center %.2f --->", ibin, xorig);
      for(Int_t irnd = 0; irnd < ndraw; irnd++){
      	 
      	 Double_t xnew = xorig + h2->GetRandom();
      	 if(posx && xnew < 0) continue;
      	 //Printf("%.2f", xnew);
      	 hNew->Fill(xnew, content/(Float_t)ndraw);
      }
   }
   return hNew;
}

//________________________________________________________________________________________________
TH2F* Convolution2DRandom(TH2F *h1, TH2F *h2, Int_t ndraw, Bool_t posx, Bool_t posy){
   if(!h1 || !h2){
      Printf("Need valid input!");
      return 0x0;
   }
   
   Int_t nbinsX = h1->GetNbinsX();
   Int_t nbinsY = h1->GetNbinsY();
   TH2F *hNew = new TH2F(Form("%sNew", h1->GetName()), Form("%s;%s;%s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), nbinsX , h1->GetXaxis()->GetBinLowEdge(1), h1->GetXaxis()->GetBinLowEdge(nbinsX+1), nbinsY , h1->GetYaxis()->GetBinLowEdge(1), h1->GetYaxis()->GetBinLowEdge(nbinsY+1));
   hNew->Sumw2();
   
   //check the extracted numbers on X axis
   TCanvas *cdebug = new TCanvas(Form("cDebug%s", h1->GetName()), Form("cDebug%s", h1->GetName()));
   TH1F* hDebugRndm = new TH1F(Form("hDebugRndm%s", h1->GetName()), Form("hDebugRndm%s;#delta #it{M} (GeV/#it{c})", h1->GetName()), h2->GetNbinsX(), h2->GetXaxis()->GetBinLowEdge(1), h2->GetXaxis()->GetBinLowEdge(h2->GetNbinsX()+1));
   hDebugRndm->SetMarkerStyle(29);
   hDebugRndm->Sumw2();
   Bool_t debug=kTRUE;
   TLegend *leg = new TLegend(0.5, 0.4, 0.9, 0.6, "Must be compatible");
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hDebugRndm, "Extracted", "PL");
   
   for(Int_t ibinX = 0 ; ibinX<nbinsX; ibinX++){
      for(Int_t ibinY = 0 ; ibinY<nbinsY; ibinY++){
      	 Int_t globorig = h1->GetBin(ibinX, ibinY);
      	 Float_t xorig = h1->GetXaxis()->GetBinCenter(ibinX);
      	 Float_t yorig = h1->GetYaxis()->GetBinCenter(ibinY);
      	 Float_t content = h1->GetBinContent(ibinX, ibinY);
      	 if(content > 1e-6){
      	    //Printf("Debug -> bin %d center %.2f, %2f, content %e --->", globorig, xorig, yorig, content);
      	    for(Int_t irnd = 0; irnd < ndraw; irnd++){
      	       Double_t rndX, rndY;
      	       h2->GetRandom2(rndX, rndY);
      	       if(debug) {
      	       	  hDebugRndm->Fill(rndX);
      	       	  
      	       }
      	       Double_t xnew = xorig + rndX;
      	       Double_t ynew = yorig + rndY;
      	       
      	       // the following request makes the folded distribution comply with the definition of mass, M > 0 (and also pT here, need to check) 
      	       
      	       if(posx && xnew < 0 ) continue;
      	       if(posy && ynew < 0 ) continue;
      	       //Printf("Fill %.2f, %.2f with %e", xnew, ynew, content/(Float_t)ndraw);
      	       hNew->Fill(xnew, ynew, content/(Float_t)ndraw);
      	       //hNew->Fill(xnew, ynew, content);
      	    }
      	    debug=kFALSE;
      	 }
      }
   }
   cdebug->cd();
   hDebugRndm->GetXaxis()->SetRangeUser(-5, 20);
   //hDebugRndm->Scale(1./hDebugRndm->Integral());
   hDebugRndm->Draw("E");
   
   TH1D *htmp = h2->ProjectionX(Form("hdMproj%s", h1->GetName()));
   htmp->Sumw2();
   htmp->SetMarkerStyle(30);
   htmp->Scale(ndraw/htmp->Integral());
   leg->AddEntry(htmp, "Projection 2D input", "PL");
   htmp->Draw("Esames");
   leg->Draw();
   return hNew;
}

//________________________________________________________________________________________________

TH2D* Folding2D(THnSparse* hResp/*4D response matrix*/, TH2D* hTrue, Int_t axMm, Int_t axpTm, Int_t axMt, Int_t axpTt){
   if(!hResp || !hTrue){
      Printf("Need valid input!");
      return 0x0;
   }
   Int_t iax = axMm; //mass meas
   Int_t jax = axpTm; //pt meas
   Int_t kax = axMt; //mass true
   Int_t lax = axpTt; //pT true
   Int_t nbinsi = hResp->GetAxis(iax)->GetNbins();
   Int_t nbinsj = hResp->GetAxis(jax)->GetNbins();
   
   Int_t nbinsk = hResp->GetAxis(kax)->GetNbins();
   Int_t nbinsl = hResp->GetAxis(lax)->GetNbins();
   
   TH2D *hFold = (TH2D*)hTrue->Clone("hFold");
   
   for(int i=0;i<nbinsi;i++){
      for(int j=0;j<nbinsj;j++){
      	 double effects=0;
      	 for(int k=0;k<nbinsk;k++){
      	    for(int l=0;l<nbinsl;l++){
      	       Int_t indexes[4] = {i+1,j+1,k+1,l+1};
      	       Int_t globIndx = hResp->GetBin(indexes);
       	       effects=effects+hTrue->GetBinContent(k+1,l+1)*hResp->GetBinContent(globIndx);
      	       
       	    }
      	 }
      	 hFold->SetBinContent(i+1,j+1,effects);
      }
   }
   return hFold;
}

//________________________________________________________________________________________________

TH2F* Convolution2D(TH2F *h1, TH2F *h2){
   
   if(!h1 || !h2){
      Printf("Need valid input!");
      return 0x0;
   }
   
   TH2 *h1bis = 0x0, *h2bis = 0x0;
   
   SetToSameBins(h1, h2, h1bis, h2bis);
   
   Int_t commonNBinsi = h1bis->GetNbinsX(); //mass
   Int_t commonNBinsj = h1bis->GetNbinsY(); //pT
   Int_t commonlowi   = h1bis->GetXaxis()->GetBinLowEdge(1), commonupi   = h1bis->GetXaxis()->GetBinLowEdge(commonNBinsi+1); //mass
   Int_t commonlowj   = h1bis->GetXaxis()->GetBinLowEdge(1), commonupj   = h1bis->GetXaxis()->GetBinLowEdge(commonNBinsj+1); //pT
  
   
   //h2bis->Scale(1./h2bis->Integral());
   h2->Scale(1./h2->Integral());
   
   TH2F *h1Modif[commonNBinsi*commonNBinsj];
   TH2F *h1ModifTot=0x0;
   TCanvas *cModifbins = new TCanvas(Form("cModifbins%s", h1->GetName()), Form("Modified %s", h1->GetName()));
   for(Int_t ib2 = 0; ib2 < commonNBinsi ; ib2++ ){
      for(Int_t jb2 = 0; jb2 < commonNBinsj ; jb2++ ){
//for(Int_t ib2 = 0; ib2 < 1 ; ib2++ ){
      //for(Int_t jb2 = 0; jb2 < 1 ; jb2++ ){
      	 //this gives the probability of being in the jb2 bin (can it be also GetBinContent?)
      	 Double_t fract = h2bis->GetBinContent(ib2, jb2);
      	 Printf("Multiply yield by %f", fract);
      	 Double_t xValue = h2bis->GetXaxis()->GetBinCenter(ib2);
      	 Double_t yValue = h2bis->GetYaxis()->GetBinCenter(jb2);
      	 Printf("xValue = %f, yValue = %f", xValue, yValue);
      	 Int_t glb2 = h2bis->GetBin(ib2, jb2);
      	 h1Modif[glb2] = new TH2F(Form("%sModif-Bin2%d", h1->GetName(), glb2), "; Mass; pT; Entries", commonNBinsi, commonlowi, commonupi, commonNBinsj, commonlowj, commonupj);
      	 h1Modif[glb2]->SetLineColor(1);//(colors[jb2]);
      	 h1Modif[glb2]->SetMarkerColor(1);//(colors[jb2]);
      	 
      	 for(Int_t ib1 = 0; ib1<commonNBinsi; ib1++){
      	    for(Int_t jb1 = 0; jb1<commonNBinsj; jb1++){
      	       //smear the X of the histogram to be changed with the corresponding X of the smearing histogram
      	       Double_t newX = h1bis->GetXaxis()->GetBinCenter(ib1+1) + xValue;
      	       Double_t newY = h1bis->GetYaxis()->GetBinCenter(jb1+1) + yValue;
      	       Printf("Bin x, y  = %f, %f", newX, newY);
      	       Int_t newB = h1bis->FindBin(newX, newY);
      	       //Fill the newX with the original content multiplied by the probability of the smearing
      	       Double_t fillWith =  h1bis->GetBinContent(ib1+1,jb1+1)  * fract;
      	       h1Modif[jb2]->SetBinContent(newB,fillWith);
      	       h1Modif[jb2]->SetBinError(newB, h1bis->GetBinError(ib1+1,jb1+1) * fract);
      	       
      	    }
      	 }
      	 
      	 //sum up all the new X smeared
      	 if(!h1ModifTot) {
      	    h1ModifTot = (TH2F*)h1Modif[glb2]->Clone(Form("%sModif", h1->GetName()));
      	    h1ModifTot->SetLineColor(colors[2]);
      	    h1ModifTot->SetMarkerStyle(24);
      	    h1ModifTot->SetMarkerColor(colors[2]);
      	 }
      	 else h1ModifTot->Add(h1Modif[glb2]);
      	 
      	 cModifbins->cd();
      	 if(jb2 == 0) h1Modif[glb2]->GetYaxis()->SetRangeUser(0,1e-2);
      	 h1Modif[glb2]->Draw("sames");
      	 
      }
   }
   return h1ModifTot;
}

//________________________________________________________________________________________________

void SetToSameBins(TH2* h1, TH2* h2, TH2*& h1bis, TH2*& h2bis){
   
   //check compatibility X axis
   Int_t nBinsX1 = h1->GetNbinsX(), nBinsX2 = h2->GetNbinsX();
   Double_t binWX1 = h1->GetXaxis()->GetBinWidth(3), binWX2 = h2->GetXaxis()->GetBinWidth(3);
   Double_t lowX1 = h1->GetXaxis()->GetBinLowEdge(1), lowX2 = h2->GetXaxis()->GetBinLowEdge(1);
   Double_t upX1 = h1->GetXaxis()->GetBinLowEdge(nBinsX1+1), upX2 = h2->GetXaxis()->GetBinLowEdge(nBinsX2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBinsX1, binWX1, lowX1, upX1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBinsX2, binWX2, lowX2, upX2);
   
   Double_t commonlowX = TMath::Min(lowX1, lowX2), commonupX = TMath::Max(upX1, upX2);
   Double_t commonWX = TMath::Max(binWX1, binWX2);
   Int_t commonNBinsX = (commonupX - commonlowX)/commonWX;
   Printf("Convol) N bins = %d, Width = %f, Range = %f - %f", commonNBinsX, commonWX, commonlowX, commonupX);
   Double_t minWidthX = TMath::Min(binWX1, binWX2);
   
   Int_t rebinX = (Int_t)commonWX/minWidthX;
   Printf("%f/%f = %f or %d", commonWX, minWidthX, commonWX/minWidthX, rebinX);
   
   if(minWidthX == binWX1) h1->RebinX(rebinX);
   else h2->RebinX(rebinX);
   
   Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsX(), h2->GetName(), h2->GetNbinsX());

   //check compatibility Y axis
   Int_t nBinsY1 = h1->GetNbinsY(), nBinsY2 = h2->GetNbinsY();
   Double_t binWY1 = h1->GetYaxis()->GetBinWidth(3), binWY2 = h2->GetYaxis()->GetBinWidth(3);
   Double_t lowY1 = h1->GetYaxis()->GetBinLowEdge(1), lowY2 = h2->GetYaxis()->GetBinLowEdge(1);
   Double_t upY1 = h1->GetYaxis()->GetBinLowEdge(nBinsY1+1), upY2 = h2->GetYaxis()->GetBinLowEdge(nBinsY2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBinsY1, binWY1, lowY1, upY1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBinsY2, binWY2, lowY2, upY2);
   
   Double_t commonlowY = TMath::Min(lowY1, lowY2), commonupY = TMath::Max(upY1, upY2);
   Double_t commonWY = TMath::Max(binWY1, binWY2);
   Int_t commonNBinsY = (commonupY - commonlowY)/commonWY;
   Printf("Convol) N bins = %d, Width = %f, Range = %f - %f", commonNBinsY, commonWY, commonlowY, commonupY);
   Double_t minWidthY = TMath::Min(binWY1, binWY2);
   
   Int_t rebinY = (Int_t)commonWY/minWidthY;
   Printf("%f/%f = %f or %d", commonWY, minWidthY, commonWY/minWidthY, rebinY);
   
   if(minWidthY == binWY1) h1->RebinY(rebinY);
   else h2->RebinY(rebinY);
   
   Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsY(), h2->GetName(), h2->GetNbinsY());
   
   h1bis = new TH2F(Form("%sb", h1->GetName()), Form("%s bis", h1->GetTitle()), commonNBinsX, commonlowX, commonupX, commonNBinsY, commonlowY, commonupY);
   h2bis = new TH2F(Form("%sb", h2->GetName()), Form("%s bis", h2->GetTitle()), commonNBinsX, commonlowX, commonupX, commonNBinsY, commonlowY, commonupY);
   
   return;
}

//_________________________________________________________________________________________________

TH2D* SmearOfdxdy(TH2D* h2, Double_t dx, Double_t dy, Double_t f, Int_t uid){
   //TH2D *h2out = (TH2D*)h2->Clone(Form("%s%d", h2->GetName(), uid));
   
   Double_t minx = -5,  maxx = h2->GetXaxis()->GetBinLowEdge(h2->GetNbinsX()+1);
   Double_t miny = -30, maxy = h2->GetYaxis()->GetBinLowEdge(h2->GetNbinsY()+1);
   
   Double_t binwx = h2->GetXaxis()->GetBinWidth(0), binwy = h2->GetYaxis()->GetBinWidth(0);
   
   //assuming that the original histogram starts from 0
   Int_t  nbinsx = TMath::Abs(minx) / binwx + h2->GetNbinsX();
   Int_t  nbinsy = TMath::Abs(miny) / binwy + h2->GetNbinsY();
   //Printf("N bins orig %f * new bins ", h2->GetYaxis()->GetBinWidth(0));
   TH2D *h2out = new TH2D(Form("%s%d", h2->GetName(), uid), Form("Smeared number %d;#it{M}_{smear} (GeV); #it{p}_{T, smear} (GeV/c)", uid), nbinsx, minx, maxx, nbinsy, miny, maxy);
   h2out->Sumw2();
   
   
   for(Int_t i = 0; i < h2->GetNbinsX(); i++){
   //for(Int_t i = 0; i < 10; i++){
      Double_t xvalue = h2->GetXaxis()->GetBinCenter(i+1);
      Double_t newXval = xvalue + dx;
      //Printf("DX is %e (%.3f bins), DY is %e (%.3f bins)", dx, dx/binwx, dy, dy/binwy);
      Int_t binx = h2out->GetXaxis()->FindBin(newXval);
      for(Int_t j = 0; j < h2->GetNbinsY(); j++){
      	 //for(Int_t j = 0; j < 16; j++){
      	 Double_t yvalue = h2->GetYaxis()->GetBinCenter(j+1);
      	 Double_t newYval = yvalue + dy;
      	 Int_t biny = h2out->GetYaxis()->FindBin(newYval);
      	 Double_t cont = h2->GetBinContent(i+1, j+1);
      	 Double_t err  = h2->GetBinError(i+1, j+1);
      	 //Printf("Bin error = %e, rel err = %e", err, err/cont);
      	 Double_t newCont = cont * f ;
      	 Double_t newContE= err  * f ;
      	 //Printf("NEW Bin error = %e, rel err = %e", newContE, newContE/newCont);
      	 //Printf("Old Bin (%d, %d) moves to Bin (%d, %d) -> (%.2f, %.2f), counts reset to %e * %e %.2e (rel er %.1f)", i+1, j+1, binx, biny, newXval, newYval, cont, f, newCont, newContE/newCont);
      	 h2out->SetBinContent(binx, biny, newCont);
      	 //h2out->SetBinError(binx, biny, newContE);
      } //loop on y bins
   
   } //loop on x bins
   
   return h2out;
}

//_________________________________________________________________________________________________

TH2D* LoopOnFluctuations(TH2D *hdMdpT, TH2D *hMpTpart, Int_t iid){
   
   TH2D *h2dSmeared = 0x0; //contains the final output
   
   if(iid != 0) hMpTpart->SetName(Form("%sBin%d", hMpTpart->GetName(), iid));
   
   PrintCompareTH2(hdMdpT); 
   PrintCompareTH2(hMpTpart);
   //histogram for debugging
   TH1D* hdMinput = new TH1D(Form("hdMinput_B%d", iid), "#delta M actually used; #delta M; Number of iterations", hdMdpT->GetNbinsY(), hdMdpT->GetYaxis()->GetBinLowEdge(1), hdMdpT->GetYaxis()->GetBinLowEdge(hdMdpT->GetNbinsY()));
   TH1D* hdpTinput = new TH1D(Form("hdpTinput_B%d", iid), "#delta #it{p}_{T} actually used; #delta #it{p}_{T}; Number of iterations", hdMdpT->GetNbinsX(), hdMdpT->GetXaxis()->GetBinLowEdge(1), hdMdpT->GetXaxis()->GetBinLowEdge(hdMdpT->GetNbinsX()));
   
   
   Int_t id = 0; //this give a unique identification to each bin of the loop
   Bool_t stoploop = kFALSE;
   for(Int_t ibin = 0; ibin < hdMdpT->GetNbinsX(); ibin++){ // loop on dMdpT X bins (dpT)
   //for(Int_t ibin = 27; ibin < 28; ibin++){
      Double_t dx = hdMdpT->GetXaxis()->GetBinCenter(ibin+1);
      
      for(Int_t jbin = 0; jbin < hdMdpT->GetNbinsY(); jbin++){ // loop on dMdpT Y bins (dM)
      	 //if(hdMdpT->GetBinContent(ibin+1, jbin+1) < 1.e-15) continue;
      	 //if(stoploop) break;
      	 Double_t dy = hdMdpT->GetYaxis()->GetBinCenter(jbin+1);
      	 Double_t f  = hdMdpT->GetBinContent(ibin+1, jbin+1);
      	 hdMinput->Fill(dy);
      	 hdpTinput->Fill(dx);
      	 //if(jbin == 24) {
      	    Printf("Bin (%d, %d) -> (%.2f, %.2f), counts in bin %.2e", ibin+1, jbin+1, dx, dy, f);
      	 //}
      	 
      	 TH2D* h2id = SmearOfdxdy(hMpTpart, dy, dx, f, id); //note that x and y are inverted wrt hdMdpT
      	 if(!h2dSmeared) {
      	    h2dSmeared = (TH2D*)h2id->Clone(Form("h2dSmeared_B%d", iid));
      	    h2dSmeared->SetTitle("Smeared Total");
      	    Printf("Smeared histogram bins pT %d, %f", h2dSmeared->GetNbinsY(), h2dSmeared->GetYaxis()->GetBinWidth(0));
      	    h2dSmeared->GetYaxis()->SetTitleOffset(1.2);
      	 }
      	 else h2dSmeared->Add(h2id);
      	 
      	 id++;
      	 
      }// end loop on dMdpT Y bins (dM)
      //Printf("\n");
      //if(id>2) stoploop = kTRUE;
   }// end loop on dMdpT X bins (dpT)
   
   if(id == 0) Printf("Returning empty histogram, sorry");
   
   TCanvas *cdM = new TCanvas("cdM", "Delta M input", 600, 600);
   cdM->cd();
   hdMinput->Draw();
   TCanvas *cdpt = new TCanvas("cdpt", "Delta pt input", 600, 600);
   cdpt->cd();
   hdpTinput->Draw();
   SaveCv(cdM);
   SaveCv(cdpt);

   return h2dSmeared;
}
//_________________________________________________________________________________________________

void GetAndLoop(Bool_t onedMdpt = kTRUE, Bool_t brecoLev = kTRUE, TString file1 = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train976/analysis/UnfoldingThnSparse.root", TString file2 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/UnfoldingMatrix.root"){
   
   TFile *fin1 = new TFile(file1);
   if(!fin1->IsOpen()){
      Printf("%s not found", file1.Data());
      return;
   
   }
   TH2D *hdMdpT = 0x0;
   THnSparseF *hdMdpTMparpTpar = (THnSparseF*)fin1->Get("fhnDeltaMass_0_proj_0_1_3_5");
   if(onedMdpt || !hdMdpTMparpTpar){
      
      Printf("No Thnsparse dM,dpT,Mpart,pTpart found, looking for the 2D distribution");
      hdMdpT = (TH2D*)fin1->Get("fhnDeltaMass_0_proj_1_0");
      if(!hdMdpT){
      	 Printf("dMdpT not found");
      	 return;
      }
      hdMdpT->Sumw2();
      hdMdpT->Rebin2D(2, 2);
      hdMdpT->Scale(1./hdMdpT->Integral());
        
      Printf("BinWidth dMdpT = %f, %f, Integral = %.2f", hdMdpT->GetXaxis()->GetBinWidth(2), hdMdpT->GetYaxis()->GetBinWidth(2), hdMdpT->Integral());

   }
   
   
    
   TFile *fin2 = new TFile(file2);
   if(!fin2->IsOpen()){
      Printf("%s not found", file2.Data());
      return;
      
   }
   TString nameResp = "fhnMassResponse_proj_1_3";
   if(brecoLev) nameResp = "fhnMassResponse_proj_0_2";
   TH2D *hMpTpart = (TH2D*)fin2->Get(nameResp);
   if(!hMpTpart){
      Printf("particle level pT vs M not found");
      return;
   }
   hMpTpart->Sumw2();
   hMpTpart->GetXaxis()->SetTitleOffset(0.9);
   hMpTpart->GetYaxis()->SetTitleOffset(0.95);
   hMpTpart->Rebin2D(2, 2); //note that x and y are inverted wrt hdMdpT
   Printf("BinWidth X - %s = %f, Y - %s = %f, Integral = %.2f", hMpTpart->GetXaxis()->GetTitle(), hMpTpart->GetXaxis()->GetBinWidth(2), hMpTpart->GetYaxis()->GetTitle(), hMpTpart->GetYaxis()->GetBinWidth(2), hMpTpart->Integral());
   hMpTpart->GetZaxis()->SetRangeUser(0., 3);
   TCanvas *cinput = new TCanvas("cinput", "Input", 1000, 500);
   cinput->Divide(2, 1);
   cinput->cd(1);
   hMpTpart->Draw("colz");
   cinput->cd(2);
   if(hdMdpT) hdMdpT->Draw("colz");
   
   
   TCanvas *c2did = new TCanvas("c2did", "Smeared 2D histogram", 600, 600);
   
   TCanvas *cproj2did = new TCanvas("cproj2did", "Projections of the smeared 2D histograms per each dMdpT bin considered", 1000, 500);
   cproj2did->Divide(2, 1);
   
   Int_t binmass = 4;
   Double_t massbins[binmass+1] = {0., 5., 10., 15., 20.};
   TH2D *h2dSmeared = 0x0;
   if(!hdMdpT) {
      if( !hdMdpTMparpTpar ) {
      	 Printf("Fluctuation thnsparse not found");
      	 return;
      }
      Printf("Axis 2 is %s", hdMdpTMparpTpar->GetAxis(2)->GetTitle());
      Printf("axis 3 is %s;", hdMdpTMparpTpar->GetAxis(3)->GetTitle());
      Printf("X axis pTpartMpart is %s", hMpTpart->GetXaxis()->GetTitle());
      
      Int_t rebina[4] = {2,2,1,1};
      //Printf("BinWidth Before 1 - %s = %f, 2 - %s = %f, 3 - %s = %f, 4 - %s = %f,Integral = %.2f", hdMdpTMparpTpar->GetAxis(0)->GetTitle(), hdMdpTMparpTpar->GetAxis(0)->GetBinWidth(2), hdMdpTMparpTpar->GetAxis(1)->GetTitle(), hdMdpTMparpTpar->GetAxis(1)->GetBinWidth(2), hdMdpTMparpTpar->GetAxis(2)->GetTitle(), hdMdpTMparpTpar->GetAxis(2)->GetBinWidth(2), hdMdpTMparpTpar->GetAxis(3)->GetTitle(), hdMdpTMparpTpar->GetAxis(3)->GetBinWidth(2), hdMdpTMparpTpar->ComputeIntegral());
      
      hdMdpTMparpTpar = (THnSparseF*)hdMdpTMparpTpar->RebinBase(rebina);
      //hdMdpTMparpTpar->Scale(1./hdMdpTMparpTpar->ComputeIntegral());
      
      //Printf("BinWidth After 1 - %s = %f, 2 - %s = %f, 3 - %s = %f, 4 - %s = %f, Integral = %.2f", hdMdpTMparpTpar->GetAxis(0)->GetTitle(), hdMdpTMparpTpar->GetAxis(0)->GetBinWidth(2), hdMdpTMparpTpar->GetAxis(1)->GetTitle(), hdMdpTMparpTpar->GetAxis(1)->GetBinWidth(2), hdMdpTMparpTpar->GetAxis(2)->GetTitle(), hdMdpTMparpTpar->GetAxis(2)->GetBinWidth(2), hdMdpTMparpTpar->GetAxis(3)->GetTitle(), hdMdpTMparpTpar->GetAxis(3)->GetBinWidth(2), hdMdpTMparpTpar->ComputeIntegral());
      // this method still has problems with the normalization and has to be refined
      
      TH2D* hdMdpTAll = hdMdpTMparpTpar->Projection(0, 1);
      hdMdpTAll->SetName("hdMdpTAll");
      hdMdpTAll->Scale(1./hdMdpTAll->Integral());
      cinput->cd(2);
      hdMdpTAll->Draw("colz");
      Printf("Integral should sum up to %f",hdMdpTAll->Integral());
      Double_t integrals[binmass];
      TH2D* hMpTRecFluc = hdMdpTMparpTpar->Projection(2, 3); //M, pT
      hMpTRecFluc->SetName("hMpTRecFluc");
      hMpTRecFluc->Scale(1./hMpTRecFluc->Integral());
      
      Double_t totalInt = 0;
      TCanvas *cdMdpTPart = new TCanvas("cdMdpTPart", "DeltaMDeltapT in bins on M", 800, 500);
      cdMdpTPart->Divide(2,2);
      
      for(Int_t im = 0; im < binmass; im++ ){
      	 Printf("Integral RANGE %d - %d",hMpTRecFluc->GetYaxis()->FindBin(massbins[im]), hMpTRecFluc->GetYaxis()->FindBin(massbins[im+1]-1));
      	 integrals[im] = hMpTRecFluc->Integral(0, -1, hMpTRecFluc->GetYaxis()->FindBin(massbins[im]), hMpTRecFluc->GetYaxis()->FindBin(massbins[im+1]-1));
      	 Printf("Normalize to %f", integrals[im]);
      	 
      	 Int_t bin1M[2] = {hdMdpTMparpTpar->GetAxis(2)->FindBin(massbins[im]), hdMdpTMparpTpar->GetAxis(2)->FindBin(massbins[im+1]) -1 };
      	 
      	 hdMdpTMparpTpar->GetAxis(2)->SetRange(bin1M[0], bin1M[1]);
      	 hdMdpT = (TH2D*) hdMdpTMparpTpar->Projection(0, 1);
      	 hdMdpT->SetName(Form("%s_B%d", hdMdpT->GetName(), im));
      	 hdMdpT->Scale(1./hdMdpT->Integral()*integrals[im]);
      	 //hdMdpT->Rebin2D(2, 2);
      	 cdMdpTPart->cd(im+1);
      	 hdMdpT->Draw("colz");
      	 
      	 if(im == 0) Printf("BinWidth X - %s = %f, Y - %s = %f, Integral = %.2f", hdMdpT->GetXaxis()->GetTitle(), hdMdpT->GetXaxis()->GetBinWidth(2), hdMdpT->GetYaxis()->GetTitle(), hdMdpT->GetYaxis()->GetBinWidth(2), hdMdpT->Integral());
      	 Double_t partInt = hdMdpT->Integral();
      	 totalInt += partInt;
      	 Printf("@@@@@@ PartInt = %f", partInt);
      	 Int_t bin2M[2] = {hMpTpart->GetXaxis()->FindBin(massbins[im]), hMpTpart->GetXaxis()->FindBin(massbins[im+1]) -1 };
      	 hMpTpart->GetXaxis()->SetRange(bin2M[0], bin2M[1]);
      	 
      	 if(!h2dSmeared) h2dSmeared = LoopOnFluctuations(hdMdpT, hMpTpart, im); //contains the final output
      	 else h2dSmeared->Add(LoopOnFluctuations(hdMdpT, hMpTpart, im));
      }
      Printf("------- Tot Int = %f", totalInt);
      hMpTpart->GetXaxis()->SetRange(0, -1);
   } // the following method uses the integrated over pT and M dMdpT and works more or less well, but needs to be checked better because the shape of the final distribution is bumpy 
   else h2dSmeared = LoopOnFluctuations(hdMdpT, hMpTpart); //contains the final output
   
 
   c2did->cd();
   h2dSmeared->Draw("colz");
   
   cproj2did->cd(1);
   gPad->SetLogy();
   
   TH1D *hMpart = hMpTpart->ProjectionX("hMpart");
   
   TH1D *hMsmea = h2dSmeared->ProjectionX("hMsmea");
   hMsmea->SetLineColor(colors[1]);
   hMsmea->Draw();
   hMpart->Draw("sames");

   cproj2did->cd(2);
   gPad->SetLogy();
   
   TH1D *hpTpart = hMpTpart->ProjectionY("hpTpart");
   
   TH1D *hpTsmea = h2dSmeared->ProjectionY("hpTsmea");
   hpTsmea->SetLineColor(colors[1]);
   hpTsmea->Draw();
   hpTpart->Draw("sames");
   
   Printf("Integral of the output %.2f", h2dSmeared->Integral());
      	 //if(id == 0) h2id->ProjectionX()->Draw();
      	 //else if(id%10 == 0) h2id->ProjectionX()->Draw("sames");
      	 //cproj2did->cd(2);
      	 //if(id == 0) h2id->ProjectionY()->Draw();
      	 //else if(id%10 == 0) h2id->ProjectionY()->Draw("sames");
      	 
   // debugging
   
  
   TCanvas *cMinpTbins = new TCanvas("cMinpTbins", "Mass in pT bins",  1000, 500);
   cMinpTbins->Divide(2, 1);
   const Int_t ptbins = 9;
   Double_t ptlims[ptbins+1] = {-30, -20, -10, 0, 10, 20, 30, 40, 50, 100};
   TLegend *legpt = new TLegend(0.2, 0.55, 0.5, 0.9, "Scaled to 1");
   legpt->SetBorderSize(0);
   legpt->SetFillStyle(0);      
      
   for(Int_t ipt = 0; ipt < ptbins; ipt++){
      Int_t binpt[2] = {h2dSmeared->GetYaxis()->FindBin(ptlims[ipt]), h2dSmeared->GetYaxis()->FindBin(ptlims[ipt+1]) - 1};
      
      TH1D *hMbin = h2dSmeared->ProjectionX(Form("hM_pt%d", ipt), binpt[0], binpt[1]);
      hMbin->SetLineColor(colors[ipt]);
      hMbin->SetLineWidth(2);
      Double_t integral = hMbin->Integral();
      if(integral>0) hMbin->Scale(1./integral);
      legpt->AddEntry(hMbin, Form("%.0f < #it{p}_{T} < %.0f", ptlims[ipt], ptlims[ipt+1]));
      cMinpTbins->cd(1);
      if(ipt == 0) {
      	 //hMbin->GetYaxis()->SetRangeUser(0.,10);
      	 hMbin->Draw();
      	 
      }
      else hMbin->Draw("sames");
      
      binpt[0] = hMpTpart->FindBin(ptlims[ipt]);
      binpt[1] = hMpTpart->FindBin(ptlims[ipt+1])-1;
      
      TH1D *hMPartbin = hMpTpart->ProjectionX(Form("hMpart_pt%d", ipt), binpt[0], binpt[1]);
      hMPartbin->SetLineColor(colors[ipt]);
      hMPartbin->SetLineWidth(2);
      integral = hMPartbin->Integral();
      
      if(integral>0) hMPartbin->Scale(1./integral);
      else continue;
      cMinpTbins->cd(2);
      if(ipt == 0) {
      	 //hMbin->GetYaxis()->SetRangeUser(0.,10);
      	 hMPartbin->Draw();
      	 
      }
      else hMPartbin->Draw("sames");
   }
   cMinpTbins->cd(1);
   legpt->Draw();
   
   //save Cv
   SaveCv(cinput);
   SaveCv(c2did);
   SaveCv(cproj2did);
   SaveCv(cMinpTbins);
}

//________________________________________________________________________________________________

void SmearWithRandom(Int_t N = 10 /*number of random numbers drawn for the smearing*/, TString file1 = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train976/output/runlist1/AnalysisResults.root", TString strLst1 = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", TString file2 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root", TString strLst2 = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged"){
   
   TList* list1 = ReadFile(file1, strLst1);
   
   THnSparseF *hspRespFluc = (THnSparseF*) list1->FindObject("fhnDeltaMass_0");
   if(!hspRespFluc){
      Printf("THnSparse embedding not found");
      list1->ls();
      return;
   }
   // dM, dpT, Mdet, Mpart, pTdet, pTpart, pTdetAreaSub (the other is det sub)
   
   Int_t pTpart = 5, Mpart = 3, pTdet = 4, Mdet = 2, dM = 0, dpT = 1;
   //printout
   PrintTHnSparse(hspRespFluc);
   //make a projection on the (pT, M)fluc at each bin of (pT,M)rec
   
   Int_t rebinFlucSp[7] = {1, 1, 2, 2, 5, 5, 1}; 
   hspRespFluc = (THnSparseF*)hspRespFluc->RebinBase(rebinFlucSp);
   PrintTHnSparse(hspRespFluc);
   TCanvas *cSlicesFluc = new TCanvas("cSlicesFluc", "Slices of the fluctuation response matrix", 1500, 1000);
   cSlicesFluc->Divide(10, 5);
   Int_t nbinspTp = hspRespFluc->GetAxis(pTpart)->GetNbins();
   Int_t nbinsMp  = hspRespFluc->GetAxis(Mpart) ->GetNbins();
   
   TH2D *h2Fluc[nbinspTp][nbinsMp];
   
   Int_t count = 0;
   for(Int_t ipt = 0; ipt<nbinspTp; ipt++){
      for(Int_t im = 0; im<nbinsMp; im++){
      	 
      	 hspRespFluc->GetAxis(pTpart)->SetRange(ipt+1, ipt+2);
      	 hspRespFluc->GetAxis(Mpart)->SetRange(im+1, im+2);
      	 //projection on the detector level response
      	 //h2Fluc[ipt][im] = (TH2D*)hspRespFluc->Projection(pTdet,Mdet);
      	 //projection on the dM dpT
      	 h2Fluc[ipt][im] = (TH2D*)hspRespFluc->Projection(dpT,dM);
      	 h2Fluc[ipt][im]->SetName(Form("h2Fluc_%d_%d", ipt, im));
      	 
      	 if(h2Fluc[ipt][im]->GetEntries()<1e-10) {
      	    h2Fluc[ipt][im] = 0x0;
      	    continue;
      	 }
      	 count++;
      	 Double_t integralFluc = h2Fluc[ipt][im]->Integral();
      	 //if(integralFluc>0) Printf("%d, %d) Integral = %e", ipt, im, integralFluc);
      	 if(count > 50) continue;
      	 cSlicesFluc->cd(count);
      	 h2Fluc[ipt][im]->Draw("colz");
      	 
      }
      
   }
   
   //read the THnSparse with the detector response
   TList* list2 = ReadFile(file2, strLst2);
   THnSparseF *hspResp = (THnSparseF*) list2->FindObject("fhnMassResponse");
   if(!hspResp){
      Printf("THnSparse PYTHIA detector response not found");
      list2->ls();
      return;
   }
   PrintTHnSparse(hspResp);
   Int_t pTpartR = 3, MpartR = 1, pTdetR = 2, MdetR = 0;
   Int_t reduced[4] = {MdetR, MpartR, pTdetR, pTpartR}; //check thiiisss in fhnMassResponse
   hspResp = (THnSparseF*)hspResp->ProjectionND(4,reduced);
   Int_t rebinResp[4] = {2, 2, 5, 5};
   hspResp = (THnSparseF*)hspResp->RebinBase(rebinResp);
   PrintTHnSparse(hspResp, "Reduced");
   
   nbinspTp = hspResp->GetAxis(pTpartR)->GetNbins();
   nbinsMp  = hspResp->GetAxis(MpartR) ->GetNbins();
   Int_t nbinspTr = hspResp->GetAxis(pTdetR)->GetNbins();
   Int_t nbinsMr  = hspResp->GetAxis(MdetR) ->GetNbins();
   //Printf("Input DetResponse Entries = %f, Nbins pTp = %d, W = %f, Mp = %d, W = %f", hspResp->GetEntries(), nbinspTp, hspResp->GetAxis(pTpartR)->GetBinWidth(1), nbinsMp, hspResp->GetAxis(MpartR) ->GetBinWidth(1));
   
   Int_t axesFullResp[4] = {pTpart, Mpart, pTdet, Mdet};
   THnSparse *hspFullResponse = (THnSparse*) hspRespFluc->ProjectionND(4, axesFullResp);
   hspFullResponse->SetName("hspFullResponse");
   hspFullResponse->Reset();

   //loop on the response pTpart,Mpart
   //for(Int_t iptp = 0; iptp < nbinspTp; iptp++){
   for(Int_t iptp = 20; iptp < 21; iptp++){
      //for(Int_t imp = 0; imp < nbinsMp; imp++){
      for(Int_t imp = 10; imp < 11; imp++){
      	 Double_t ptp = hspResp->GetAxis(pTpartR)->GetBinCenter(iptp);
      	 Double_t mp  = hspResp->GetAxis(MpartR) ->GetBinCenter(imp);
      	 //loop on the pTrec,Mrec corresponding bins
      	 for(Int_t iptr = 0; iptr < nbinspTr; iptr++){
      	    for(Int_t imr = 0; imr < nbinsMr; imr++){
      	       //take fluctuation response of a specific (pT,M)rec bin and extract a pair of random numbers
      	       Double_t ptFl, mFl;
      	       
       	       for(Int_t i = 0; i<N; i++){
      	       	  if(!h2Fluc[iptr][imr]){
      	       	     ptFl = hspResp->GetAxis(pTdetR)->GetBinCenter(iptr);
      	       	     mFl  = hspResp->GetAxis(MdetR) ->GetBinCenter(imr);
      	       	  }
      	       	  else h2Fluc[iptr][imr]->GetRandom2(ptFl,mFl);
      	       	  Double_t arrsmear[4] = { mFl,  mp, ptFl, ptp};
      	       	  //fill a 2D histogram with the two random numbers: indeces [pTpart],[Mpart],[pTrec],[Mrec]
      	       	  hspFullResponse->Fill(arrsmear);
      	       }
      	    }
      	 }
      }
   }
   TCanvas *cFullResp = new TCanvas("cFullResp", "Full response");
   
   TFile *fout = new TFile("FullResponse.root", "recreate");
   fout->cd();
   hspFullResponse->Write();
}

//________________________________________________________________________________________________

void PrintTHnSparse(THnSparse *hsp, TString appendString){
   Printf("%s", appendString.Data());
   Printf("Nentries = %.0f", hsp->GetEntries());
   Int_t ndim = hsp->GetNdimensions();
   Printf("Ndim = %d", ndim);
   for(Int_t i = 0 ; i< ndim; i++){
      Int_t nbins = hsp->GetAxis(i)->GetNbins();
      Printf("Axis %d (%s) -> (%.3f, %.3f), N bins = %d, Bin W = %e", i, hsp->GetAxis(i)->GetTitle(), hsp->GetAxis(i)->GetBinLowEdge(1), hsp->GetAxis(i)->GetBinLowEdge(nbins), nbins, hsp->GetAxis(i)->GetBinWidth(1));
   
   }

}

//________________________________________________________________________________________________

TH2D* MakeGaussiandMdpT(Int_t mode){
   
   Int_t nbinx = 50, nbiny = 50;
   TH2D *hdMdpTGaus = new TH2D("hdMdpTGaus", "Gaus ; #delta #it{p}_{T} (GeV/c); #delta M (GeV)", nbinx, -50, 50, nbiny, -25, 25);
   Int_t fixedx = 36, fixedy = 36;
   for(Int_t ix = 0; ix< nbinx; ix++){
      for(Int_t iy = 0; iy< nbiny; iy++){
      	 Printf("Low edge %f ", hdMdpTGaus->GetBinLowEdge(iy+1));
      	 switch (mode){
      	 case 1:
      	    
      	    if(ix==iy) hdMdpTGaus->SetBinContent(ix+1, iy+1, 1);
      	    else hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    if(hdMdpTGaus->GetXaxis()->GetBinLowEdge(ix+1) < 0 || hdMdpTGaus->GetYaxis()->GetBinLowEdge(iy+1) < 0 ) hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    break;
      	    
      	 case 2:
      	    
      	    if(iy == fixedy) hdMdpTGaus->SetBinContent(ix+1, iy+1, 1);
      	    else hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    if(hdMdpTGaus->GetXaxis()->GetBinLowEdge(ix+1) < 0 || hdMdpTGaus->GetYaxis()->GetBinLowEdge(iy+1) < 0 ) hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    break;
      	    
      	 case 3:
      	    
      	    if(ix == fixedx) hdMdpTGaus->SetBinContent(ix+1, iy+1, 1);
      	    else hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    if(hdMdpTGaus->GetXaxis()->GetBinLowEdge(ix+1) < 0 || hdMdpTGaus->GetYaxis()->GetBinLowEdge(iy+1) < 0 ) hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    break;
      	    
      	 case 4:
      	    
      	    if(ix==iy) hdMdpTGaus->SetBinContent(ix+1, iy+1, 1);
      	    else hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    if(hdMdpTGaus->GetXaxis()->GetBinLowEdge(ix+1) < 0 || hdMdpTGaus->GetYaxis()->GetBinLowEdge(iy+1) < 0 ) hdMdpTGaus->SetBinContent(ix+1, iy+1, 0);
      	    if(TMath::Abs(ix-iy) < 10) {
      	       hdMdpTGaus->SetBinContent(ix+1, iy+1, 1-0.1*(TMath::Abs(ix-iy)));
      	    
      	    }
      	    break;
      	 }
      	 
      }
   }
   hdMdpTGaus->Scale(1./hdMdpTGaus->Integral());
   return hdMdpTGaus;
   
}

//________________________________________________________________________________________________
void RunSmearingDiagonale(Int_t mode, Bool_t brecoLev = kTRUE, TString file2 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/UnfoldingMatrix.root"){
   
   TH2D *hdMdpT = MakeGaussiandMdpT(mode);
   if(!hdMdpT) return;
   TFile *fin2 = new TFile(file2);
   if(!fin2->IsOpen()){
      Printf("%s not found", file2.Data());
      return;
      
   }
   TString nameResp = "fhnMassResponse_proj_1_3";
   if(brecoLev) nameResp = "fhnMassResponse_proj_0_2";
   TH2D *hMpTpart = (TH2D*)fin2->Get(nameResp);
   if(!hMpTpart){
      Printf("particle level pT vs M not found");
      return;
   }
   hMpTpart->Sumw2();
   hMpTpart->GetXaxis()->SetTitleOffset(0.9);
   hMpTpart->GetYaxis()->SetTitleOffset(0.95);
   hMpTpart->Rebin2D(2, 2); //note that x and y are inverted wrt hdMdpT
   Printf("BinWidth X - %s = %f, Y - %s = %f, Integral = %.2f", hMpTpart->GetXaxis()->GetTitle(), hMpTpart->GetXaxis()->GetBinWidth(2), hMpTpart->GetYaxis()->GetTitle(), hMpTpart->GetYaxis()->GetBinWidth(2), hMpTpart->Integral());
   hMpTpart->GetZaxis()->SetRangeUser(0., 3);
   TCanvas *cinput = new TCanvas("cinput", "Input", 1000, 500);
   cinput->Divide(2, 1);
   cinput->cd(1);
   hMpTpart->Draw("colz");
   cinput->cd(2);
   if(hdMdpT) hdMdpT->Draw("colz");
   
   TH2D *h2dSmeared = LoopOnFluctuations(hdMdpT, hMpTpart);
   
   TCanvas *c2did = new TCanvas("c2did", "Smeared 2D histogram", 600, 600);
   
   TCanvas *cproj2did = new TCanvas("cproj2did", "Projections of the smeared 2D histograms per each dMdpT bin considered", 1000, 500);
   cproj2did->Divide(2, 1);
   
   c2did->cd();
   h2dSmeared->Draw("colz");
   
   cproj2did->cd(1);
   gPad->SetLogy();
   
   TH1D *hMpart = hMpTpart->ProjectionX("hMpart");
   
   TH1D *hMsmea = h2dSmeared->ProjectionX("hMsmea");
   hMsmea->SetLineColor(colors[1]);
   hMsmea->Draw();
   hMpart->Draw("sames");

   cproj2did->cd(2);
   gPad->SetLogy();
   
   TH1D *hpTpart = hMpTpart->ProjectionY("hpTpart");
   
   TH1D *hpTsmea = h2dSmeared->ProjectionY("hpTsmea");
   hpTsmea->SetLineColor(colors[1]);
   hpTsmea->Draw();
   hpTpart->Draw("sames");
   
}

void PrintCompareTH2(TH2* h){
   Printf("%s \n X axis %s Nbins %d, W = %e, (%.3f, %.3f)",h->GetName(), h->GetXaxis()->GetTitle(), h->GetXaxis()->GetNbins(), h->GetXaxis()->GetBinWidth(1), h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetNbins()+1));
   Printf("Y axis %s Nbins %d, W = %e, (%.3f, %.3f)", h->GetYaxis()->GetTitle(), h->GetYaxis()->GetNbins(), h->GetYaxis()->GetBinWidth(1), h->GetYaxis()->GetBinLowEdge(1), h->GetYaxis()->GetBinLowEdge(h->GetYaxis()->GetNbins()+1));
   
}
