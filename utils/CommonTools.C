#ifndef CommonTools_C
#define CommonTools_C
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH2D.h>
#include <TClass.h>
#include <TGraphAsymmErrors.h>
#include <TPaveText.h>
#include <THnSparse.h>

Int_t colorSeriesBlue[] = {kCyan-3, kCyan-7, kCyan-6, kCyan-2, kCyan+4, kBlue, kBlue-7, kBlue-1, kBlue-6, kBlue+3};
Int_t colorSeriesGreen[] = {kGreen-9, kGreen+4, kGreen+1, kGreen-8, kGreen-6, kGreen+2, kGreen+3, kGreen+4 };
void SaveCv(TCanvas* c, TString suffix = "", Int_t format = 3);   
TList* ReadFile(TString strIn, TString strLst);
TList* ReadFile(TString strIn, TString strDir, TString strLst);
TObject* ReadObjInFile(TString strIn, TString strLst, TString objname);
TH1* CompareDistributions(TVirtualPad* pad, TH1* h1, TH1* h2);
void CalculatePads(Int_t n, Int_t&nx, Int_t&ny, Int_t&dx, Int_t&dy, Int_t perrow = 1, Int_t stdd = 400);
Int_t UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type = "TH1F", Bool_t onlycheck = kFALSE);
Int_t FindBinInArray(Double_t q, Double_t array[], Int_t dim, Bool_t verbose = kFALSE);
TH1F* ConvertTH1DinF(TH1D *h1d);
TH2F* ConvertTH2DinF(TH2D *h2d);
TH2F* TransformAxisRanges(TH2F *h2orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew, Int_t ynbinsnew, Float_t yminnew, Float_t ymaxnew);
TH1F* TransformAxisRanges(TH1F *h1orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew);
void DefineNewAxes(TH2 *hOrig, TH2 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Float_t deltaDownY, Float_t deltaUpY, Int_t& newNbinsX, Int_t& newNbinsY, Float_t& minX, Float_t& maxX, Float_t& minY, Float_t& maxY );
void DefineNewAxes(TH1 *hOrig, TH1 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Int_t& newNbinsX, Float_t& minX, Float_t& maxX );
TH1D*** DrawProjFrom2D(const Int_t n2dh, TH2D** hist, Int_t axpj, const Int_t axrange, Int_t nbinsrange, Double_t rangelims[], const char* namenew, const Int_t cls[], const Int_t mks[], Int_t nametype = 0, Double_t eps = 0.001);

void PrintHisto2(TH2* h1);

TPaveText** GetPavePtBins(Int_t nptbins, Double_t *lims, Int_t x1 = 0.3, Int_t y1 = 0.8, Int_t x2 = 0.8, Int_t y2 = 0.9);

TH1D** GetPythiaOrThnSpaseProjections(Int_t axrange, Int_t axproj, Double_t binWaxproj, const Int_t nbins, Double_t binlims[], TString finame, TString listname, TString hname, TString pjbasename);

TH2D* RebinVariableWidth2D(const Int_t nbinsvarwX, Double_t varwlimsX[], const Int_t nbinsvarwY, Double_t varwlimsY[], TH2D* hinput, Double_t ex = 0.001, Double_t ey = 0.001);

//method implementation

void SaveCv(TCanvas* c,TString suffix, Int_t format){
   if(format > 0) c->SaveAs(Form("%s%s.png",c->GetName(),suffix.Data()));
   if(format > 1) c->SaveAs(Form("%s%s.pdf",c->GetName(),suffix.Data()));
   if(format > 2) c->SaveAs(Form("%s%s.eps",c->GetName(),suffix.Data()));
   if(format > 3 || format==0) c->SaveAs(Form("%s%s.root",c->GetName(),suffix.Data()));
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TList* ReadFile(TString strIn, TString strLst){
   
   TFile *f = new TFile(strIn.Data());
   if(!f->IsOpen()){
      Printf("File %s not found", strIn.Data());
      return 0x0;
   }
   TList *lst = static_cast<TList*>(f->Get(strLst.Data()));
   if(!lst){
      Printf("Error, list %s not found", strLst.Data());
      f->ls();
      return 0x0;
   }
   //Printf("Returning list %s, %p", strLst.Data(), lst);
   f->Close();
   return lst;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TList* ReadFile(TString strIn, TString strDir, TString strLst){
   
   TFile *f = new TFile(strIn.Data());
   if(!f->IsOpen()){
      Printf("File %s not found", strIn.Data());
      return 0x0;
   }
   TDirectoryFile *df = (TDirectoryFile*) f->Get(strDir);
   if(!df){
      Printf("Directory %s not found", strDir.Data());
      return 0x0;
   }
   
   TList *lst = static_cast<TList*>(df->Get(strLst.Data()));
   if(!lst){
      Printf("Error, list %s not found", strLst.Data());
      return 0x0;
   }
   
   f->Close();
   
   return lst;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TObject* ReadObjInFile(TString strIn, TString strLst, TString objname){
   TList* l = ReadFile(strIn, strLst);
   if(!l) return 0x0;
   TObject *obj = l->FindObject(objname);
   if(!obj){
      Printf("%s not found", objname.Data());
   } else {
      Printf("%s is a %s", objname.Data(), obj->IsA()->GetName());
   }
   return obj;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------

TH1* CompareDistributions(TVirtualPad* pad, TH1* h1, TH1* h2){
   h1->SetName(Form("%s1", h1->GetName()));
   h2->SetName(Form("%s2", h2->GetName()));
   Printf("Compare %s and %s", h1->GetName(), h2->GetName());
   h1->Scale(1./h1->Integral("width"));
   h2->Scale(1./h2->Integral("width"));
   
   Int_t nbins1 = h1->GetNbinsX(), nbins2 = h1->GetNbinsX();
   Double_t edges1[2] = {h1->GetBinLowEdge(1), h1->GetBinLowEdge(nbins1+1)}, edges2[2] = {h2->GetBinLowEdge(1), h2->GetBinLowEdge(nbins2+1)};
   
   Int_t diffbins = nbins1 - nbins2;
   if(diffbins!=0){
      
      
   }
   TH1* hr=(TH1*)h1->Clone("hr");
   hr->Divide(h2);
   
   pad->cd();
   h1->Draw("sames");
   h2->Draw("sames");
   
   return hr;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
void CalculatePads(Int_t n, Int_t&nx, Int_t&ny, Int_t&dx, Int_t&dy, Int_t perrow, Int_t stdd){
   
   Int_t percolumn = n/perrow;
   if(n%perrow > 0) percolumn++;
   nx = percolumn;
   ny = perrow;
   dx = stdd*nx;
   dy= stdd*ny;
   
   return;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
void SaveGraphInCanvasToFile(TCanvas *corigin, TString graphname, TString newname, TString fileoutname){
   
   if(!corigin){
      Printf("Canvas of origin is null");
      return;
   }
   corigin->ls();
   
   if(graphname.IsNull() && fileoutname.IsNull()) return;
   
   TFile *fout=new TFile(fileoutname, "recreate");
   TClass* cl = corigin->FindObject(graphname)->IsA();
   TString type = cl ->GetName();
   
   Printf("Got a %s", type.Data());
   Printf("Types implemented: TGraphAsymmErrors"); //update
   //TCanvas *myc = new TCanvas("myc", "Draw what you have");
   if(type = "TGraphAsymmErrors"){
      TGraphAsymmErrors* g = (TGraphAsymmErrors*)(corigin->FindObject(graphname));
      if(!g) {
      	 Printf(" %s, type %s not found", graphname.Data(), type.Data());
      	 return;
      	 
      }
      //myc->cd();
      g->SetName(newname);
      
      fout->cd();
      g->Write();
      
      TH1F* h = g->GetHistogram();
      h->SetName(Form("h%s", newname.Data()));
      fout->cd();
      h->Write();
   }
   /*
   TH1F *htest = new TH1F("fanculo","Fanculo", 100, 1, 5);
   htest->FillRandom("pol0");
   fout->cd();
   htest->Write();
   */
   Printf("Content");
   fout->ls();
   
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type, Bool_t onlycheck){
   
   /// returns -1 if input is not valid, 0 if no change in the histograms is needed, 1 if histograms are redefined
   
   if(!h1 || !h2){
      Printf("Need valid input!");
      return -1;
   }
   Printf("UniformTH1FForDivide:: THINK: if the histograms are ratios this method does not do what you want");
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h2->GetNbinsX();
   Double_t binW1 = h1->GetBinWidth(3), binW2 = h2->GetBinWidth(3);
   Double_t low1 = h1->GetBinLowEdge(1), low2 = h2->GetBinLowEdge(1);
   Double_t up1 = h1->GetBinLowEdge(nBins1+1), up2 = h2->GetBinLowEdge(nBins2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBins2, binW2, low2, up2);
   if(((nBins1 == nBins2) && ((binW1-binW2) < 1e-6) &&  ((low1-low2) < 1e-6) && ((up1-up2) < 1e-6)) || onlycheck) {
      if(onlycheck) Printf("Check performed");
      else Printf("Histograms can be divided");
      h1bis = (TH1*)h1->Clone(Form("%sb", h1->GetName()));
      h2bis = (TH1*)h2->Clone(Form("%sb", h2->GetName()));
      return 0;
   }
   
   Double_t commonlow = TMath::Min(low1, low2), commonup = TMath::Max(up1, up2);
   Double_t commonW = TMath::Max(binW1, binW2); // the wider bin width it the one that will be used for both
   Int_t commonNBins = (commonup - commonlow)/commonW;
   Printf("Common) N bins = %d, Width = %f, Range = %f - %f", commonNBins, commonW, commonlow, commonup);
   Double_t minWidth = TMath::Min(binW1, binW2); // pick the smaller among the two bin widths...
      
   Int_t rebin = (Int_t)(commonW/minWidth); //... and calculate the number of bins to be combined
   Printf("%f/%f = %f or %d", commonW, minWidth, commonW/minWidth, rebin);
   
   if(minWidth == binW1) h1->Rebin(rebin); //if the smaller width was for h1, rebin h1
   else h2->Rebin(rebin); //otherwise rebin h2
   
   
   //TCanvas *ctest = new TCanvas(Form("c%s%s", h1->GetName(), h2->GetName() ), Form("canvas %s%s", h1->GetName(), h2->GetName() ));
   //ctest->cd();
   //h1->Draw();
   //h2->Draw("sames");
   
   //Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsX(), h2->GetName(), h2->GetNbinsX());
   //create the new histograms with correct range
   if(type == "TH1F"){
      h1bis = new TH1F(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      h2bis = new TH1F(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
   } else {
      if(type == "TH1D"){
      	 h1bis = new TH1D(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      	 h2bis = new TH1D(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle()), commonNBins, commonlow, commonup);
      } else {
      	 Printf("%s Not defined!", type.Data());
      	 return -1;
      }
   }
   h1bis->Sumw2();
   h1bis->SetLineColor(h1->GetLineColor());
   h1bis->SetLineWidth(h1->GetLineWidth());
   h1bis->SetMarkerColor(h1->GetMarkerColor());
   h1bis->SetMarkerStyle(h1->GetMarkerStyle());
   
   h2bis->Sumw2();
   h2bis->SetLineColor(h2->GetLineColor());
   h2bis->SetLineWidth(h2->GetLineWidth());
   h2bis->SetMarkerColor(h2->GetMarkerColor());
   h2bis->SetMarkerStyle(h2->GetMarkerStyle());
   
   //refill the original histograms using the same binning
   for(Int_t i = 0; i < commonNBins; i++){
      
      //First histogram
      Float_t bincentre = h1bis ->GetBinCenter(i+1);
      Int_t fillingbin = h1->FindBin(bincentre);
      Double_t fillWith = 0, error = 0;
      if(fillingbin>0 && fillingbin <= nBins1) { // exclude over- / under-flow bins
      	 //Printf("Here, debug %d", nBins1);
      	 fillWith = (Double_t)h1->GetBinContent(fillingbin);
      	 error = (Double_t)h1->GetBinError(fillingbin);
      }
      h1bis ->SetBinContent(i+1, fillWith);
      h1bis ->SetBinError(i+1, error);
      if(i%10 == 0) Printf("Filling %s x = %.3f new bin %d with old bin %d content %e [...]", h1->GetName(), bincentre, i+1, fillingbin, fillWith);
      
      //second histogram
      bincentre = h2bis ->GetBinCenter(i+1);
      fillingbin = h2->FindBin(bincentre);
      fillWith = 0; error = 0;
      if(fillingbin>0 && fillingbin <= nBins2) {
      	 fillWith = h2->GetBinContent(fillingbin);
      	 error = h2->GetBinError(fillingbin);
      }
      h2bis ->SetBinContent(i+1, fillWith);
      h2bis ->SetBinError(i+1, error);
       if(i%10 == 0) Printf("Filling %s x = %.3f new bin %d with old bin %d content %e [...]", h2->GetName(), bincentre, i+1, fillingbin, fillWith);
   }
   return 1;
   
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t FindBinInArray(Double_t q, Double_t array[], Int_t dim, Bool_t verbose){
   
   if((q > array[dim-1]) || (q<array[0])) {
      if(verbose) Printf("Out of bound");
      return -1;
   }
   for(Int_t i = 0; i<dim-1; i++ ){
   	   //Printf("i = %d", i);
      if(q >= array[i] && q < array[i+1]) return i;
      
   }
   return -2;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TH1F* ConvertTH1DinF(TH1D *h1d){
   //deletes the original histogram
   Int_t nbins = h1d->GetNbinsX();
   TH1F *h1f = new TH1F(Form("%sF", h1d->GetName()), Form("%s;%s;%s", h1d->GetTitle(), h1d->GetXaxis()->GetTitle(), h1d->GetYaxis()->GetTitle()), nbins, h1d->GetBinLowEdge(1), h1d->GetBinLowEdge(nbins+1));
   
   h1f->Sumw2();
   
   
   for(Int_t i = 0; i<nbins; i++){
      h1f->SetBinContent(i+1, h1d->GetBinContent(i+1));
      h1f->SetBinError(i+1, h1d->GetBinError(i+1));
   }
   
   //delete h1d;
   //h1d = 0;
   return h1f;
   
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TH2F* ConvertTH2DinF(TH2D *h2d){
   //deletes the original histogram
   Int_t nbinsX = h2d->GetNbinsX();
   Int_t nbinsY = h2d->GetNbinsY();
   TH2F *h2f = new TH2F(Form("%sF", h2d->GetName()), Form("%s;%s;%s", h2d->GetTitle(), h2d->GetXaxis()->GetTitle(), h2d->GetYaxis()->GetTitle()), nbinsX, h2d->GetXaxis()->GetBinLowEdge(1), h2d->GetXaxis()->GetBinLowEdge(nbinsX+1), nbinsY, h2d->GetYaxis()->GetBinLowEdge(1), h2d->GetYaxis()->GetBinLowEdge(nbinsY+1));
   
   h2f->Sumw2();
   
   
   for(Int_t i = 0; i<nbinsX; i++){
      for(Int_t j = 0; j<nbinsY; j++){
      	 Double_t content = h2d->GetBinContent(i+1, j+1);
      	 h2f->SetBinContent(i+1, j+1, content);
      	 //Printf("Bin X %d, Y %d Content = %f ---> Bin X %d, Y %d Content = %f", i+1, j+1, content, i+1, j+1, h2f->GetBinContent(i+1, j+1));
      	 h2f->SetBinError(i+1, j+1, h2d->GetBinError(i+1, j+1));
      }
   }
   //delete h2d;
   //h2d = 0;
   return h2f;
   
   
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TH2F* TransformAxisRanges(TH2F *h2orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew, Int_t ynbinsnew, Float_t yminnew, Float_t ymaxnew){
   
   TH2F* h2new = new TH2F(Form("%sNew", h2orig->GetName()), Form("%sNew;%s;%s", h2orig->GetTitle(),  h2orig->GetXaxis()->GetTitle(),  h2orig->GetYaxis()->GetTitle()), xnbinsnew, xminnew, xmaxnew, ynbinsnew, yminnew, ymaxnew);
   h2new->Sumw2();
   
   //Define nbinx nbiny number of bins original histos
   Int_t nbinx = h2orig->GetNbinsX(), nbiny = h2orig->GetNbinsY();
   
   //fill the new histo with the original histogram using the same binning
   for(Int_t i = 0; i < xnbinsnew; i++){
      for(Int_t j = 0; j < ynbinsnew; j++){
      	 
      	 Float_t bincentrex = h2new ->GetXaxis()->GetBinCenter(i+1);
      	 Float_t bincentrey = h2new ->GetYaxis()->GetBinCenter(j+1);
      	 Int_t fillingbin = h2orig->FindBin(bincentrex, bincentrey);
      	 Int_t binx = h2orig->GetXaxis()->FindBin(bincentrex);
      	 Int_t biny = h2orig->GetYaxis()->FindBin(bincentrey);
      	 Double_t fillWith = 0, error = 0;
      	 //check that we're within range of the original histogram, otherwise it will be filled with zero
      	 if(fillingbin>0 && binx <= nbinx && biny <= nbiny) {
      	    //Printf("Here, debug X bins = %d, Y bins = %d", nbinx, nbiny);
      	    fillWith = (Double_t)h2orig->GetBinContent(fillingbin);
      	    error = (Double_t)h2orig->GetBinError(fillingbin);
      	 }
      	 h2new ->SetBinContent(i+1,j+1, fillWith);
      	 h2new ->SetBinError(i+1, j+1, error);
      	 //Printf("Filling %s x = %.3f y = %.3f new bin %d, %d with old bin glob %d, %d, %d content %e", h2orig->GetName(), bincentrex, bincentrey, i+1, j+1, fillingbin, binx, biny, fillWith);
      	 
      }
   }
   return h2new;

}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TH1F* TransformAxisRanges(TH1F *h1orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew){
   
   TH1F* h1new = new TH1F(Form("%sNew", h1orig->GetName()), Form("%sNew;%s;%s", h1orig->GetTitle(),  h1orig->GetXaxis()->GetTitle(),  h1orig->GetYaxis()->GetTitle()), xnbinsnew, xminnew, xmaxnew);
   h1new->Sumw2();
   
   //Define nbinx nbiny number of bins original histos
   Int_t nbinx = h1orig->GetNbinsX();
   
   //fill the new histo with the original histogram using the same binning
   for(Int_t i = 0; i < xnbinsnew; i++){
      	 
      	 Float_t bincentrex = h1new->GetBinCenter(i+1);
      	 Int_t fillingbin = h1orig->FindBin(bincentrex);
      	 Double_t fillWith = 0, error = 0;
      	 //check that we're within range of the original histogram, otherwise it will be filled with zero
      	 if(fillingbin>0 && fillingbin <= nbinx) {
      	    //Printf("Here, debug X bins = %d, Y bins = %d", nbinx, nbiny);
      	    fillWith = (Double_t)h1orig->GetBinContent(fillingbin);
      	    error = (Double_t)h1orig->GetBinError(fillingbin);
      	 }
      	 h1new ->SetBinContent(i+1, fillWith);
      	 h1new ->SetBinError(i+1,  error);
      	 //Printf("Filling %s x = %.3f y = %.3f new bin %d, %d with old bin glob %d, %d, %d content %e", h2orig->GetName(), bincentrex, bincentrey, i+1, j+1, fillingbin, binx, biny, fillWith);
      	 
   }
   return h1new;

}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
void DefineNewAxes(TH2 *hOrig, TH2 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Float_t deltaDownY, Float_t deltaUpY, Int_t& newNbinsX, Int_t& newNbinsY, Float_t& minX, Float_t& maxX, Float_t& minY, Float_t& maxY ){
         //lower edge
      minX = deltaDownX;
      minY = deltaDownY;
      Int_t minimumBx = hdelta->GetXaxis()->FindBin(deltaDownX); 
      Int_t zeroBx = hdelta->GetXaxis()->FindBin(0.);
      Int_t minimumBy = hdelta->GetYaxis()->FindBin(deltaDownY);
      Int_t zeroBy = hdelta->GetYaxis()->FindBin(0.);
      Int_t dBinXD = zeroBx - minimumBx;
      Int_t dBinYD = zeroBy - minimumBy;
      Int_t maximumBx = hdelta->GetXaxis()->FindBin(deltaUpX); 
      Int_t maximumBy = hdelta->GetYaxis()->FindBin(deltaUpY);
      Int_t dBinXU = maximumBx - zeroBx;
      Int_t dBinYU = maximumBy - zeroBy; //x = Delta pT, y = delta M
      Printf("Minimum of delta M delta %d-%d = %d, pT  %d-%d = %d", minimumBy, zeroBy, dBinYD, minimumBx, zeroBx, dBinXD);
   
      //original nbins
      Int_t binmaxX = hOrig->GetNbinsX();
      Int_t binmaxY = hOrig->GetNbinsY();
      // upper edge
      maxX = hOrig->GetXaxis()->GetBinLowEdge(binmaxX+1) + deltaUpX;
      maxY = hOrig->GetYaxis()->GetBinLowEdge(binmaxY+1) + deltaUpY;

      newNbinsX = binmaxX + dBinXD + dBinXU;
      newNbinsY = binmaxY + dBinYD + dBinYU;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
void DefineNewAxes(TH1 *hOrig, TH1 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Int_t& newNbinsX, Float_t& minX, Float_t& maxX ){
         //lower edge
      minX = deltaDownX;
      Int_t minimumBx = hdelta->GetXaxis()->FindBin(deltaDownX); 
      Int_t zeroBx = hdelta->GetXaxis()->FindBin(0.);
      Int_t dBinXD = zeroBx - minimumBx;
      Int_t maximumBx = hdelta->GetXaxis()->FindBin(deltaUpX); 
      Int_t dBinXU = maximumBx - zeroBx;
      Printf("Minimum of delta M delta %d-%d = %d",  minimumBx, zeroBx, dBinXD);
   
      //original nbins
      Int_t binmaxX = hOrig->GetNbinsX();
      // upper edge
      maxX = hOrig->GetXaxis()->GetBinLowEdge(binmaxX+1) + deltaUpX;

      newNbinsX = binmaxX + dBinXD + dBinXU;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
void PrintHisto2(TH2* h1){
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h1->GetNbinsY();
   Double_t binW1 = h1->GetXaxis()->GetBinWidth(3), binW2 = h1->GetYaxis()->GetBinWidth(3);
   Double_t low1 = h1->GetXaxis()->GetBinLowEdge(1), low2 = h1->GetYaxis()->GetBinLowEdge(1);
   Double_t up1 = h1->GetXaxis()->GetBinLowEdge(nBins1+1), up2 = h1->GetYaxis()->GetBinLowEdge(nBins2+1);
   Printf("X) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("Y) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins2, binW2, low2, up2);
   
   return;

}

//---------------------------------------------------------------------------------------------------------------------------------------------------------
TPaveText** GetPavePtBins(const Int_t nptbins, Double_t *lims, Int_t x1, Int_t y1, Int_t x2, Int_t y2){
   
   TPaveText** pvtxt = new TPaveText*[nptbins];
   
   for(Int_t ip = 0; ip<nptbins; ip++){
      pvtxt[ip] = new TPaveText(x1, y1, x2, y2, "NDC");
      pvtxt[ip]->SetBorderSize(0);
      pvtxt[ip]->SetFillStyle(0);
      pvtxt[ip]->AddText(Form("%.0f < #it{p}_{T} < %.0f GeV/#it{c}", lims[ip], lims[ip+1]));
      Printf("pave %p", pvtxt[ip]);
   }
   
   return pvtxt;
}
//--------------------------------------------------------------------

TH1D** GetPythiaOrThnSpaseProjections(Int_t axrange, Int_t axproj, Double_t binWaxproj, const Int_t nbins, Double_t binlims[], TString finame, TString listname, TString hname, TString pjbasename){
	
	TList *lin = ReadFile(finame, listname);
	if(!lin){
		Printf("List %s not found", listname.Data());
		return 0x0;
	}
	
	THnSparseF *hsph = (THnSparseF*)lin->FindObject(hname.Data());
	if(!hsph) {
		Printf("%s not found", hname.Data());
		lin->ls();
		return 0x0;
	}
	Int_t ndim = hsph->GetNdimensions();
	
	if(axrange < 0 || axrange > ndim-1){
		Printf("Check axis for binning %d (ndims %d) failed", axrange, ndim);
		return 0x0;
	}
	
	if(axproj < 0 || axproj > ndim-1){
		Printf("Check axis for projection %d (ndims %d) failed", axproj, ndim);
		return 0x0;
	}
	Printf("Binning on %s, projecting on %s", hsph->GetAxis(axrange)->GetTitle(), hsph->GetAxis(axproj)->GetTitle());
	
	Int_t rebinar[ndim];
	for(Int_t id = 0; id<ndim; id++){
		if(id != axproj) rebinar[id] = 1;
		else{
			Double_t currentW = hsph->GetAxis(axproj)->GetBinWidth(1);
			rebinar[id] = (Int_t)(binWaxproj/currentW);
			Printf("Rebinning from %.3f to %.3f -> reb = %d", currentW, binWaxproj, rebinar[id]);
		}
	}
	hsph = (THnSparseF*)hsph->Rebin(rebinar);
	
	TH1D** hproj = new TH1D*[nbins];
	for(Int_t ipt = 0; ipt < nbins; ipt++){
		hproj[ipt] = 0x0;
		Int_t binrange[2] = {hsph->GetAxis(axrange)->FindBin(binlims[ipt]), hsph->GetAxis(axrange)->FindBin(binlims[ipt+1] - 0.1)};
		
		Printf("PtBin = %d, Bin range %f-%f -> %d - %d", ipt, binlims[ipt], binlims[ipt+1] - 0.1, binrange[0], binrange[1]);
		
		hsph->GetAxis(axrange)->SetRange(binrange[0], binrange[1]);
		
		hproj[ipt] = hsph->Projection(axproj);
		hproj[ipt]->SetName(Form("%s_%d", pjbasename.Data(), ipt));
		Printf("Pointer %s %p", hproj[ipt]->GetName(), hproj[ipt]);
	}
	
	return hproj;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------

TH1D** GetProjectionsTH3(Int_t axrange, Int_t axproj, Double_t binWaxproj, const Int_t nbins, Double_t binlims[], TString filename, TString listname, TString hname, TString pjbasename){
	
	Printf("Reading %s, List %s", filename.Data(), listname.Data());
	TList *list = ReadFile(filename, listname);
	
	
	if(!list){
		return 0;
	}
	
	TObject* obj = list->FindObject(hname.Data());
	if(!obj) {
		Printf("No obj %s found", hname.Data());
		
	}
	TH3F *h3 = 0x0;
	TH1D **hM = new TH1D*[nbins];
	
	TClass *objclas = obj->IsA();
	TString clname = objclas->GetName();
	Printf("TClass is %s", clname.Data());
	if(clname == "TH3D"){
		Printf("Reading TH3D");
		h3 = (TH3F*)list->FindObject(hname.Data());
		if(!h3){
			list->ls();
			return 0;
		} else {
			Printf("Info, TH3F %s found", hname.Data());
		}
	} else return 0;
	
	Int_t rebinar[3];
	for(Int_t id = 0; id<3; id++){
		if(id != axproj) rebinar[id] = 1;
		else{
			Double_t currentW = h3->GetXaxis()->GetBinWidth(1);
			if(axproj == 1) currentW = h3->GetYaxis()->GetBinWidth(1);
			if(axproj == 2) currentW = h3->GetZaxis()->GetBinWidth(1);
			rebinar[id] = (Int_t)(binWaxproj/currentW);
			Printf("Rebinning from %.3f to %.3f -> reb = %d", currentW, binWaxproj, rebinar[id]);
		}
	}
	h3 = (TH3F*)h3->Rebin3D(rebinar[0], rebinar[1], rebinar[2]);
	
	for(Int_t ipt = 0; ipt < nbins; ipt++){
		if(axproj == 0) {
			if(axrange == 1){
				Int_t binrange[2] = {h3->GetYaxis()->FindBin(binlims[ipt]), h3->GetYaxis()->FindBin(binlims[ipt+1] - 0.001)};
				
				hM[ipt] = h3->ProjectionX(Form("%sPj%d_%.0f_%.0f", pjbasename.Data(), axproj, binlims[ipt], binlims[ipt+1]), binrange[0], binrange[1], 0,  -1);
				
			} else 
			if(axrange == 2){
				
				Int_t binrange[2] = {h3->GetZaxis()->FindBin(binlims[ipt]), h3->GetZaxis()->FindBin(binlims[ipt+1] - 0.001)};
				
				hM[ipt] = h3->ProjectionX(Form("%sPj%d_%.0f_%.0f", pjbasename.Data(), axproj, binlims[ipt], binlims[ipt+1]), 0, -1, binrange[0], binrange[1]);
				
			} else {
				Printf("Error range and projection on the same axis!");
				return 0;	
			}
			
			
		}
		if(axproj == 1)  {
			if(axrange == 0){
				Int_t binrange[2] = {h3->GetXaxis()->FindBin(binlims[ipt]), h3->GetXaxis()->FindBin(binlims[ipt+1] - 0.001)};
				
				hM[ipt] = h3->ProjectionY(Form("%sPj%d_%.0f_%.0f", pjbasename.Data(), axproj, binlims[ipt], binlims[ipt+1]), binrange[0], binrange[1], 0,  -1);
				
			} else 
			if(axrange == 2){
				
				Int_t binrange[2] = {h3->GetZaxis()->FindBin(binlims[ipt]), h3->GetZaxis()->FindBin(binlims[ipt+1] - 0.001)};
				
				hM[ipt] = h3->ProjectionY(Form("%sPj%d_%.0f_%.0f", pjbasename.Data(), axproj, binlims[ipt], binlims[ipt+1]), 0, -1, binrange[0], binrange[1]);
				
			} else {
				Printf("Error range and projection on the same axis!");
				return 0;	
			}
			
		}
		if(axproj == 2) {
			if(axrange == 0){
				Int_t binrange[2] = {h3->GetXaxis()->FindBin(binlims[ipt]), h3->GetXaxis()->FindBin(binlims[ipt+1] - 0.001)};
				
				hM[ipt] = h3->ProjectionZ(Form("%sPj%d_%.0f_%.0f", pjbasename.Data(), axproj, binlims[ipt], binlims[ipt+1]), binrange[0], binrange[1], 0,  -1);
				
			} else 
			if(axrange == 1){
				
				Int_t binrange[2] = {h3->GetYaxis()->FindBin(binlims[ipt]), h3->GetYaxis()->FindBin(binlims[ipt+1] - 0.001)};
				
				hM[ipt] = h3->ProjectionZ(Form("%sPj%d_%.0f_%.0f", pjbasename.Data(), axproj, binlims[ipt], binlims[ipt+1]), 0, -1, binrange[0], binrange[1]);
				
			} else {
				Printf("Error range and projection on the same axis!");
				return 0;	
			}
			
		}
	}
	return hM;
	
}

TH2D* RebinVariableWidth2D(const Int_t nbinsvarwX, Double_t varwlimsX[], const Int_t nbinsvarwY, Double_t varwlimsY[], TH2D* hinput, Double_t ex, Double_t ey){
	
	//checks
	Int_t nbinsXorig = hinput->GetNbinsX();
	if(nbinsvarwX > nbinsXorig){
		Printf("New number of bins requested %d > than the original number of bins %d X", nbinsvarwX, nbinsXorig);
		return 0x0;
	}
	Double_t minXorig = hinput->GetXaxis()->GetBinLowEdge(1);
	if(minXorig > varwlimsX[nbinsvarwX-1]) {
		Printf("Minimum of the original histogram %f > than the maximum of the new X axis %f", minXorig, varwlimsX[nbinsvarwX-1]);
		return 0x0;
	}
	Double_t maxXorig = hinput->GetXaxis()->GetBinLowEdge(nbinsXorig+1);
	if(maxXorig < varwlimsX[0]) {
		Printf("Maximum of the original histogram %f < than the minimum of the new X axis %f", maxXorig, varwlimsX[0]);
		return 0x0;
	}
	
	Int_t nbinsYorig = hinput->GetNbinsY();
	if(nbinsvarwY > nbinsYorig){
		Printf("New number of bins requested %d > than the original number of bins %d Y", nbinsvarwY, nbinsYorig);
		return 0x0;
	}
	Double_t minYorig = hinput->GetYaxis()->GetBinLowEdge(1);
	if(minYorig > varwlimsY[nbinsvarwY-1]) {
		Printf("Minimum of the original histogram %f > than the maximum of the new Y axis %f", minYorig, varwlimsY[nbinsvarwY-1]);
		return 0x0;
	}
	Double_t maxYorig = hinput->GetYaxis()->GetBinLowEdge(nbinsYorig+1);
	if(maxYorig < varwlimsY[0]) {
		Printf("Maximum of the original histogram %f < than the minimum of the new Y axis %f", maxYorig, varwlimsY[0]);
		return 0x0;
	}
	
	TH2D* houtput = new TH2D(Form("%sRb", hinput->GetName()), Form("%sRb; %s; %s",  hinput->GetTitle(),  hinput->GetXaxis()->GetTitle(), hinput->GetYaxis()->GetTitle()), nbinsvarwX, varwlimsX, nbinsvarwY, varwlimsY);
	
	houtput->Sumw2();
	
	for(Int_t i = 0; i< nbinsvarwX; i++){
		for(Int_t j = 0; j< nbinsvarwY; j++){
			
			//determine how many bins of the original histograms have to be integrated
			Int_t binrangeOrigX[2] = {hinput->GetXaxis()->FindBin(varwlimsX[i]), hinput->GetXaxis()->FindBin(varwlimsX[i+1] - ex)};
			Int_t binrangeOrigY[2] = {hinput->GetYaxis()->FindBin(varwlimsY[j]), hinput->GetYaxis()->FindBin(varwlimsY[j+1] - ey)};

			Printf("Range orig X = %d-%d, Y = %d-%d", binrangeOrigX[0], binrangeOrigX[1], binrangeOrigY[0], binrangeOrigY[1]);
			Double_t fillwith = 0, fillerr = 0;
			fillwith = hinput->IntegralAndError(binrangeOrigX[0], binrangeOrigX[1], binrangeOrigY[0], binrangeOrigY[1], fillerr);
			
			houtput->SetBinContent(i+1, j+1, fillwith);
			houtput->SetBinError(i+1, j+1, fillerr);
			
			Printf("Filled bin %d - %d, range %.0f-%.0f, %.0f-%.0f", i+1, j+1, houtput->GetXaxis()->GetBinLowEdge(i+1), houtput->GetXaxis()->GetBinLowEdge(i+2), houtput->GetYaxis()->GetBinLowEdge(j+1), houtput->GetYaxis()->GetBinLowEdge(j+2));
		}
	
	}
	
	return houtput;

}


TH1D*** DrawProjFrom2D(const Int_t n2dh, TH2D** hist, Int_t axpj, const Int_t axrange, Int_t nbinsrange, Double_t rangelims[], const char* namenew, const Int_t cls[], const Int_t mks[], Int_t nametype, Double_t eps){
	Bool_t doratio = kTRUE;
	Int_t nx, ny, dx, dy;
	
	if(nbinsrange>3) CalculatePads(nbinsrange, nx, ny, dx, dy, 2);
	else CalculatePads(nbinsrange, nx, ny, dx, dy);
	
	
	TCanvas *cProj = new TCanvas(Form("cProj%s", namenew), "Projections", dx, dy);
	cProj->Divide(nx, ny);
	TCanvas *cRati = new TCanvas(Form("cRat%s", namenew), "Ratios", dx, dy);
	cRati->Divide(nx, ny);
	
	TH1D*** hproj = new TH1D**[n2dh];
	Bool_t drawfirst = kTRUE;
	TH1D*** hRati = new TH1D**[n2dh];
	
	for(Int_t ih = 0; ih<n2dh; ih++){
		if(!hist[ih]){
			Printf("2D histo %d not found", ih);
			continue;
		}
		hproj[ih] = new TH1D*[nbinsrange];
		hRati[ih] = new TH1D*[nbinsrange];
		Printf("Color %d", cls[ih]);
	
		for(Int_t irb = 0; irb < nbinsrange; irb++){
			hproj[ih][irb] = 0x0;
			hRati[ih][irb] = 0x0;
			TString currentname = namenew;
			if(nametype == 0) currentname = Form("%s%d_R%d", namenew, ih,  irb);
			if(nametype == 1) currentname = Form("%s%d_R%.0f_%.0f", namenew, ih, rangelims[irb], rangelims[irb+1]);
			
			if(axrange == 0 && axpj == 1){ //projection on Y axis
				Int_t bin[2] = {hist[ih]->GetXaxis()->FindBin(rangelims[irb]), hist[ih]->GetXaxis()->FindBin(rangelims[irb+1] - eps)};
				
				
				
				
				hproj[ih][irb] = hist[ih]->ProjectionY(currentname, bin[0], bin[1]);
				
			}
			
			if(axrange == 1 && axpj == 0){ //projection on X axis
				Int_t bin[2] = {hist[ih]->GetYaxis()->FindBin(rangelims[irb]), hist[ih]->GetYaxis()->FindBin(rangelims[irb+1] - eps)};
				
				
				hproj[ih][irb] = hist[ih]->ProjectionX(currentname, bin[0], bin[1]);
				
				
			}
			hproj[ih][irb]->SetLineColor(cls[ih]);
			hproj[ih][irb]->SetMarkerColor(hproj[ih][irb]->GetLineColor());
			hproj[ih][irb]->SetMarkerStyle(mks[ih]);
			//hproj[ih][irb]->SetLineColor(ih);
			//hproj[ih][irb]->SetMarkerColor(hproj[ih][irb]->GetLineColor());
			//hproj[ih][irb]->SetMarkerStyle(20);
			hproj[ih][irb]->SetLineWidth(2);
			Double_t integral = hproj[ih][irb]->Integral();
			if(integral > 0) hproj[ih][irb]->Scale(1./integral);
			
			// ratio
			if(!hRati[ih][irb]) hRati[ih][irb] = (TH1D*)hproj[ih][irb]->Clone(Form("hRatio%s", currentname.Data()));
			
			
			cProj->cd(irb+1);
			if(drawfirst){
				hproj[ih][irb]->Draw();
				drawfirst = kFALSE;
			}
			else hproj[ih][irb]->Draw("sames");
		}
	}
	
	SaveCv(cProj);
	
	if(doratio){
	for(Int_t irb = 0; irb < nbinsrange; irb++){
		TH1 *h1bis;
		TH1 *h2bis;
		
		drawfirst=kTRUE;
		for(Int_t ih = 1; ih<n2dh; ih++){
			Int_t out = UniformTH1FForDivide(hRati[ih][irb], hproj[0][irb], h1bis, h2bis, "TH1D", kFALSE);
			
			hRati[ih][irb] = (TH1D*)h1bis;
			hproj[0][irb]  = (TH1D*)h2bis;
			
			hRati[ih][irb]->Divide(hproj[0][irb]);
			
			cRati->cd(irb+1);
			if(drawfirst){
				hRati[ih][irb]->GetYaxis()->SetRangeUser(0, 3.);
				hRati[ih][irb]->Draw();
				drawfirst = kFALSE;
			}
			else hRati[ih][irb]->Draw("sames");
		}
	}
	SaveCv(cRati);
	}
	
	
	return hproj;
}

#endif
