#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TClass.h>
#include <TGraphAsymmErrors.h>
#include <TPaveText.h>


TPaveText** CommonTools::GetPavePtBins(Int_t nptbins, Double_t *lims, Int_t x1 = 0.3, Int_t y1 = 0.8, Int_t x2 = 0.8, Int_t y2 = 0.9);

void SaveCv(TCanvas* c,TString suffix, Int_t format){
   if(format > 0) c->SaveAs(Form("%s%s.png",c->GetName(),suffix.Data()));
   if(format > 1) c->SaveAs(Form("%s%s.pdf",c->GetName(),suffix.Data()));
   if(format > 2) c->SaveAs(Form("%s%s.eps",c->GetName(),suffix.Data()));
   if(format > 3 || format==0) c->SaveAs(Form("%s%s.root",c->GetName(),suffix.Data()));
}


TList* CommonTools::ReadFile(TString strIn, TString strLst){
   
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
   return lst;
}

TList* CommonTools::ReadFile(TString strIn, TString strDir, TString strLst){
   
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
   return lst;
}


TObject* CommonTools::ReadObjInFile(TString strIn, TString strLst, TString objname){
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
//-------------------------------------------------------------------------------------------

TH1* CommonTools::CompareDistributions(TVirtualPad* pad, TH1* h1, TH1* h2){
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

void CommonTools::CalculatePads(Int_t n, Int_t&nx, Int_t&ny, Int_t&dx, Int_t&dy, Int_t perrow, Int_t stdd){
   
   Int_t percolumn = n/perrow;
   if(n%perrow > 0) percolumn++;
   nx = percolumn;
   ny = perrow;
   dx = stdd*nx;
   dy= stdd*ny;
   
   return;
}

void CommonTools::SaveGraphInCanvasToFile(TCanvas *corigin, TString graphname="", TString newname = "", TString fileoutname=""){
   
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

void CommonTools::UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type){
   if(!h1 || !h2){
      Printf("Need valid input!");
      return;
   }
   
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h2->GetNbinsX();
   Double_t binW1 = h1->GetBinWidth(3), binW2 = h2->GetBinWidth(3);
   Double_t low1 = h1->GetBinLowEdge(1), low2 = h2->GetBinLowEdge(1);
   Double_t up1 = h1->GetBinLowEdge(nBins1+1), up2 = h2->GetBinLowEdge(nBins2+1);
   Printf("1) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("2) %s N bins = %d, Width = %f, Range = %f - %f", h2->GetName(), nBins2, binW2, low2, up2);
   if((nBins1 == nBins2) && ((binW1-binW2) < 1e-6) &&  ((low1-low2) < 1e-6) && ((up1-up2) < 1e-6)) {
      Printf("Histograms can be divided");
      h1bis = h1;
      h2bis = h2;
      return;
   }
   Double_t commonlow = TMath::Min(low1, low2), commonup = TMath::Max(up1, up2);
   Double_t commonW = TMath::Max(binW1, binW2);
   Int_t commonNBins = (commonup - commonlow)/commonW;
   Printf("Common) N bins = %d, Width = %f, Range = %f - %f", commonNBins, commonW, commonlow, commonup);
   Double_t minWidth = TMath::Min(binW1, binW2);
   
   Int_t rebin = (Int_t)(commonW/minWidth);
   Printf("%f/%f = %f or %d", commonW, minWidth, commonW/minWidth, rebin);
   
   if(minWidth == binW1) h1->Rebin(rebin);
   else h2->Rebin(rebin);
   
   Printf("N bins : %s: %d , for %s: %d ", h1->GetName(), h1->GetNbinsX(), h2->GetName(), h2->GetNbinsX());
   if(type == "TH1F"){
      h1bis = new TH1F(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetXaxis()->GetTitle()), commonNBins, commonlow, commonup);
      h2bis = new TH1F(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetXaxis()->GetTitle()), commonNBins, commonlow, commonup);
   } else{
      if(type == "TH1D"){
      	 h1bis = new TH1D(Form("%sb", h1->GetName()), Form("%s bis; %s; %s", h1->GetTitle(), h1->GetXaxis()->GetTitle(), h1->GetXaxis()->GetTitle()), commonNBins, commonlow, commonup);
      	 h2bis = new TH1D(Form("%sb", h2->GetName()), Form("%s bis; %s; %s", h2->GetTitle(), h2->GetXaxis()->GetTitle(), h2->GetXaxis()->GetTitle()), commonNBins, commonlow, commonup);
      } else {
      	 Printf("%s Not defined!", type.Data());
      	 return;
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
      if(fillingbin>0 && fillingbin <= nBins1) {
      	 Printf("Here, debug %d", nBins1);
      	 fillWith = (Double_t)h1->GetBinContent(fillingbin);
      	 error = (Double_t)h1->GetBinError(fillingbin);
      }
      h1bis ->SetBinContent(i+1, fillWith);
      h1bis ->SetBinError(i+1, error);
      Printf("Filling %s x = %.3f new bin %d with old bin %d content %e", h1->GetName(), bincentre, i+1, fillingbin, fillWith);
      
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
      Printf("Filling %s x = %.3f new bin %d with old bin %d content %e", h2->GetName(), bincentre, i+1, fillingbin, fillWith);
   }
   
   
}

Int_t CommonTools::FindBinInArray(Double_t q, Double_t array[], Int_t dim, Bool_t verbose){
   
   if((q > array[dim-1]) || (q<array[0])) {
      if(verbose) Printf("Out of bound");
      return -1;
   }
   for(Int_t i = 0; i<dim-1; i++ ){
      if(q > array[i] && q < array[i+1]) return i;
      
   }
   return -2;
}


TH1F* CommonTools::ConvertTH1DinF(TH1D *h1d){
   //deletes the original histogram
   Int_t nbins = h1d->GetNbinsX();
   TH1F *h1f = new TH1F(Form("%s", h1d->GetName()), Form("%s;%s;%s", h1d->GetTitle(), h1d->GetXaxis()->GetTitle(), h1d->GetYaxis()->GetTitle()), nbins, h1d->GetBinLowEdge(1), h1d->GetBinLowEdge(nbins+1));
   
   h1f->Sumw2();
   
   
   for(Int_t i = 0; i<nbins; i++){
      h1f->SetBinContent(i+1, h1d->GetBinContent(i+1));
      h1f->SetBinError(i+1, h1d->GetBinError(i+1));
   }
   
   //delete h1d;
   //h1d = 0;
   return h1f;
   
}

TH2F* CommonTools::TransformAxisRanges(TH2F *h2orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew, Int_t ynbinsnew, Float_t yminnew, Float_t ymaxnew){
   
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

TH1F* CommonTools::TransformAxisRanges(TH1F *h1orig, Int_t xnbinsnew, Float_t xminnew, Float_t xmaxnew){
   
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

void CommonTools::DefineNewAxes(TH2 *hOrig, TH2 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Float_t deltaDownY, Float_t deltaUpY, Int_t& newNbinsX, Int_t& newNbinsY, Float_t& minX, Float_t& maxX, Float_t& minY, Float_t& maxY ){
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

void CommonTools::DefineNewAxes(TH1 *hOrig, TH1 *hdelta, Float_t deltaDownX, Float_t deltaUpX, Int_t& newNbinsX, Float_t& minX, Float_t& maxX ){
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

void CommonTools::PrintHisto2(TH2* h1){
   Int_t nBins1 = h1->GetNbinsX(), nBins2 = h1->GetNbinsY();
   Double_t binW1 = h1->GetXaxis()->GetBinWidth(3), binW2 = h1->GetYaxis()->GetBinWidth(3);
   Double_t low1 = h1->GetXaxis()->GetBinLowEdge(1), low2 = h1->GetYaxis()->GetBinLowEdge(1);
   Double_t up1 = h1->GetXaxis()->GetBinLowEdge(nBins1+1), up2 = h1->GetYaxis()->GetBinLowEdge(nBins2+1);
   Printf("X) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins1, binW1, low1, up1);
   Printf("Y) %s N bins = %d, Width = %f, Range = %f - %f", h1->GetName(), nBins2, binW2, low2, up2);
   
   return;

}

TPaveText** CommonTools::GetPavePtBins(const Int_t nptbins, Double_t *lims, Int_t x1, Int_t y1, Int_t x2, Int_t y2){
   
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
