#include <TF1.h>
#include <TList.h>
#include <THnSparse.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TKey.h>
#include <TLegend.h>
//#include </data/macros/LoadALICEFigures.C>
//#include </data/Work/MyCodeJetMass/utils/CommonTools.C>
#include </data/Work/MyCodeJetMass/macros/Convolution.C>

void Draw2DProjectionsinBins(THnSparse *hsp, Int_t ax2Dproj[], TString title, Int_t axSel, TString axSeltitle, Int_t npTbins, Double_t pTlims[]);
void Draw1DProjectionsinBins(THnSparse *hsp, Int_t ax1Dproj, TString title, Int_t axSel, TString axSeltitle, Int_t npTbins, Double_t pTlims[], Bool_t logy=kTRUE, Bool_t debug=kFALSE);
void ResetAxisSel(THnSparse *hps);
void plotOverlap(Int_t nfiles, TString paths[], TString listnames[], TString names[]);
void CompareRandomConeAndSingleTrackEmbedding(TString filename, TString dirnameRC, TString dirnameConst);
void Compare1DProjectionsinBins(THnSparse *hspB, THnSparse *hspR, Int_t ax1Dproj, TString title, Int_t axSel, TString axSeltitle, Int_t npTbins, Double_t pTlims[], Bool_t logy=kTRUE, Bool_t debug=kFALSE, Double_t normB= 1., Double_t normR = 1., TString legB = "", TString legR = "", Bool_t makeRatio = kFALSE);

void DrawCompare1DProjectionsinBins(Int_t nfiles, TString pathfiles[], TString legT[], TString listname, TString hspname, Int_t ax1Dproj, Int_t axSel, Int_t npTbins, Double_t pTlims[], Int_t norm = 2 /*0 = don't scale the histograms , 1 = scale by the number of events, 2 = scale to integral 1*/, Int_t rebinFactor = 1, Bool_t logy = kFALSE, Bool_t makeRatio = kTRUE);



Int_t saveLv=3; //2 standard, 3 also eps

void plotOutputTaskJetShapeDeriv(TString filename = "AnalysisResults.root", TString listname = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", Bool_t buseths = kTRUE){
   //
   // draws the mass distribution of jets matched with a single track embedded
   // in different pT bins
   //
   // buseths=kTRUE to be used when embedding of massive track because the 3D histogram contains mass, not delta M
   
   TList *list = 0x0;
   TH1F *hMassSub = 0x0;
   THnSparse *hnsp = 0x0;
   TH3F *h3MSubPtRawDRMatch = 0x0;
   Int_t axisP = 4; //this is for the normal input, namely the JetShapeTask
   if(!listname.IsNull()){  
      list = ReadFile(filename, listname);
      
      if(!list) return;
      TString h2name = "fh2MSubPtRawAll_0";
      TH2F *h2MSubPtRaw = (TH2F*)list->FindObject(h2name);
      TCanvas *cEventLook = new TCanvas("cEventLook", "pPb + embedded tracks", 800, 800);
      cEventLook->cd();
      gPad->SetLogz();
      h2MSubPtRaw->Draw("colz");
      SaveCv(cEventLook, "", saveLv);
      
      
      TString h3name = "fh3MSubPtRawDRMatch_0";
      
      h3MSubPtRawDRMatch = (TH3F*)list->FindObject(h3name);
      TH1F *hDR = (TH1F*)h3MSubPtRawDRMatch->Project3D("z");
      hMassSub = (TH1F*)h3MSubPtRawDRMatch->Project3D("x");
      Printf("Bin mass width %.2f", h3MSubPtRawDRMatch->GetXaxis()->GetBinWidth(3));
      
      TH2F* hMPt = (TH2F*)h3MSubPtRawDRMatch->Project3D("yx");
      TCanvas *cMpTMatched = new TCanvas("cMpTMatched", "Mass Vs pT, jets matched with embedded track", 800, 800);
      hMPt->Draw("colz");
      SaveCv(cMpTMatched, "", saveLv);
      
      TCanvas *cDR = new TCanvas("cDR", "Distance embedded track jet", 800, 800);
      cDR->cd();
      hDR->Draw("sames");
   }
   
   if(buseths) {
      if(list) hnsp = dynamic_cast<THnSparse*>(list->FindObject("fhnDeltaMass_0"));
      if(!hnsp) Printf("Something wrong, fhnDeltaMass_0 not found");
      hMassSub = (TH1F*)hnsp->Projection(0);
   }
   hMassSub->SetName("hMassSub");
   hMassSub->SetLineWidth(2);
   
   TCanvas * cmass= new TCanvas("cmass", "Mass subtracted", 800, 800);
   cmass->cd();
   hMassSub->Scale(1./hMassSub->Integral());
   hMassSub->Draw();
   TPaveText *pvmean = new TPaveText(0.1, 0.8, 0.5, 0.9, "NDC");
   pvmean->SetFillStyle(0);
   pvmean->SetBorderSize(0);
   pvmean->SetTextColor(hMassSub->GetLineColor());
   pvmean->AddText(Form("Mean = %.3f #pm %.0e; #sigma = %.2f", hMassSub->GetMean(), hMassSub->GetMeanError(), hMassSub->GetRMS()));
   cmass->cd();
   pvmean->Draw();
   TString namefout = "MassFromBkgFluctuations.root";
   TFile *fout = new TFile(namefout, "recreate");
   Printf("This macro writes an output in a ROOT file (%s).", namefout.Data());
   gDirectory->cd(0);
   
   Double_t pTlims[4] = {40., 60., 80., 100.};
   for(Int_t ipt = 0; ipt < 3; ipt++){
      Int_t biny1= 1;
      Int_t biny2= 1;
      TH1F *hMassSubpTSel = 0x0;
      if(h3MSubPtRawDRMatch){
      	 biny1 = h3MSubPtRawDRMatch->GetYaxis()->FindBin(pTlims[ipt]);
      	 biny2 = h3MSubPtRawDRMatch->GetYaxis()->FindBin(pTlims[ipt+1]) -1;
      
      	 hMassSubpTSel = (TH1F*)h3MSubPtRawDRMatch->ProjectionX(Form("hMassSubpT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]), biny1, biny2, 0, -1);
      }
      
      if(buseths && hnsp) {
      	 biny1 = hnsp->GetAxis(axisP)->FindBin(pTlims[ipt]); //pTdet level
      	 biny2 = hnsp->GetAxis(axisP)->FindBin(pTlims[ipt+1]) -1;
      	 hnsp->GetAxis(axisP)->SetRange(biny1, biny2);
      	 //pTdet level -> currently it's flat distr
      	 hMassSubpTSel = (TH1F*)hnsp->Projection(0); //projection on dM
      	 hMassSubpTSel->SetName(Form("hMassSubpT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]));
      }
      
      hMassSubpTSel->SetLineColor(colors[ipt]);
      cmass->cd();
      hMassSubpTSel->Scale(1./hMassSubpTSel->Integral());
      hMassSubpTSel->Draw("sames");
      
      TPaveText *pvmeanpT = new TPaveText(0.1, 0.4+ (0.1*ipt), 0.5, 0.4+(0.1*(ipt+1)), "NDC");
      pvmeanpT->SetFillStyle(0);
      pvmeanpT->SetBorderSize(0);
      pvmeanpT->SetTextColor(hMassSubpTSel->GetLineColor());
      pvmeanpT->AddText(Form("%.0f-%.0f GeV/#it{c}",pTlims[ipt], pTlims[ipt+1]));
      pvmeanpT->AddText(Form("Mean = %.3f #pm %.0e; #sigma = %.2f", hMassSubpTSel->GetMean(), hMassSubpTSel->GetMeanError(), hMassSubpTSel->GetRMS()));
      cmass->cd();
      pvmeanpT->Draw();
      fout->cd();
      hMassSubpTSel->Write();
   }
   SaveCv(cmass, "", saveLv);
   /*
   THnSparseF *hsp = (THnSparseF*) list->FindObject("fhnMassResponse_0");
   Int_t npTBins = hsp->GetAxis(3)->GetNbins();
   TH2F *h2pT = hsp->Projection(3,2);
   
   for(Int_t iX = 0; iX < npTBins; iX++){
     	for(Int_t iY = 0; iY < npTBins; iY++){
     	   
     	   Double_t deltapT = h2pT->GetXaxis()->GetBinCenter(iX+1) -  h2pT->GetYaxis()->GetBinCenter(iY+1);
     	   Double_t binContent = h2pT->GetBinContent(iX+1,iY+1);
     	   hDeltapT->Fill(deltapT,binContent);
     	   
     	  
     	   
     	}
     }
     */
}

void PerformConvolution(TString inputFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root",  TString listname = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged", TString correctionFile = "MassFromBkgFluctuations.root", Int_t ndraw = 0, Bool_t xpos = kTRUE){
   
   //
   // Convolute the mass distribution obtained above with a MC particle level mass distribution
   // The aim is to see the effect of the background fluctuations, encoded in the mass of jets matched with a high-pT embedded single track, on the original mass distribution
   //
   Bool_t analytic = kTRUE;
   if (ndraw > 0) {
      analytic = kFALSE;
      Printf("Use MC method with %d iterations", ndraw);  
   }
   
   //SetStyle(0);
   TList *l=ReadFile(inputFile, listname);

   TString hname = "fhnMassResponse";
   
   THnSparseF *hsp = (THnSparseF*) l->FindObject(hname);
   if(!hsp){
      Printf("THnSparse not found");
      return;
   }
   TFile *f = new TFile(correctionFile);
   TString namefout = Form("ModifiedMassBkgFluct%s.root", analytic ? "" : "Rndm");
   TFile *fout = new TFile(namefout, "recreate");
   Printf("This macro writes an output in a ROOT file (%s).", namefout.Data());
   Double_t pTlims[4] = {40., 60., 80., 100.};
   
   Int_t nx, ny, dx, dy;
   CalculatePads(4, nx, ny, dx, dy);
   TCanvas *cR = new TCanvas("cR", "Ratios", dx, dy);
   cR->Divide(nx, ny);
   TCanvas *cdMdpt = new TCanvas("cdMdpt", "Delta M delta pT used as input", dx, dy);
   cdMdpt->Divide(nx, ny);
   
   for(Int_t ipt = 0; ipt < 3; ipt++){
      TH1F *h = (TH1F*)f->Get(Form("hMassSubpT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]));
      if(!h){
      	 Printf("%d) not found", ipt);
      	 continue;
      }
      cdMdpt->cd(ipt+1);
      h->Draw();
      
      Int_t binpt[2] = {hsp->GetAxis(3)->FindBin(pTlims[ipt]) , hsp->GetAxis(3)->FindBin(pTlims[ipt+1])-1 }; //3 = pt particle level
      hsp->GetAxis(3)->SetRange(binpt[0], binpt[1]);
      
      TH1F* hMass = (TH1F*)hsp->Projection(1); //1 = mass particle level
      hMass->SetName(Form("hMasspT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]));
      Printf("Projection %d", ipt);
      hMass->SetMarkerStyle(20);
      Bool_t drawlog = kTRUE;
      TH1F* hMassNew = 0x0;
      TCanvas *cMNew = 0x0; 
      if(analytic) hMassNew = ConvolutionH1XH2(hMass, h, drawlog);
      else {
      	 Float_t deltaDownX = -15, deltaUpX = 20;
      	 Int_t newNbinsX;
      	 Float_t minX, maxX;
      	 DefineNewAxes(hMass, h, deltaDownX, deltaUpX, newNbinsX, minX, maxX );
      	 TH1F* hMassBis = TransformAxisRanges(hMass, newNbinsX, minX, maxX);
      	 hMassBis->SetName(hMass->GetName());
      	 hMassNew = Convolution1DRandom(hMassBis, h, ndraw, xpos);
      	 if(!hMassNew) {
      	    Printf("Convolution result is null");
      	    return;  
      	 }
      	 hMassNew->SetMarkerStyle(24);
      	 hMassNew->SetMarkerColor(kBlue);
      	 cMNew = new TCanvas(Form("cMNew%d", ipt), Form("Mass Folded pT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]), 600, 600);
      	 cMNew->cd();
      	 if(drawlog) gPad->SetLogy();
      	 hMassNew->Draw("E");
      	 hMass->Draw("Esames");
      	 
      }
      TPaveText *pvInfo = new TPaveText(0.3, 0.2, 0.5, 0.4, "NDC");
      pvInfo->SetFillStyle(0);
      pvInfo->SetBorderSize(0);
      Double_t muorig = hMass->GetMean(), sigmaorig = hMass->GetRMS();
      pvInfo->AddText(Form("#delta#mu %.1f%s, #delta#sigma %.1f%s", (hMassNew->GetMean()-muorig)/muorig *100., "%", (hMassNew->GetRMS()-sigmaorig)/sigmaorig *100., "%"));
      
      //ratios
      TH1* hratioMass2 = 0x0;
      if(hMass && hMassNew){ 
      	 TH1F *hratioMass = (TH1F*)(hMass->Clone(Form("hMassRatiopT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1])));
      	 
      	 hratioMass->Scale(hMassNew->Integral()/hMass->Integral());
      	 gPad->cd();
      	 hratioMass->SetMarkerStyle(24);
      	 hratioMass->DrawClone("sames");
      	 
      	 pvInfo->Draw();
      	 TH1 *hMassNew2 = 0x0;
      	 
      	 UniformTH1FForDivide(hratioMass, hMassNew, hratioMass2, hMassNew2);
      	 hratioMass2->SetTitle("Ratio Mass Part/Convoluted; #it{M} (GeV); Ratio Part/Convoluted");
      	 hratioMass2->SetName(hratioMass->GetName());  	 
      	 hratioMass2->Divide(hMassNew2);
      	 hratioMass2->GetXaxis()->SetRangeUser(-0.5, 20);
      	 cR->cd(ipt+1);
      	 hratioMass2->Draw();
      }
      if(cMNew) SaveCv(cMNew);
      fout->cd();
      hMassNew->Write();
      if(hratioMass2) hratioMass2->Write();
      
   }
   SaveCv(cR);
}

void CompareFluctResults(TString file1 = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train976/analysis/Convolution/ModifiedMassBkgFluct.root", TString leg1 = "Analytic", TString file2 = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train976/analysis/Convolution/MCMethodm/ModifiedMassBkgFluctRndm.root", TString leg2 = "MCmethod", Int_t type = 1){

   //Type = 1 compare PerformConvolution output
   //Type = 2 compare PerformConvolution2D output: note the order is important, first 1D, then 2D
   
   Printf("Drawing files \n%s\n%s", file1.Data(), file2.Data());
   TFile *file[2];
   file[0] = new TFile(file1);
   file[1] = new TFile(file2);
   TLegend *leg = new TLegend(0.4, 0.3, 0.7, 0.5);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   TString legtext[2] = {leg1, leg2};
   
   TString *hbasename = 0x0;
   TString basesuff = "";
   Int_t   nhist = 1;
   
   if(type == 1){ //
      nhist = 2; // when changing this, need to modify the code 
      hbasename = new TString[nhist];
      hbasename[0] = "hMasspT";
      hbasename[1] = "hMassRatiopT" ; //TH1F , TH1D
      basesuff = "Modif"; //to be added to the first
   }
   if(type == 2){
      nhist = 1; // when changing this, need to modify the code
      hbasename = new TString[2];
      hbasename[0] = "hMasspT";
      hbasename[1] = "hMNewpTNew_";
      basesuff = "Modif";
   }
   const Int_t nhistFixed = nhist;
   
   const Int_t npt = 3;
   Double_t pTlims[npt+1] = {40., 60., 80., 100.};
   Int_t nx, ny, dx, dy;
   CalculatePads(npt, nx, ny, dx, dy);
   
   //create the needed canvas
   TCanvas *c[nhistFixed];
   for(Int_t ic = 0; ic<nhistFixed; ic++){
      c[ic] = new TCanvas(Form("c%s", hbasename[ic].Data()), Form("Canvas %s", hbasename[ic].Data()), dx, dy);
      c[ic]->Divide(nx,ny);
      
   }
   
   //loop over the 2 files to be compared
   for(Int_t ifl = 0; ifl < 2; ifl++){
      if(!file[ifl]->IsOpen()) {
      	 Printf("File %d not found", ifl);
      	 continue;
      }
      
      //loop over pt bins
      for(Int_t ipt = 0; ipt<npt; ipt++){
      	 
      	 Printf("Bin %d", ipt);
      	 
      	 TString hfullname = Form("%s%.0f-%.0f%s", hbasename[0].Data(), pTlims[ipt], pTlims[ipt+1], basesuff.Data());
      	 if(type == 2) hfullname = Form("%s%.0f%s%.0f%s", hbasename[ifl].Data(), pTlims[ipt], ifl==0 ? "-" : "", pTlims[ipt+1], ifl==0 ? basesuff.Data() : "");

      	 TH1 *hMass = 0x0;
      	 hMass = ((TH1F*)file[ifl]->Get(hfullname));
      	 //if(type == 2 ) hMass = dynamic_cast<TH1D*>(file[ifl]->Get(hfullname));
      	 if(hMass){
      	    hMass->Scale(1./hMass->Integral());
      	    hMass->SetMarkerColor(colors[ifl]);
      	    c[0]->cd(ipt+1);
      	    gPad->SetLogy();
      	    if(ifl==0) hMass->Draw();
      	    else hMass->Draw("sames");
      	 } else {
      	    Printf("%s not found in file %d", hfullname.Data(), ifl);
      	    file[ifl]->ls();
      	 }
      	 if(ipt == 0) {
      	    leg->AddEntry(hMass, legtext[ifl], "P");
      	 }
      	 if(type == 2) continue; //No ratio available for the time being
      	 TH1F *hMassRatio = dynamic_cast<TH1F*>(file[ifl]->Get(Form("%s%.0f-%.0f", hbasename[1].Data(), pTlims[ipt], pTlims[ipt+1])));
      	 if(!hMassRatio){
      	    hMassRatio = (TH1F*)file[ifl]->Get(Form("%s%.0f-%.0f", hbasename[1].Data(), pTlims[ipt], pTlims[ipt+1]));
      	 }
      	 if(hMassRatio){
      	    hMassRatio->Scale(1./hMassRatio->Integral());
      	    hMassRatio->SetMarkerColor(colors[ifl]);
      	    c[1]->cd(ipt+1);
      	    if(ifl==0) hMassRatio->Draw();
      	    else hMassRatio->Draw("sames");
      	 }  else Printf("hMassRatio not found in file %d", ifl);
      	 
      }
   }
   for(Int_t ic = 0; ic<nhistFixed; ic++){
      c[ic]->cd(npt);
      leg->Draw();
      SaveCv(c[ic]);
   }
}

void PerformConvolution2D(TString inputFilePythia = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root",  TString listname = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged", TString inputFileFluct = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train976/output/runlist1/AnalysisResults.root", TString listnameFluc = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", Int_t ndraw = 1000, Bool_t xpos = kTRUE, Bool_t ypos = kTRUE){
   
   //
   // Convolute the mass, pt distribution obtained above with a MC particle level mass distribution
   // The aim is to see the effect of the background fluctuations, encoded in the mass of jets matched with a high-pT embedded single track, on the original mass distribution
   //
   // The dMdpT input is not yet final, this has flat pT embedded distribution
   
   //SetStyle(0);
   //read PYTHIA, particle level to be convoluted
   
   
   const Int_t nptbins = 3;
   Double_t pTlims[nptbins+1] = {40., 60., 80., 100.};
   TString hname = "";
   THnSparseF *hsp = 0x0;
   if(!listname.IsNull()){
      TList *lP = ReadFile(inputFilePythia, listname);
      if(!lP) return;
      hname = "fhnMassResponse";
      hsp = (THnSparseF*) lP->FindObject(hname);
      if(!hsp){
      	 Printf("THnSparse not found");
      	 return;
      }
   } else{
      TFile *fin = new TFile(inputFilePythia);
      if(!fin->IsOpen()) return;
      hsp = dynamic_cast<THnSparseF*>(fin->Get("hShapeRespTree"));
      hname = "hMassRespTree";
      //axisP = 2; //see DrawTreesOutputTask
      if(!hsp){
      	 Printf("THnSparse not found");
      	 return;
      }
      
      
   }
   
   Int_t binpt[2] = {hsp->GetAxis(3)->FindBin(pTlims[0]) , hsp->GetAxis(3)->FindBin(pTlims[nptbins])-1 }; //3 = pt particle level
   hsp->GetAxis(3)->SetRange(binpt[0], binpt[1]);
   //used to sum up the different pT bins below
   TH2F* hMasspTAllPtOrig = (TH2F*)hsp->Projection(1, 3);
   hMasspTAllPtOrig->Reset();
   hMasspTAllPtOrig->SetName(Form("hMasspT%.0f-%.0f", pTlims[0], pTlims[nptbins]));
   
         
   TH2F* hMasspTAllPtNew = 0x0;
   
   //read fluctuations -> dMdpT
   TList *lF = ReadFile(inputFileFluct, listnameFluc);
   if(!lF) return;
   TString hnameF = "fhnDeltaMass_0";
   THnSparseF *hspF = (THnSparseF*) lF->FindObject(hnameF);
   if(!hspF){
      Printf("THnSparse Fluc not found");
      return;
   }
   Printf("Nentries %.0f", hspF->GetEntries());
   TString namefout = "ModifiedMassBgkFluct.root";
   TFile *fout = new TFile(namefout, "recreate");
   Printf("This macro writes an output in a ROOT file (%s).", namefout.Data());
   
   TH2F *hdMdpTAll = (TH2F*) hspF->Projection(0,1);
   Int_t nbinX, nbinY;
   Float_t minValX, minValY, maxValX, maxValY;
   //x = Delta pT, y = delta M
   // these number are set by hand looking ad the dM, dpT distribution
   Float_t deltaUpX = 50., deltaUpY = 15, deltaDownX = -50, deltaDownY = -15;
   DefineNewAxes(hMasspTAllPtOrig, hdMdpTAll, deltaDownX, deltaUpX, deltaDownY, deltaUpY, nbinX, nbinY, minValX, maxValX, minValY, maxValY );
   
   TH2F *hMasspTAllPt =  TransformAxisRanges(hMasspTAllPtOrig, nbinX,  minValX, maxValX, nbinY,  minValY, maxValY);
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cR = new TCanvas("cR", "Ratios", dx, dy);
   cR->Divide(nx, ny);
   TCanvas *cdMdpt = new TCanvas("cdMdpt", "Delta M delta pT used as input", dx, dy);
   cdMdpt->Divide(nx, ny);
   TCanvas *cM2D = new TCanvas("cM2D", Form("Mass Particle Level Input"), dx, dy);
   cM2D->Divide(nx, ny);
   
   Int_t axisptP = 5;
   Int_t dims = hspF->GetNdimensions();
   if(axisptP >= dims) Printf("Error, axis %d out of bound (max %d)", axisptP, dims);
   TH2F *hdMdpT[nptbins];
   
   Bool_t drawlog = kTRUE;
   TH2F* hMasspTNew = 0x0;
   TCanvas *cM2DNew = cM2DNew = new TCanvas("cM2DNew", "Mass Folded in pT bins", dx, dy);
   cM2DNew->Divide(nx, ny);
   
   for(Int_t ipt = 0; ipt < nptbins; ipt++){
      //prepare fluctuations
      Int_t binF[2] = { hspF->GetAxis(axisptP)->FindBin(pTlims[ipt]), hspF->GetAxis(axisptP)->FindBin(pTlims[ipt+1]-1)};
      
      hspF->GetAxis(axisptP)->SetRange(binF[0], binF[1]);
      hdMdpT[ipt] = (TH2F*) hspF->Projection(0,1); //0 = dM, 1 = dpT
      if(!hdMdpT[ipt]){
      	 Printf("%d) not found", ipt);
      	 continue;
      }
      hdMdpT[ipt]->SetName(Form("hdMdpT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]));
      Printf("Now entries %.0f bins %d, %d",  hdMdpT[ipt]->GetEntries(), binF[0], binF[1]);
      cdMdpt->cd(ipt+1);
      hdMdpT[ipt]->Draw("colz");
   

     
      //prepare Pythia particle level
      Int_t binpt[2] = {hsp->GetAxis(3)->FindBin(pTlims[ipt]) , hsp->GetAxis(3)->FindBin(pTlims[ipt+1])-1 }; //3 = pt particle level
      hsp->GetAxis(3)->SetRange(binpt[0], binpt[1]);
      
      TH2F* hMasspTOrig = (TH2F*)hsp->Projection(1, 3); //1 = mass particle level, 3 = particle level pT
      //hMasspTOrig->SetName(Form("hMasspT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]));
      hMasspTOrig->SetName(Form("hMasspT%.0f-%.0f", pTlims[ipt], pTlims[ipt + 1]));
            
      
      
      Printf("Transforming");
      PrintHisto2(hMasspTOrig);
      Printf("Into \nX) N bins %d, min %.3f, max %.3f", nbinX,  minValX, maxValX);
      Printf("Y) N bins %d, min %.3f, max %.3f", nbinY,  minValY, maxValY);
      
      TH2F* hMasspT = TransformAxisRanges(hMasspTOrig, nbinX,  minValX, maxValX, nbinY,  minValY, maxValY);
      //TH2F* hMasspT = TransformAxisRanges(hMasspTOrig, nbinX,  binmaxX, binmaxX, nbinY,  binmaxY, binmaxY);

      cM2D->cd(ipt+1);
      hMasspT->Draw("colz");
      
      
      hMasspTAllPt->Add(hMasspT);
      
      //here comes the convolution
      
      hMasspTNew = Convolution2DRandom(hMasspT, hdMdpT[ipt], ndraw);
      if(!hMasspTNew) {
      	 Printf("Convolution result is null");
      	 return;  
      }
      hMasspTNew->GetYaxis()->SetTitle("#it{M} (GeV)");
      hMasspTNew->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      
      cM2DNew->cd(ipt+1);
      if(drawlog) gPad->SetLogz();
      hMasspTNew->Draw("colz");
      
      if(!hMasspTAllPtNew) hMasspTAllPtNew = dynamic_cast<TH2F*>(hMasspTNew->Clone(Form("hMasspTNew%.0f-%.0f", pTlims[0], pTlims[nptbins + 1])));
      else hMasspTAllPtNew->Add(hMasspTNew);
      
   }
   //save input and partial folded
   SaveCv(cdMdpt);
   SaveCv(cM2D);
   SaveCv(cM2DNew);
   
   TCanvas *cM2DAllPt = new TCanvas(Form("cM2DAllPt"), Form("Mass before and after pT%.0f-%.0f", pTlims[0], pTlims[nptbins]), 900, 500);
   cM2DAllPt->Divide(2, 1);
   cM2DAllPt->cd(1);
   hMasspTAllPt->DrawClone("colz");
   cM2DAllPt->cd(2);
   hMasspTAllPtNew->DrawClone("colz");
   SaveCv(cM2DAllPt);
   
   //Printf("Projection %d", ipt);
   TH1D* hMass = hMasspTAllPt->ProjectionY(Form("hMass%.0f-%.0f", pTlims[0], pTlims[nptbins]));
   hMass->SetMarkerStyle(20);
   hMass->Scale(1./hMass->Integral("width"));
   TH1D* hpT = hMasspTAllPt->ProjectionX(Form("hpT%.0f-%.0f", pTlims[0], pTlims[nptbins]));
   hpT->SetMarkerStyle(20);
   hpT->Scale(1./hpT->Integral("width"));
   
   //project onto the two axes and draw
   TH1D *hMassNew = hMasspTAllPtNew->ProjectionY(Form("hMass%.0f-%.0fNew", pTlims[0], pTlims[nptbins]));
   hMassNew->Scale(1./hMassNew->Integral("width"));
   hMassNew->SetMarkerStyle(24);
   
   TH1D *hpTNew = hMasspTAllPtNew->ProjectionX(Form("hpT%.0f-%.0fNew", pTlims[0], pTlims[nptbins]));
   hpTNew->Scale(1./hpTNew->Integral("width"));
   hpTNew->SetMarkerStyle(24);
   
   TCanvas *cMPjNew = new TCanvas(Form("cMpJNew%d", 0), Form("Mass Folded pT%.0f-%.0f", pTlims[0], pTlims[nptbins]), 1000, 500);
   cMPjNew->Divide(2,1);
   TLegend *legPj = new TLegend(0.5, 0.6, 0.9, 0.8);
   legPj->SetBorderSize(0);
   legPj->SetFillStyle(0);
   legPj->AddEntry(hMassNew, "Convoluted", "P");
   legPj->AddEntry(hMass, "Particle lev", "P");
   
   cMPjNew->cd(1);
   if(drawlog) gPad->SetLogy();
   hMassNew->Draw("E");
   hMass->Draw("Esames"); //
   legPj->Draw();
   
   cMPjNew->cd(2);
   if(drawlog) gPad->SetLogy();
   hpTNew->Draw("E");
   hpT->Draw("Esames");
   
   TPaveText *pvInfo = new TPaveText(0.2, 0.2, 0.6, 0.4, "NDC");
   pvInfo->SetFillStyle(0);
   pvInfo->SetBorderSize(0);
   Double_t muorig = hMass->GetMean(), sigmaorig = hMass->GetRMS();
   pvInfo->AddText(Form("#delta#mu %.1f%s, #delta#sigma %.1f%s", (hMassNew->GetMean()-muorig)/muorig *100., "%", (hMassNew->GetRMS()-sigmaorig)/sigmaorig *100., "%"));
   cMPjNew->cd(1);
   pvInfo->Draw();
   
   SaveCv(cMPjNew);
   
   TCanvas *cMpjNewpTbinNew = new TCanvas("cMpjNewpTbinNew", "Mass folded in bins of folded pT", dx, dy);
   cMpjNewpTbinNew->Divide(nx, ny);
   fout->cd();
   hMasspTAllPt->Write();
   hMasspTAllPtNew->Write();
   for(Int_t ipt=0; ipt<nptbins; ipt++){
      Int_t binptNew[2] = { hMasspTAllPtNew->GetXaxis()->FindBin(pTlims[ipt]), hMasspTAllPtNew->GetXaxis()->FindBin(pTlims[ipt+1]-1)};
      hMasspTAllPtNew->GetXaxis()->SetRange(binptNew[0], binptNew[1]);
      
      TH1D *hMNew = hMasspTAllPtNew->ProjectionY(Form("hMNewpTNew_%.0f%.0f", pTlims[ipt], pTlims[ipt+1]));
      hMNew->SetMarkerStyle(24);
      hMNew->Sumw2();
      
      cMpjNewpTbinNew->cd(ipt+1);
      if(drawlog) gPad->SetLogy();
      hMNew->Draw();
      fout->cd();
      hMNew->Write();
      
      TPaveText *pvptDet = new TPaveText(0.4, 0.7, 0.8, 0.85, "NDC");
      pvptDet->SetBorderSize(0);
      pvptDet->SetFillStyle(0);
      pvptDet->AddText(Form("%.0f < #it{p}_{T,det} < %.0f GeV/#it{c}", pTlims[ipt], pTlims[ipt+1]));
      cMpjNewpTbinNew->cd(ipt+1);
      pvptDet->Draw();
   }
   SaveCv(cMpjNewpTbinNew);
   
   //ratios (TO BE CHECKED)
   TH1* hratioMass2 = 0x0;
   /*
   if(hMass && hMassNew){ 
      TH1F *hratioMass = (TH1F*)(hMass->Clone(Form("hMassRatiopT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1])));
      
      //hratioMass->Scale(hMassNew->Integral()/hMass->Integral());
      cMPjNew->cd(ipt+1);
      hratioMass->SetMarkerStyle(25);
      hratioMass->DrawClone("sames");
      
      
      TH1 *hMassNew2 = 0x0;
      //hMassNew2 = hMassNew;
      UniformTH1FForDivide(hratioMass, hMassNew, hratioMass2, hMassNew2);
      hratioMass2->SetTitle("Ratio Mass Part/Convoluted; #it{M} (GeV); Ratio Part/Convoluted");  	 
      hratioMass2->Divide(hMassNew2);
      hratioMass2->GetXaxis()->SetRangeUser(-0.5, 20);
      cR->cd(ipt+1);
      hratioMass2->Draw();
      
   }
   
   fout->cd();
   hMassNew->Write();
   //if(hratioMass2) hratioMass2->Write();
   
   SaveCv(cR);
   */
}

void CompareFoldedWithEmbedding(TString inputFold = "ModifiedMassBgkFluct.root", TString inputEmb = "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train1061/output/AnalysisResults.root", TString listEmb = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracksEmb_pT0150_E_scheme_TC"){
   
   // the Particle-level folded distribution should be similar to the data embedded distribution, or better to the data only if both detector and fluctuations are considered
   
   TFile *finFold = new TFile(inputFold);
   if(!finFold->IsOpen()){
      Printf("%s not found", inputFold.Data());
      return;
   }
   TString hbaseFoldName = "hMNewpTNew_";
   
   TString hname = "fhnMassResponse_0";
   THnSparseF *hsp = (THnSparseF*)ReadObjInFile(inputEmb, listEmb, hname);
   if(!hsp) return;
   
   const Int_t nptbins = 3;
   Double_t pTlims[nptbins+1] = {40., 60., 80., 100.};
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   
   TCanvas *cMCmp = new TCanvas("cMCmp", "Mass from folding and from Embedding response in bins of folded/detector pT ", dx, dy);
   cMCmp->Divide(nx, ny);
   
   TLegend *leg = new TLegend(0.5, 0.3, 0.9, 0.5);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   
   for(Int_t ipt=0; ipt<nptbins; ipt++){
      
      //folded
      TH1D *hMFold = dynamic_cast<TH1D*>(finFold->Get(Form("%s%.0f%.0f", hbaseFoldName.Data(), pTlims[ipt], pTlims[ipt+1])));
      if(!hMFold){
      	 Printf("%s%.0f%.0f not found", hbaseFoldName.Data(), pTlims[ipt], pTlims[ipt+1]);
      	 finFold->ls();
      	 continue;
      }
      hMFold->Scale(1./hMFold->Integral());
      
      //embedding
      Int_t binpt[2] = {hsp->GetAxis(2)->FindBin(pTlims[ipt]) , hsp->GetAxis(2)->FindBin(pTlims[ipt+1])-1 }; //2 = pt detector level
      hsp->GetAxis(2)->SetRange(binpt[0], binpt[1]);
      
      TH1D* hMassEmb = hsp->Projection(0); //0 = mass detector level
      hMassEmb->SetName(Form("hMassEmb_%.0f%.0f", pTlims[ipt], pTlims[ipt+1]));
      hMassEmb->Sumw2();
      hMassEmb->Scale(1./hMassEmb->Integral());
      hMassEmb->SetMarkerStyle(20);
      hMassEmb->SetMarkerColor(colors[3]);
   
      if(ipt == 0) {
      	 leg->AddEntry(hMFold, "Particle Folded", "P");
      	 leg->AddEntry(hMassEmb, "Embedding", "P");
      }
      cMCmp->cd(ipt+1);
      gPad->SetLogy();
      hMassEmb->Draw("E");
      hMFold->Draw("sames");
      leg->Draw();
   }
   SaveCv(cMCmp);

}
//_______________________________________________________________________________________________
void WriteMassvspTDistribution(TString inputFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root", TString listname = "JetMass_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCMerged", TString hname = "fh3PtJet1VsMassVsLeadPtTagged_0", Bool_t isThnSparse = kFALSE, Int_t pax1 = 0, Int_t pax2 =0, TString suffix = ""){
   
   // save the Mass vs pT distribution to do the unfolding/refolding
   
   // list names 
   // "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged", "fhnMassResponse"
   
   SetStyle(0);
   TList *l=ReadFile(inputFile, listname);
   TH2D *hMasspT = 0x0;
   if(!l){
      Printf("List %s NOT found", listname.Data());
      return;
   }
   if(isThnSparse){
      THnSparseF *hsp = 0x0;
      hsp = (THnSparseF*) l->FindObject(hname);
      if(!hsp){
      	 Printf("THnSparse not found");
      	 l->ls();
      	 return;
      }
      hMasspT = (TH2D*) hsp->Projection(pax1, pax2);
      
   } else {
      TH3D *h3D = (TH3D*) l->FindObject(hname);
      if(!h3D){
      	 Printf("TH3D not found");
      	 return;
      }
      Int_t axesReco[2] = {};
      hMasspT = (TH2D*)h3D->Project3D("yx");
      hMasspT->SetName("hMasspT");
      
   }
   TFile *fout = new TFile(Form("MassvspT%s.root", suffix.Data()), "recreate");
   fout->cd();
   hMasspT->Write();
   Printf("Written file MassvspT%s.root", suffix.Data());

   TCanvas *cshow = new TCanvas("cshow", Form("Histogram saved %s", hMasspT->GetName()), 600, 600);
   hMasspT->Draw("colz");
}

void PlotRatioMassInpTbins(TString path1, TString suffix1, TString path2 = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576", TString suffix2 = ""){
   
   //first file
   TString filename = Form("%s/MassvspT%s.root", path1.Data(), suffix1.Data());
   TFile *fin = new TFile(filename.Data());
   if(!fin->IsOpen()){
      Printf("File %s not found", filename.Data());
      return;
   }
   
   TString histoname = fin->GetListOfKeys()->At(0)->GetName();
   Printf("%s", histoname.Data());
   TH2D *hpTMass = (TH2D*)fin->Get(histoname);
   if(!hpTMass) {
      Printf("%s not found", histoname.Data());
      return;
   
   }
   
   //second file (Reference)
   TString filenameR = Form("%s/MassvspT%s.root", path2.Data(), suffix2.Data());
   TFile *finR = new TFile(filenameR.Data());
   if(!finR->IsOpen()){
      Printf("File %s not found", filenameR.Data());
      return;
   }
   
   TString histonameR = finR->GetListOfKeys()->At(0)->GetName();
   Printf("%s", histonameR.Data());
   TH2D *hpTMassR = (TH2D*)finR->Get(histonameR);
   if(!hpTMassR) {
      Printf("%s not found", histonameR.Data());
      return;
   
   }
   
   TCanvas *c2D = new TCanvas("c2D", "Mass vs pT", 600, 600);
   hpTMass->Draw("colz");
   TCanvas *c2DR = new TCanvas("c2DR", "Mass vs pT (reference)", 600, 600);
   hpTMassR->Draw("colz");
   
   Int_t nptbins = 2;
   Double_t ptlims[nptbins+1] = {20., 40., 100.};
   
   //properties of histogram 1
   Int_t massXAxis = 0;
   TString titleX = hpTMass->GetXaxis()->GetTitle(), titleBins = "";
   if(titleX.Contains("M")) {
      massXAxis = 1;
      titleBins = hpTMass->GetYaxis()->GetTitle();
   }
   else {
      titleBins = titleX;
      Printf("X axis = %s, Y axis is mass", titleX.Data());
   }
   
   //properties of reference histogram
   Int_t massXAxisR = 0;
   TString titleXR = hpTMassR->GetXaxis()->GetTitle(), titleBinsR = "";
   if(titleXR.Contains("M")) {
      massXAxisR = 1;
      titleBinsR = hpTMassR->GetYaxis()->GetTitle();
   }
   else {
      titleBinsR = titleXR;
      Printf("Reference X axis = %s, Y axis is mass", titleXR.Data());
   }
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cMvsPt = new TCanvas("cMvsPt", "Mass in pT bins", dx, dy);
   cMvsPt->Divide(nx, ny);
   TCanvas *cMvsPtBis = new TCanvas("cMvsPtBis", "Rebinned and resized histograms", dx, dy);
   cMvsPtBis->Divide(nx, ny);
   TCanvas *cMvsPtRatios = new TCanvas("cMvsPtRatios", "Mass ratio with reference in pT bins", dx, dy);
   cMvsPtRatios->Divide(nx, ny);
   Printf("Loop on %d pT bins ", nptbins);
   
   for(Int_t i=0; i<nptbins; i++){
      TH1D *hM = 0x0, *hMR = 0x0;
      Int_t ptb[2] = {0, 0};
      if(massXAxis) {
      	 ptb[0] = hpTMass->GetYaxis()->FindBin(ptlims[i]);
      	 ptb[1] = hpTMass->GetYaxis()->FindBin(ptlims[i+1])-1;
      	 hM = hpTMass->ProjectionX(Form("hM_%d", i), ptb[0], ptb[1]);
      	 
      }
      else {
      	 ptb[0] = hpTMass->GetXaxis()->FindBin(ptlims[i]);
      	 ptb[1] = hpTMass->GetXaxis()->FindBin(ptlims[i+1])-1;
      	 hM = hpTMass->ProjectionY(Form("hM_%d", i), ptb[0], ptb[1]);
      }
      hM->SetLineWidth(2);
      
      if(massXAxisR) {
      	 ptb[0] = hpTMassR->GetYaxis()->FindBin(ptlims[i]);
      	 ptb[1] = hpTMassR->GetYaxis()->FindBin(ptlims[i+1])-1;
      	 hMR = hpTMassR->ProjectionX(Form("hMR_%d", i), ptb[0], ptb[1]);
      	 
      }
      else {
      	 ptb[0] = hpTMassR->GetXaxis()->FindBin(ptlims[i]);
      	 ptb[1] = hpTMassR->GetXaxis()->FindBin(ptlims[i+1])-1;
      	 hMR = hpTMassR->ProjectionY(Form("hMR_%d", i), ptb[0], ptb[1]);
      }
      hMR->SetLineWidth(2);
      hMR->SetMarkerStyle(20);
      hMR->SetLineColor(kGray);
      hMR->SetMarkerColor(kGray);
      
      hM->Scale(1./hM->Integral());
      hMR->Scale(1./hMR->Integral());
      TH1 *hMDividable = 0x0, *hMRDividable = 0x0;
      /* 
      Double_t binW = hM->GetBinWidth(5), binWR = hMR->GetBinWidth(5);
      Double_t binWRat = binW/binWR;
      
      if(binWRat<1) {
      	 binWRat = 1./binWRat;
      	 hM->Rebin(binWRat);
      } else hMR->Rebin(binWRat);
      
     
      if(binWRat == 1) {
      	 hMDividable = hM;
      	 hMRDividable = hMR;
      }
      else {
      	 Int_t nnewBin = TMath::Min(hM->GetNbinsX(), hMR->GetNbinsX());
      	 Int_t minnew = TMath::Max(hM->GetBinLowEdge(1), hMR->GetBinLowEdge(1));
      	 Int_t maxnew = TMath::Min(hM->GetBinLowEdge(hM->GetNbinsX()), hMR->GetBinLowEdge(hMR->GetNbinsX()));
      	 Printf("Modify histograms to have %d bins, range (%f, %f)", nnewBin, minnew, maxnew);
      	 hMDividable = new TH1F();
      }
      //information on binning
      if(i==0){
      	 Printf("Nbins: 1 -> %d, 2 -> %d, Range: 1 -> (%f, %f), 2 -> (%f, %f), BinWidth: 1 -> %f, 2 -> %f", hM->GetNbinsX(), hMR->GetNbinsX(), hM->GetBinLowEdge(1), hM->GetBinLowEdge(hM->GetNbinsX()+1), hMR->GetBinLowEdge(1), hMR->GetBinLowEdge(hMR->GetNbinsX()+1), hM->GetBinWidth(5), hMR->GetBinWidth(5));
      }
      */
      UniformTH1FForDivide((TH1F*)hM, (TH1F*)hMR, hMDividable, hMRDividable);
      
      hMDividable->SetLineWidth(2);
      hMDividable->SetLineColor(kBlack);
      hMRDividable->SetMarkerColor(kBlack);
      
      hMRDividable->SetMarkerStyle(20);
      hMRDividable->SetLineColor(kGray);
      hMRDividable->SetMarkerColor(kGray);
      
      
      TPaveText *pvt = new TPaveText(0.2, 0.75, 0.5, 0.9,"NDC");
      pvt->SetBorderSize(0);
      pvt->SetFillStyle(0);      
      pvt->AddText(Form("%.0f < %s < %.0f", ptlims[i], titleBins.Data() , ptlims[i+1]));
      	 
      cMvsPt->cd(i+1);
      hM->Draw();
      hMR->Draw("sames");
      pvt->Draw();
      
      cMvsPtBis->cd(i+1);
      hMDividable->Draw();
      hMRDividable->Draw("sames");
      pvt->Draw();
      
      
      TH1D *hMRatio = (TH1D*) hMDividable->Clone("hMRatio");
      hMRatio->Divide(hMRDividable);
      cMvsPtRatios->cd(i+1);
      hMRatio->Draw();
      pvt->Draw();
   }

   SaveCv(cMvsPtBis);
   SaveCv(cMvsPtRatios);
}
//_______________________________________________________________________________________________

void SuperimposeModifiedMassFixedInput(){
   Int_t nfiles = 2;
   TString paths[nfiles] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/Train855/analysis/", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Train856/analysis/"};
   TString filename = "ModifiedMassBgkFluct.root";
   TString leg[nfiles] = {"Single Track", "Single Track weighted by Ncoll"};
   TFile *fin[nfiles];

   //SuperimposeModifiedMass();
   Int_t nentries[nfiles];
   TH1F* hcompare[nfiles][20];
   for(Int_t i = 0; i<nfiles; i++){
      fin[i] = new TFile(Form("%s%s", paths[i].Data(), filename.Data()));
      if(!fin[i]->IsOpen()){
      	 Printf("Check paths or run plotOutputTaskJetShapeDeriv + PerformConvolution first");
      	 continue;
      }
      nentries[i] = fin[i]->GetNkeys();
      if(i>0 && nentries[i] != nentries[i-1]) {
      	 Printf("Problem!! different number of keys in file %d (%d instead of %d)", i, nentries[i],  nentries[i-1]);
      	 return;
      }
      TIter next(fin[i]->GetListOfKeys());
      TKey *key;
      Int_t count = 0;
      while(key = (TKey*)next()){
      	 TString hname = key->GetName();
      	 Printf ("%s", hname.Data());
      	 hcompare[i][count] = (TH1F*)fin[i]->Get(hname);
      	 count++;
      }
   }
   Bool_t notok=kFALSE;
   for(Int_t i = 0; i<nentries[0]; i++){
      
      for(Int_t j = 0; j< nfiles; j++){
      	 
      	 if(!hcompare[i][j]) {
      	    Printf("Histogram not found");
      	    notok = kTRUE;
      	    continue;
      	 }
      }
      TCanvas *c=new TCanvas(Form("c%d", i), Form("Plot bin %d", i));
      if(!notok) CompareDistributions(c->cd(), hcompare[0][i], hcompare[1][i]); //can handle only two!!
   
   }
   
}


//--------------------------------------------------------------------------------------------
void plotOutputTaskJetShapeDeriv2(TString filename = "AnalysisResults.root", TString listname = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", Int_t nsuperimp = 0, TString suplists[]=0){
   
   TCanvas * cmass= new TCanvas("cmass", "Mass subtracted", 800, 800);
   for(Int_t s = -1; s<nsuperimp; s++){
      
      TList * list = 0x0;
      if(s<0) {
      	 list = ReadFile(filename, listname);
      	 
      	 
      	 TString h2name = "fh2MSubPtRawAll_0";
      	 TH2F *h2MSubPtRaw = (TH2F*)list->FindObject(h2name);
      	 TCanvas *cEventLook = new TCanvas("cEventLook", "pPb + embedded tracks", 800, 800);
      	 cEventLook->cd();
      	 gPad->SetLogz();
      	 h2MSubPtRaw->Draw("colz");
      	 SaveCv(cEventLook, "", saveLv);
      } else list = ReadFile(filename, suplists[s]);
      
      TString h3name = "fh3MSubPtRawDRMatch_0";
      
      TH3F *h3MSubPtRawDRMatch = (TH3F*)list->FindObject(h3name);
      
      TH1F *hDR = (TH1F*)h3MSubPtRawDRMatch->Project3D("z");
      TH1F *hMassSub = (TH1F*)h3MSubPtRawDRMatch->Project3D("x");
      Printf("Bin mass width %.2f", h3MSubPtRawDRMatch->GetXaxis()->GetBinWidth(3));
      hMassSub->SetName("hMassSub");
      hMassSub->SetLineWidth(2);
      hMassSub->SetLineStyle(s+2);
      TCanvas *cDR = new TCanvas("cDR", "Distance embedded track jet", 800, 800);
      cDR->cd();
      hDR->SetLineStyle(s+2);
      hDR->Draw("sames");
      
      
      cmass->cd();
      hMassSub->Scale(1./hMassSub->Integral());
      hMassSub->Draw("sames");
      TPaveText *pvmean = new TPaveText(0.4, 0.5, 0.8, 0.6, "NDC");
      pvmean->SetFillStyle(0);
      pvmean->SetBorderSize(0);
      pvmean->SetTextColor(hMassSub->GetLineColor());
      pvmean->AddText(Form("Mean = %.3f #pm %.3f", hMassSub->GetMean(), hMassSub->GetMeanError()));
      cmass->cd();
      pvmean->Draw();
      
      Double_t pTlims[4] = {40., 60., 80., 100.};
      for(Int_t ipt = 0; ipt < 3; ipt++){
      	 Int_t biny1 = h3MSubPtRawDRMatch->GetYaxis()->FindBin(pTlims[ipt]);
      	 Int_t biny2 = h3MSubPtRawDRMatch->GetYaxis()->FindBin(pTlims[ipt+1]) -1;
      	 TH1F *hMassSubpTSel = (TH1F*)h3MSubPtRawDRMatch->ProjectionX(Form("hMassSubpT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]), biny1, biny2, 0, -1);
      	 
      	 
      	 hMassSubpTSel->SetLineColor(colors[ipt]);
      	 hMassSubpTSel->SetLineStyle(s+2);
      	 cmass->cd();
      	 hMassSubpTSel->Scale(1./hMassSubpTSel->Integral());
      	 hMassSubpTSel->Draw("sames");
      	 
      	 TPaveText *pvmeanpT = new TPaveText(0.4, 0.1+ (0.1*ipt), 0.8, 0.1+(0.1*(ipt+1)), "NDC");
      	 pvmeanpT->SetFillStyle(0);
      	 pvmeanpT->SetBorderSize(0);
      	 pvmeanpT->SetTextColor(hMassSubpTSel->GetLineColor());
      	 pvmeanpT->AddText(Form("Mean (%.0f-%.0f GeV/#it{c}) = %.3f #pm %.3f",pTlims[ipt], pTlims[ipt+1], hMassSubpTSel->GetMean(), hMassSubpTSel->GetMeanError()));
      	 cmass->cd();
      	 pvmeanpT->Draw();
      	 
      }
   }
   SaveCv(cmass, "", saveLv);
   
}

//-------------------------------------------------------------------------------------------

void MatrixDeltaMDeltapT(TString filename = "AnalysisResults.root", TString listname = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", TString suffix = ""){
   //
   // Draw several projections of THnSparse
   // dM, dpT, pTpart/det, rho, rho_m...
   //
   
   
   //Int_t nrhobins = 9;
   //Double_t rholims[nrhobins+1] = {0, 0.4, 0.8, 1.6, 2, 3.2, 4.8, 8, 11.2, 14.8}; //some troubles with the last bin
   Int_t nrhobins = 8;
   Double_t rholims[nrhobins+1] = {0, 0.4, 0.8, 1.6, 2, 3.2, 4.8, 8, 14.8}; 
   Int_t nrhombins = 7;
   Double_t rhomlims[nrhombins+1] = {0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.3, 0.5}; 

   TList * list = ReadFile(filename, listname);
   TString spname = "fhnDeltaMass_0";
   if(!list){
      Printf("Exit");
      return;
   }
   THnSparse *hsp = (THnSparse*)list->FindObject(spname);
   if(!hsp){
      Printf("%s  not found", spname.Data());
      return;
   }
   
   //save unfolding matrix (this must be taken from PYTHIA not track embedding) (tbc)
   TFile *funf = new TFile(Form("UnfoldingThnSparse%s.root", suffix.Data()), "recreate");
   Int_t axunf[4] = {2, 3, 4, 5}; //Mrec, Mpart, pTrec, pTpart 
   THnSparse *h4Dunf = (THnSparse*)hsp->ProjectionND(4, axunf);
   funf->cd();
   h4Dunf->Write();
   
   //integrated 2D projection dM vs dpT
   
   TH2D *hdMdpT =(TH2D*) hsp->Projection(0,1);
   Int_t axfluc[4] = {0, 1, 3, 5}; // dM, dpT, pTpart, Mpart
   THnSparse *h4DFluc = (THnSparse*)hsp->ProjectionND(4, axfluc);
   //save fluctuations (this comes from track embedding on data) (tbc)
   funf->cd();
   hdMdpT->Write();
   h4DFluc->Write();
   TCanvas *c = new TCanvas("c", "dM dpT", 800,800);
   
   hdMdpT->Draw("colz");
   
   SaveCv(c);
   Int_t naxes = hsp->GetNdimensions();
   Int_t nx, ny, dx, dy;
   CalculatePads(naxes, nx, ny, dx, dy);
   TCanvas *cpj = new TCanvas("cpj", "Projections", dx, dy);
   cpj->Divide(nx, ny);
   TString title[naxes];
   for(Int_t i=0; i<naxes; i++){
      title[i] = hsp->GetAxis(i)->GetTitle();
      
      if (i==naxes-1){
      	 title[i] = "#it{N}_{const}^{part}";//"#it{p}_{T,part}";
      	 hsp->GetAxis(i)->SetTitle(title[i]);
      }
      Printf("Axis %d = %s", i, title[i].Data());
      //1D projections
      cpj->cd(i+1);
      TH1D* hp1d = hsp->Projection(i);
      if(i == 0 || i == 2) hp1d->GetXaxis()->SetRangeUser(-10, 10);
      hp1d->Scale(1./hp1d->Integral("width"));
      hp1d->Draw();
      if(i == 2) {
      	 cpj->cd(1);
      	 hp1d->Draw("sames");
      }
      //Printf("End");
   }
   //Printf("Out");
   //SaveCv(cpj);
   //Printf("Saved");
   //return;
   Int_t ax2Dproj[2] = {0, 1}; //y, x (dM, dpT)
   // 2D projection dM vs dpT in bins of pTrec
   Int_t npTbins = 5;
   Double_t pTlims[npTbins+1] = {10, 20, 40, 60, 80, 100};
   Int_t nMbins = 3;
   Double_t pMlims[nMbins+1] = {0, 5, 10, 15};
   Int_t ndMbins = 3;
   //Double_t pdMlims[ndMbins+1] = {-10, -5, 5, 10};
   Double_t pdMlims[ndMbins+1] = {0, 5, 10, 15};
   Int_t ndNbins = 4;
   Double_t dNlims[ndNbins+1] = {-5, 0, 1, 2, 5};
   Int_t nNbins = 3;
   Double_t Nlims[ndNbins+1] = {0, 1, 2, 5};
   
   Int_t axSel = 5; //pTTrue 4; //pTrec
   TString title2D = Form("%sVs%s", title[ax2Dproj[0]].Data(), title[ax2Dproj[1]].Data());
   Int_t ax1Dproj;
   
   Draw2DProjectionsinBins(hsp, ax2Dproj, title2D, axSel, title[axSel], npTbins, pTlims);
   Printf("0");   
   ax2Dproj[0] = 2; //M_det
   ax2Dproj[1] = 3; //M_part
   title2D = Form("%sVs%s", title[ax2Dproj[0]].Data(), title[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hsp, ax2Dproj, title2D, axSel, title[axSel], npTbins, pTlims);
   Printf("1");
   ax2Dproj[0] = 4; //pT_det
   ax2Dproj[1] = 5; //pT_part
   axSel = 3; //M_part
   title2D = Form("%sVs%s", title[ax2Dproj[0]].Data(), title[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hsp, ax2Dproj, title2D, axSel, title[axSel], nMbins, pMlims);  
   Printf("2");
   ax2Dproj[0] = 2; //M_det
   ax2Dproj[1] = 3; //M_part
   axSel = 7; //dN
   title2D = Form("%sVs%s", title[ax2Dproj[0]].Data(), title[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hsp, ax2Dproj, title2D, axSel, title[axSel], ndNbins, dNlims);
   Printf("3");
   ax2Dproj[0] = 2; //M_det
   ax2Dproj[1] = 4; //pT_det
   axSel = -1; //no bins
   Draw2DProjectionsinBins(hsp, ax2Dproj, title2D, axSel, title[axSel], npTbins, pTlims);
   ax2Dproj[0] = 3; //M_part
   ax2Dproj[1] = 5; //pT_part
   ax1Dproj = 0; //dM
   axSel = -1; //no bins
   Draw2DProjectionsinBins(hsp, ax2Dproj, title2D, axSel, title[axSel], npTbins, pTlims);
   Printf("4");
   ax1Dproj = 1; //dpT
   axSel = 5; //pTTrue
   Draw1DProjectionsinBins(hsp, ax1Dproj, title[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kFALSE);
   Printf("5");
   ax1Dproj = 4; //pTrec
   axSel = 5; //pTTrue
   Draw1DProjectionsinBins(hsp, ax1Dproj, title[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kFALSE);
   Printf("6");
   ax1Dproj = 7; //DeltaNconst
   axSel = 5; //pTTrue
   Draw1DProjectionsinBins(hsp, ax1Dproj, title[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kFALSE);
   Printf("7");
   ax1Dproj = 0; // dM
   axSel = 7; //DeltaNconst
   Draw1DProjectionsinBins(hsp, ax1Dproj, title[ax1Dproj], axSel, title[axSel], ndNbins, dNlims, kFALSE);
   Printf("8");
   ax1Dproj = 7; //DeltaNconst
   axSel = 8; //Nconst
   Draw1DProjectionsinBins(hsp, ax1Dproj, title[ax1Dproj], axSel, title[axSel], nNbins, Nlims, kFALSE);
   
   Printf("Done hsp");
   TString spnameb = "fhnDeltaMassAndBkgInfo";
   THnSparse *hspbkg = (THnSparse*)list->FindObject(spnameb);
   if(!hspbkg){
      Printf("%s  not found", spnameb.Data());
      return;
   }
   
   Int_t naxesbkg = hspbkg->GetNdimensions();
   TString titlebkg[naxesbkg];
   for(Int_t i=0; i<naxesbkg; i++){
      titlebkg[i] = hspbkg->GetAxis(i)->GetTitle();
      
      if (i==naxesbkg-1){
      	 titlebkg[i] = "#rho_{m}";
      	 hspbkg->GetAxis(i)->SetTitle(titlebkg[i]);
      }
      Printf("Axis %d = %s", i, titlebkg[i].Data());
   }
   
   axSel = 5; //pTpart
   ax2Dproj[0] = 0; ax2Dproj[1] = 6 ; //rho -> dM vs rho
   title2D = Form("%sVs%s", titlebkg[ax2Dproj[0]].Data(), titlebkg[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hspbkg, ax2Dproj, title2D, axSel, titlebkg[axSel], npTbins, pTlims);
      
   ax2Dproj[0] = 0; ax2Dproj[1] = 7 ; //rho_m  -> dM vs rho_m
   title2D = Form("%sVs%s", titlebkg[ax2Dproj[0]].Data(), titlebkg[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hspbkg, ax2Dproj, title2D, axSel, titlebkg[axSel], npTbins, pTlims);

   axSel = 5; //pTpart
   ax2Dproj[0] = 1; ax2Dproj[1] = 6 ; //rho -> dM vs rho
   title2D = Form("%sVs%s", titlebkg[ax2Dproj[0]].Data(), titlebkg[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hspbkg, ax2Dproj, title2D, axSel, titlebkg[axSel], npTbins, pTlims);
      
   ax2Dproj[0] = 1; ax2Dproj[1] = 7 ; //rho_m  -> dM vs rho_m
   title2D = Form("%sVs%s", titlebkg[ax2Dproj[0]].Data(), titlebkg[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hspbkg, ax2Dproj, title2D, axSel, titlebkg[axSel], npTbins, pTlims);
      
   axSel = 0;
   ax2Dproj[0] = 6; ax2Dproj[1] = 7 ; //rho vs rho_m
   title2D = Form("%sVs%s", titlebkg[ax2Dproj[0]].Data(), titlebkg[ax2Dproj[1]].Data());
   Draw2DProjectionsinBins(hspbkg, ax2Dproj, title2D, axSel, titlebkg[axSel], ndMbins, pdMlims);
      
   ax1Dproj = 0; //dM
   axSel = 6; //rho
   Draw1DProjectionsinBins(hspbkg, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhobins, rholims, kFALSE);
      
   Printf("As a funciton of rhom,");
   ax1Dproj = 0; //dM
   axSel = 7; //rho_m
   Draw1DProjectionsinBins(hspbkg, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kFALSE);
      
   ax1Dproj = 1; //dpT
   axSel = 6; //rho
   Draw1DProjectionsinBins(hspbkg, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhobins, rholims, kFALSE);
       
   ax1Dproj = 1; //dpT
   axSel = 7; //rho_m
   Draw1DProjectionsinBins(hspbkg, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kFALSE);

   ax1Dproj = 3; //m unsub
   axSel = 6; //rho
   Draw1DProjectionsinBins(hspbkg, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhobins, rholims, kFALSE);
   ResetAxisSel(hspbkg);

   axSel = 7; //rho_m
   Draw1DProjectionsinBins(hspbkg, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kFALSE);
   
   ax1Dproj = 5; //pTunsub
   Draw1DProjectionsinBins(hspbkg, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kTRUE, kFALSE);
   
   
}

void CompareMatrixDeltaMDeltapT(TString filename1 = "AnalysisResults.root", TString filename2 = "AnalysisResults.root", TString listname1 = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", TString listname2 = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCnoOvl", TString leg1 = "All", TString leg2 = "OverlapExclusion", Int_t norm = 1 /*0 = don't scale the histograms , 1 = scale by the number of events, 2 = scale to integral 1*/){
   //
   // Draw several projections of THnSparse
   // dM, dpT, pTpart/det, rho, rho_m...
   //
   
   TPaveText *pvNorm = new TPaveText(0.11, 0.12, 0.31, 0.22, "NDC");
   pvNorm->SetFillStyle(0);
   pvNorm->SetBorderSize(0);
   
   TString normS = "#it{N}_{ev}";
   if(norm == 0) normS = "None";
   if(norm == 2) normS = "1";
   pvNorm->AddText(Form("Norm = %s", normS.Data()));
   Printf("Info::::::::Normalisation %s", normS.Data());
   //Int_t nrhobins = 9;
   //Double_t rholims[nrhobins+1] = {0, 0.4, 0.8, 1.6, 2, 3.2, 4.8, 8, 11.2, 14.8}; //some troubles with the last bin
   Int_t nrhobins = 8;
   Double_t rholims[nrhobins+1] = {0, 0.4, 0.8, 1.6, 2, 3.2, 4.8, 8, 14.8}; 
   Int_t nrhombins = 7;
   Double_t rhomlims[nrhombins+1] = {0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.3, 0.5}; 

   TList * list1 = ReadFile(filename1, listname1);
   if(!list1){
      Printf("Exit");
      return;
   }
   TList* list2 = ReadFile(filename2, listname2);
   if(!list2){
      Printf("Exit");
      return;
   }
   Double_t nevents1 = 1, nevents2 = 1;
   TString hneventsname = "fHistEventCount";
   TH1F *hnevents = (TH1F*) list1->FindObject(hneventsname);
   if(!hnevents) {
      Printf("Warning! Cannot properly normalize input 1");
      
   }
   nevents1 = hnevents->GetBinContent(1);
   
   hnevents = (TH1F*) list2->FindObject(hneventsname);
   if(!hnevents) {
      Printf("Warning! Cannot properly normalize input 1");
      
   }
   nevents2 = hnevents->GetBinContent(1);
   if(norm == 0) {
      nevents1 = 1;
      nevents2 = 1;
   }
   if(norm == 2) {
      nevents1 = -1;
      nevents2 = -1;
   }
   TString spname = "fhnDeltaMass_0";
   THnSparse *hsp1 = (THnSparse*)list1->FindObject(spname);
   if(!hsp1){
      Printf("%s  not found", spname.Data());
      return;
   }
   Int_t naxes = hsp1->GetNdimensions();
   TString title[naxes];
   for(Int_t i=0; i<naxes; i++){
      title[i] = hsp1->GetAxis(i)->GetTitle();
      
      if (i==naxes-1){
      	 title[i] = "#it{p}_{T,part}";;
      	 hsp1->GetAxis(i)->SetTitle(title[i]);
      }
      Printf("Axis %d = %s", i, title[i].Data());
   }
   

   THnSparse *hsp2 = (THnSparse*)list2->FindObject(spname);
   if(!hsp2){
      Printf("%s  not found", spname.Data());
      return;
   }

   //Int_t npTbins = 2;
   //Double_t pTlims[npTbins+1] = {5, 10., 100}; // 20, 40, 60, 80,
   Int_t npTbins = 3;
   Double_t pTlims[npTbins+1] = {40, 60, 80, 100}; // 20, 40, 60, 80,
   Int_t nMbins = 3;
   Double_t pMlims[nMbins+1] = {0, 5, 10, 15};
   Int_t ndMbins = 3;
   //Double_t pdMlims[ndMbins+1] = {-10, -5, 5, 10};
   Double_t pdMlims[ndMbins+1] = {0, 5, 10, 15};
   Int_t ndNbins = 4;
   Double_t dNlims[ndNbins+1] = {-5, 0, 1, 2, 5};
   Int_t nNbins = 3;
   Double_t Nlims[ndNbins+1] = {0, 1, 2, 5};
   Int_t axSel = 5; //pTTrue 4; //pTrec

   Int_t ax1Dproj = 0; //dM

   Compare1DProjectionsinBins(hsp1, hsp2, ax1Dproj, title[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kTRUE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2, kTRUE);
   gPad->cd();
   pvNorm->Draw();   
   ax1Dproj = 1; //dpT
   axSel = 5; //pTTrue
   Compare1DProjectionsinBins(hsp1, hsp2, ax1Dproj, title[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kTRUE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2, kTRUE);
      
   ax1Dproj = 4; //pTrec
   axSel = 5; //pTTrue
   Compare1DProjectionsinBins(hsp1, hsp2, ax1Dproj, title[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);

      
   ax1Dproj = 7; //DeltaNconst
   axSel = 5; //pTTrue
   Compare1DProjectionsinBins(hsp1, hsp2, ax1Dproj, title[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   
   ax1Dproj = 0; // dM
   axSel = 7; //DeltaNconst
   Compare1DProjectionsinBins(hsp1, hsp2, ax1Dproj, title[ax1Dproj], axSel, title[axSel], ndNbins, dNlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   
   ax1Dproj = 7; //DeltaNconst
   axSel = 8; //Nconst
   Compare1DProjectionsinBins(hsp1, hsp2, ax1Dproj, title[ax1Dproj], axSel, title[axSel], nNbins, Nlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   
   TString spnameb = "fhnDeltaMassAndBkgInfo";
  
   THnSparse *hspbkg1 = (THnSparse*)list1->FindObject(spnameb);
   if(!hspbkg1){
      Printf("%s  not found", spnameb.Data());
      return;
   }
   Int_t naxesbkg = hspbkg1->GetNdimensions();
   TString titlebkg[naxesbkg];
   for(Int_t i=0; i<naxesbkg; i++){
      titlebkg[i] = hspbkg1->GetAxis(i)->GetTitle();
      
      if (i==naxesbkg-1){
      	 titlebkg[i] = "#rho_{m}";
      	 hspbkg1->GetAxis(i)->SetTitle(titlebkg[i]);
      }
      Printf("Axis %d = %s", i, titlebkg[i].Data());
   }
   

   THnSparse *hspbkg2 = (THnSparse*)list2->FindObject(spnameb);
   if(!hspbkg2){
      Printf("%s  not found", spnameb.Data());
      return;
   }
   
   ax1Dproj = 0; //dM
   axSel = 5; //pTTrue
   axSel = 6; //rho
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhobins, rholims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);

   ax1Dproj = 0; //dM
   axSel = 6; //rho
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhobins, rholims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
      
   Printf("As a funciton of rhom,");
   ax1Dproj = 0; //dM
   axSel = 7; //rho_m
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
      
   ax1Dproj = 1; //dpT
   axSel = 6; //rho
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhobins, rholims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
       
   ax1Dproj = 1; //dpT
   axSel = 7; //rho_m
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);

   ax1Dproj = 3; //m unsub
   axSel = 6; //rho
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhobins, rholims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   ResetAxisSel(hspbkg1); ResetAxisSel(hspbkg2);

   axSel = 7; //rho_m
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   
   ax1Dproj = 5; //pTunsub
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], nrhombins, rhomlims, kTRUE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   
   axSel = 0;  //dM
   ax1Dproj = 6; //rho
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], ndMbins, pdMlims, kTRUE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   ax1Dproj = 7; //rho_m
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, titlebkg[axSel], ndMbins, pdMlims, kTRUE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);
   
   ax1Dproj = 3; //Mrec
   axSel = 5; //pTTrue
   Compare1DProjectionsinBins(hspbkg1, hspbkg2, ax1Dproj, titlebkg[ax1Dproj], axSel, title[axSel], npTbins, pTlims, kFALSE, kFALSE, 1./nevents1, 1./nevents2, leg1, leg2);

}
//--------------------------------------------------------------------------------------------

void DrawdMvsMassSingle(Bool_t less = kFALSE){
   Int_t nfiles = 6;
   TString base = "/data/Work/jets/JetMass/pPbJetMassAnalysis";
   TString trains[6] = {"Train921/runlist1", "Train966", "Train967", "Train968", "Train969", "Train970"};
   Double_t massSingle[6] = {0, 1, 2, 3, 5, 8};
   if(less){
      trains[0] = "Train921/runlist1";
      trains[1] = "Train967";
      trains[2] = "Train970";
      trains[3] = trains[4] = trains[5] = "";
      massSingle[0] = 0;
      massSingle[1] = 2;
      massSingle[2] = 8;
      massSingle[3] = massSingle[4] = massSingle[5] = 0;
      nfiles = 3;
   }
   
   TString pathfiles[nfiles], legT[nfiles];
   TString listname = "JetShapeConst_JetEmb_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC", 
   hspname = "fhnDeltaMass_0";
   
   
   Printf("Set up");
   Bool_t reverse = kTRUE;
   if(less) reverse = kFALSE;
   for(Int_t i = 0; i<nfiles; i++){
      Int_t j = i;
      if(reverse) {
      	 j = nfiles-1-i;
      }
      pathfiles[i] = Form("%s/%s/AnalysisResults.root", base.Data(), trains[j].Data());
      legT[i] = Form("Mass = %.0f GeV", massSingle[j]);
      Printf("%d %s, %s", i, pathfiles[i].Data(), legT[i].Data());
   }

   Int_t npTbins = 2;
   Double_t pTlims[npTbins+1] = {40, 80, 120}; // 20, 40, 60, 80,

   Int_t nMbins = 3;
   Double_t Mlims[nMbins+1] = {0, 5, 10, 15};
   
   Int_t ax1Dproj = 0, 
   axSel = 5; // dM, pTTrue
   DrawCompare1DProjectionsinBins(nfiles, pathfiles, legT, listname, hspname, ax1Dproj, axSel, npTbins, pTlims, 2, 1);
   
   ax1Dproj = 0; // dM
   axSel = 2; // M_det
   DrawCompare1DProjectionsinBins(nfiles, pathfiles, legT, listname, hspname, ax1Dproj, axSel, nMbins, Mlims, 2, 1, kTRUE);
   
   TString hspnameb = "fhnDeltaMassAndBkgInfo";
   Int_t nrhobins = 4;
   Double_t rholims[nrhobins+1] = {0,  2, 4.8, 8, 14.8}; //0.4, 0.8,1.6, 3.2, 
   
   ax1Dproj = 0; // dM
   axSel = 6; // rho
   DrawCompare1DProjectionsinBins(nfiles, pathfiles, legT, listname, hspnameb, ax1Dproj, axSel, nrhobins, rholims, 2, 1, kFALSE);
}

void DrawdMvsMassSingleToy(){
   Int_t nfiles = 5;
   TString base = "/data/Work/jets/JetMass/BkgFluctStudies/CleanEmbedding";
   TString trains[5] = {"SingleTrackInThermal20k", "SingleMassive3GeVTrackInThermal20k", "SingleMassive4GeVTrackInThermal20k", "SingleMassive5GeVTrackInThermal20k", "SingleMassive8GeVTrackInThermal20k"};
   Double_t massSingle[5] = {0, 3, 4, 5, 8};
   TString pathfiles[nfiles], legT[nfiles];
   TString listname = "JetShapeConst_Jet_AKTChargedR040_ThrmTracksSingle_pT0150_E_scheme_TCMCMatch", 
   hspname = "fhnDeltaMass_0";
   Printf("Set up");
   Bool_t reverse = kTRUE;
   for(Int_t i = 0; i<nfiles; i++){
      Int_t j = i;
      if(reverse) {
      	 j = nfiles-1-i;
      }
      pathfiles[i] = Form("%s/%s/merge/AnalysisResults.root", base.Data(), trains[j].Data());
      legT[i] = Form("Mass = %.0f GeV", massSingle[j]);
      Printf("%d %s, %s", i, pathfiles[i].Data(), legT[i].Data());
   }

   Int_t npTbins = 2;
   Double_t pTlims[npTbins+1] = {40, 80, 120}; // 20, 40, 60, 80,

   Int_t ax1Dproj = 0, 
   axSel = 5; // dM, pTTrue
   DrawCompare1DProjectionsinBins(nfiles, pathfiles, legT, listname, hspname, ax1Dproj, axSel, npTbins, pTlims, 2, 1);
}
//--------------------------------------------------------------------------------------------

void DrawCompare1DProjectionsinBins(Int_t nfiles, TString pathfiles[], TString legT[], TString listname, TString hspname, Int_t ax1Dproj, Int_t axSel, Int_t npTbins, Double_t pTlims[], Int_t norm /*0 = don't scale the histograms , 1 = scale by the number of events, 2 = scale to integral 1*/, Int_t rebinFactor , Bool_t logy, Bool_t makeRatio){
   
   TPaveText *pvNorm = new TPaveText(0.11, 0.12, 0.31, 0.22, "NDC");
   pvNorm->SetFillStyle(0);
   pvNorm->SetBorderSize(0);
   
   TString normS = "#it{N}_{ev}";
   if(norm == 0) normS = "None";
   if(norm == 2) normS = "1";
   pvNorm->AddText(Form("Norm = %s", normS.Data()));
   Printf("Info::::::::Normalisation %s", normS.Data());


   //read the first file to get general information
   
   TList * list = ReadFile(pathfiles[0], listname);
   if(!list){
      Printf("Exit");
      return;
   }
   THnSparse *hsp = (THnSparse*)list->FindObject(hspname);
   if(!hsp){
      Printf("%s  not found", hspname.Data());
      return;
   }
   
   
   ResetAxisSel(hsp);
   
   Int_t ndim = hsp->GetNdimensions();
   if(ax1Dproj >= ndim || axSel>= ndim) {
      Printf("Requested axis %d, %d out of range %d", ax1Dproj, axSel, ndim);
      return; 
   }
   TString title[ndim];
   for(Int_t i=0; i<ndim; i++){
      title[i] = hsp->GetAxis(i)->GetTitle();
      
      if (i==ndim-1){
      	 title[i] = "#it{p}_{T,part}";;
      	 hsp->GetAxis(i)->SetTitle(title[i]);
      }
      Printf("Axis %d = %s", i, title[i].Data());
   }
   
   
   Int_t nx, ny, dx, dy;
   CalculatePads(npTbins, nx, ny, dx, dy);
   TCanvas *c1DpjpTrecsel = new TCanvas(Form("c1Dpj%d_Bin%d", ax1Dproj, axSel), Form("Projection 1D %s in %s bins", title[ax1Dproj].Data(), title[axSel].Data()), dx, dy);
   if(npTbins>1) c1DpjpTrecsel->Divide(nx, ny);
   
   TCanvas *c1DpjpTrecselRatio=0x0;
   if(makeRatio) {
      c1DpjpTrecselRatio = new TCanvas(Form("c1Dpj%d_Bin%dRatio",  ax1Dproj, axSel), Form("Ratio Projection 1D %s in %s bins", title[ax1Dproj].Data(), title[axSel].Data()), dx, dy);
      if(npTbins>1) c1DpjpTrecselRatio->Divide(nx, ny);
   }
   
   TLegend *leg = new TLegend(0.58, 0.5, 0.9, 0.8);
      
   TH1F *h1First[npTbins];
   Bool_t drawfirst[npTbins];
   for(Int_t i = 0; i< npTbins; i++) {
      h1First[i] = 0;
      drawfirst[i]= kTRUE;
   }
   
   Double_t maxR = 0;
   
   for(Int_t ifile = 0; ifile < nfiles; ifile++){ // loop over the input files
      TList * list = ReadFile(pathfiles[ifile], listname);
      if(!list){
      	 Printf("Exit");
      	 return;
      }
      THnSparse *hsp = (THnSparse*)list->FindObject(hspname);
      if(!hsp){
      	 Printf("%s  not found", hspname.Data());
      	 return;
      }
      
      ResetAxisSel(hsp);
      

      for(Int_t i=0; i<npTbins; i++){
      	 ResetAxisSel(hsp);
      	 Int_t binsel[2] = {hsp->GetAxis(axSel)->FindBin(pTlims[i]), hsp->GetAxis(axSel)->FindBin(pTlims[i+1])-1};
      	 //if (debug) Printf("SetRange axis %d %d (%f) - %d (%f)", axSel, binsel[0], pTlims[i], binsel[1], pTlims[i+1]);
      	 TPaveText *pvt = new TPaveText(0.2, 0.75, 0.5, 0.9,"NDC");
      	 pvt->SetBorderSize(0);
      	 pvt->SetFillStyle(0);
      	 if(pTlims[i+1]>10)
      	    pvt->AddText(Form("%.0f < %s < %.0f", pTlims[i],title[axSel].Data(), pTlims[i+1]));
      	 else if(pTlims[i+1]>1)                             
      	    pvt->AddText(Form("%.1f < %s < %.1f", pTlims[i],title[axSel].Data(), pTlims[i+1]));
      	 else                                               
      	    pvt->AddText(Form("%.2f < %s < %.2f", pTlims[i],title[axSel].Data(), pTlims[i+1]));
      	 
      	 hsp->GetAxis(axSel)->SetRange(binsel[0], binsel[1]);
      	 TH1D *h1 = hsp->Projection(ax1Dproj);
      	 h1->SetName(Form("%s_%d_Sel%dbin%d_f%d", hsp->GetName(), ax1Dproj, axSel, i, ifile));
      	 h1->SetLineColor(colors[ifile]);
      	 h1->SetFillColor(colors[ifile]);
      	 Int_t silsty;
      	 if(!(ifile % 2)) silsty = 3004;
      	 else silsty = 3005;
      	 h1->SetFillStyle(silsty);
      	 h1->SetLineWidth(2);
      	       	 
      	 h1->Rebin(rebinFactor);
      	 Double_t integral = h1->Integral();
      	 if(norm > 0) h1->Scale(norm);
      	 else if (integral > 0) h1->Scale(1./integral);
      	 
      	 if(ax1Dproj == 0) h1->GetXaxis()->SetRangeUser(-5., 15.); //dM

      	 c1DpjpTrecsel->cd(i+1);
      	 if(logy) gPad->SetLogy();
      	 if(ifile == 0) h1->Draw();
      	 else h1->Draw("sames");
      	 pvt->Draw();
      	 if(!h1First[i]) {
      	    h1First[i] = (TH1F*)h1->Clone(Form("%s_%d_Sel%dbin%d_Den", hsp->GetName(), ax1Dproj, axSel, i));
      	 }
      	 if(makeRatio) {
      	    TH1D *h1Ratio = (TH1D*)h1->Clone(Form("%sRatio", h1->GetName()));
      	    h1Ratio->Divide(h1First[i]);
      	    h1Ratio->GetYaxis()->SetRangeUser(0.,5.);
      	    h1Ratio->SetLineColor(colors[ifile]);
      	    h1Ratio->SetFillColor(0);
      	    h1Ratio->GetYaxis()->SetRangeUser(0., 70.);
      	    if(ax1Dproj == 0) h1Ratio->GetXaxis()->SetRangeUser(-5., 20.); //dM
      	    c1DpjpTrecselRatio->cd(i+1);
      	    if(drawfirst[i]) {
      	       h1Ratio->Draw("hist");
      	       
      	    }
      	    else h1Ratio->Draw("histsames");
      	    pvt->Draw();
      	    drawfirst[i] = kFALSE;
      	 }
      	 
      	 if (i==0) {
      	    leg->AddEntry(h1, legT[ifile], "L");
      	    c1DpjpTrecsel->cd(i+1);
      	    pvNorm->Draw();
      	    c1DpjpTrecselRatio->cd(i+1);
      	    pvNorm->Draw();
      	 }
      	 if (i == npTbins - 1) {
      	    c1DpjpTrecsel->cd(i+1);
      	    leg->Draw();
      	    if(makeRatio){
      	       c1DpjpTrecselRatio->cd(i+1);
      	       leg->Draw();
      	    }
      	 }
      } // bins
      
      
   } //loop on files
   
   SaveCv(c1DpjpTrecsel);
   if(makeRatio) SaveCv(c1DpjpTrecselRatio);
}
//--------------------------------------------------------------------------------------------

void Draw2DProjectionsinBins(THnSparse *hsp, Int_t ax2Dproj[], TString title, Int_t axSel, TString axSeltitle, Int_t npTbins, Double_t pTlims[]){
   ResetAxisSel(hsp);
   Int_t ndim = hsp->GetNdimensions();
   if(ax2Dproj[0] >= ndim || ax2Dproj[1] >= ndim || axSel>= ndim) {
      Printf("Requested axis %d, %d, %d,out of range %d", ax2Dproj[0], ax2Dproj[1], axSel, ndim);
      return; 
   }

   Int_t nx, ny, dx, dy;
   CalculatePads(npTbins, nx, ny, dx, dy);
   TCanvas *c2DpjpTrecsel = new TCanvas(Form("c2D%spj%d%d_bin%d", hsp->GetName(), ax2Dproj[0], ax2Dproj[1], axSel), Form("Projection 2D %s in %s bins", title.Data(), axSeltitle.Data()), dx, dy);
   if(axSel>=0) c2DpjpTrecsel->Divide(nx, ny);
   
   for(Int_t i=0; i<npTbins; i++){
      if(axSel >= 0){
      	 Int_t binsel[2] = {hsp->GetAxis(axSel)->FindBin(pTlims[i]), hsp->GetAxis(axSel)->FindBin(pTlims[i+1])-1};
      	 hsp->GetAxis(axSel)->SetRange(binsel[0], binsel[1]);
      }
      TPaveText *pvt = new TPaveText(0.2, 0.75, 0.5, 0.9,"NDC");
      pvt->SetBorderSize(0);
      pvt->SetFillStyle(0);
      if(pTlims[i+1]>10)
      	 pvt->AddText(Form("%.0f < %s < %.0f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      else if(pTlims[i+1]>1)
      	 pvt->AddText(Form("%.1f < %s < %.1f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      else 
      	 pvt->AddText(Form("%.2f < %s < %.2f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      
      c2DpjpTrecsel->cd(i+1);
      
      TH2D *h2 = hsp->Projection(ax2Dproj[0], ax2Dproj[1]);
      h2->SetName(Form("%s_%d_%d_Sel%dbin%d", hsp->GetName(), ax2Dproj[0], ax2Dproj[1], axSel, i));
      h2->Draw("colz");
      pvt->Draw();
      if(axSel<0) break;
   }
   
   SaveCv(c2DpjpTrecsel);
   
}
//--------------------------------------------------------------------------------------------
void Draw1DProjectionsinBins(THnSparse *hsp, Int_t ax1Dproj, TString title, Int_t axSel, TString axSeltitle, Int_t npTbins, Double_t pTlims[], Bool_t logy, Bool_t debug){
   ResetAxisSel(hsp);
   Int_t ndim = hsp->GetNdimensions();
   if(ax1Dproj >= ndim || axSel>= ndim) {
      Printf("Requested axis %d, %d out of range %d", ax1Dproj, axSel, ndim);
      return; 
   }
   Int_t nx, ny, dx, dy;
   CalculatePads(npTbins, nx, ny, dx, dy);
   TCanvas *c1DpjpTrecsel = new TCanvas(Form("c1D%spj%d_Bin%d", hsp->GetName(), ax1Dproj, axSel), Form("Projection 1D %s in %s bins", title.Data(), axSeltitle.Data()), dx, dy);
   c1DpjpTrecsel->Divide(nx, ny);
   
   for(Int_t i=0; i<npTbins; i++){
      if(axSel >= 0){
      	 Int_t binsel[2] = {hsp->GetAxis(axSel)->FindBin(pTlims[i]), hsp->GetAxis(axSel)->FindBin(pTlims[i+1])-1};
      	 if (debug) Printf("SetRange axis %d %d (%f) - %d (%f)", axSel, binsel[0], pTlims[i], binsel[1], pTlims[i+1]);
      	 hsp->GetAxis(axSel)->SetRange(binsel[0], binsel[1]);
      }
      TPaveText *pvt = new TPaveText(0.2, 0.75, 0.5, 0.9,"NDC");
      pvt->SetBorderSize(0);
      pvt->SetFillStyle(0);
      if(pTlims[i+1]>10)
      	 pvt->AddText(Form("%.0f < %s < %.0f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      else if(pTlims[i+1]>1)
      	 pvt->AddText(Form("%.1f < %s < %.1f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      else 
      	 pvt->AddText(Form("%.2f < %s < %.2f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      
      
      c1DpjpTrecsel->cd(i+1);
      if(logy) gPad->SetLogy();
      TH1D *h1 = hsp->Projection(ax1Dproj);
      h1->SetName(Form("%s_%d_Sel%dbin%d", hsp->GetName(), ax1Dproj, axSel, i));
      h1->Draw("colz");
      pvt->Draw();
      if(axSel < 0) break;
   }
   SaveCv(c1DpjpTrecsel);

}

void Compare1DProjectionsinBins(THnSparse *hspB, THnSparse *hspR, Int_t ax1Dproj, TString title, Int_t axSel, TString axSeltitle, Int_t npTbins, Double_t pTlims[], Bool_t logy, Bool_t debug, Double_t normB, Double_t normR, TString legB, TString legR, Bool_t makeRatio){
   ResetAxisSel(hspB);
   ResetAxisSel(hspR);
   Int_t ndim = hspB->GetNdimensions();
   if(ax1Dproj >= ndim || axSel>= ndim) {
      Printf("Requested axis %d, %d out of range %d", ax1Dproj, axSel, ndim);
      return; 
   }
   Int_t nx, ny, dx, dy;
   CalculatePads(npTbins, nx, ny, dx, dy);
   TCanvas *c1DpjpTrecsel = new TCanvas(Form("c1D%spj%d_Bin%d", hspB->GetName(), ax1Dproj, axSel), Form("Projection 1D %s in %s bins", title.Data(), axSeltitle.Data()), dx, dy);
   if(npTbins>1) c1DpjpTrecsel->Divide(nx, ny);
   
   TCanvas *c1DpjpTrecselRatio=0x0;
   if(makeRatio) {
      c1DpjpTrecselRatio = new TCanvas(Form("c1D%spj%d_Bin%dRatio", hspB->GetName(), ax1Dproj, axSel), Form("Ratio Projection 1D %s in %s bins", title.Data(), axSeltitle.Data()), dx, dy);
      if(npTbins>1) c1DpjpTrecselRatio->Divide(nx, ny);
   }
   
   TLegend *leg = new TLegend(0.6, 0.2, 0.9, 0.35);
   for(Int_t i=0; i<npTbins; i++){
      
      Int_t binsel[2] = {hspB->GetAxis(axSel)->FindBin(pTlims[i]), hspB->GetAxis(axSel)->FindBin(pTlims[i+1])-1};
      if (debug) Printf("SetRange axis %d %d (%f) - %d (%f)", axSel, binsel[0], pTlims[i], binsel[1], pTlims[i+1]);
      TPaveText *pvt = new TPaveText(0.2, 0.75, 0.5, 0.9,"NDC");
      pvt->SetBorderSize(0);
      pvt->SetFillStyle(0);
      if(pTlims[i+1]>10)
      	 pvt->AddText(Form("%.0f < %s < %.0f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      else if(pTlims[i+1]>1)
      	 pvt->AddText(Form("%.1f < %s < %.1f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      else 
      	 pvt->AddText(Form("%.2f < %s < %.2f", pTlims[i], axSeltitle.Data(), pTlims[i+1]));
      
      hspB->GetAxis(axSel)->SetRange(binsel[0], binsel[1]);
      TH1D *h1B = hspB->Projection(ax1Dproj);
      h1B->SetName(Form("%s_%d_Sel%dbin%dB", hspB->GetName(), ax1Dproj, axSel, i));
      h1B->SetLineColor(kBlue);
      h1B->SetLineWidth(2);
      Double_t integralB = h1B->Integral();
      if(normB > 0) h1B->Scale(normB);
      else if (integralB > 0) h1B->Scale(1./integralB);
      
      hspR->GetAxis(axSel)->SetRange(binsel[0], binsel[1]);
      TH1D *h1R = hspR->Projection(ax1Dproj);
      h1R->SetName(Form("%s_%d_Sel%dbin%dR", hspR->GetName(), ax1Dproj, axSel, i));
      h1R->SetLineColor(kRed);
      h1R->SetLineWidth(2);
      Double_t integralR = h1R->Integral();
      if(normR > 0) h1R->Scale(normR);
      else if (integralR > 0) h1R->Scale(1./integralR);

 
      c1DpjpTrecsel->cd(i+1);
      if(logy) gPad->SetLogy();
      h1B->Draw();
      h1R->Draw("sames");
      pvt->Draw();
      
      if(makeRatio) {
      	 TH1D *h1Ratio = (TH1D*)h1B->Clone(Form("%sRatio", h1B->GetName()));
      	 h1Ratio->Divide(h1R);
      	 h1Ratio->GetYaxis()->SetRangeUser(0.,5.);
      	 c1DpjpTrecselRatio->cd(i+1);
      	 h1Ratio->Draw("hist");
      }
      
      if (i==0) {
      	 leg->AddEntry(h1B, legB, "L");
      	 leg->AddEntry(h1R, legR, "L");
      }
      if (i == npTbins - 1) {
      	 c1DpjpTrecsel->cd(i+1);
      	 leg->Draw();
      }
   }
   
   SaveCv(c1DpjpTrecsel);
   if(makeRatio) SaveCv(c1DpjpTrecselRatio);
}

void ResetAxisSel(THnSparse *hps){
   for(Int_t i=0; i< hps->GetNdimensions();i++){
      hps->GetAxis(i)->SetRange(0,-1);
   
   }
}
//--------------------------------------------------------------------------------------------------
void plotOverlapInput(){
   
   Int_t nfiles = 2;
   TString paths[nfiles] = {"/data/Work/jets/JetMass/BkgFluctStudies/EmbedSinglePart/ExcludeOverlap/merge/weights/AnalysisResults.root", "/data/Work/jets/JetMass/BkgFluctStudies/EmbedSinglePart/ExcludeOverlap/JetWithSingleTrack/AnalysisResults.root"};
   TString listnames[nfiles] = {"JetShapeConst_Jet_AKTChargedR040_ThrmTracksEmb_pT0150_E_scheme_TCMCMatchnoOvl", "JetShapeConst_Jet_AKTChargedR040_ThrmTracksEmb_pT0150_E_scheme_TCMCMatchnoOvl"};
   TString names[nfiles] = {"PythiaThrm", "PythiaThrmSingle"};
   plotOverlap(nfiles, paths, listnames, names);
}

void plotOverlap(Int_t nfiles, TString paths[], TString listnames[], TString names[]){
   //SetStyle();
   TString h2dname = "fRjetTrvspTj";
   TLegend *leg1D = new TLegend(0.5, 0.7, 0.7, 0.9);
   leg1D->SetBorderSize(0);
   leg1D->SetFillStyle(0);
   
   TCanvas *cproj = new  TCanvas("cproj", "Projections", 600, 1000);
   cproj->Divide(1,2);
   
   for(Int_t f = 0; f<nfiles; f++){
      TLegend *leg2D = new TLegend(0.5, 0.7, 0.7, 0.9);
      leg2D->SetBorderSize(0);
      leg2D->SetFillStyle(0);
      
      TList * list = ReadFile(paths[f], listnames[f]);
      if(!list) continue;
      TCanvas *c2d = new TCanvas(Form("c2D%s", names[f].Data()), "pT jet overlap vs R(jet, emb track)", 800, 800);
      TH2F * h2d = (TH2F*)list->FindObject(h2dname);
      h2d->SetName(Form("%s%s", h2d->GetName(), names[f].Data()));
      c2d->cd();
      h2d->Draw("colz");
      leg2D->AddEntry(h2d, names[f]);
      leg2D->Draw();
      SaveCv(c2d);
      
      TH1F *hprojR = (TH1F*)h2d->ProjectionX(Form("hprojR%s", names[f].Data()));
      TH1F *hprojpTj = (TH1F*)h2d->ProjectionY(Form("hprojpTj%s", names[f].Data()));
      hprojR->SetLineColor(colors[f]);
      hprojpTj->SetLineColor(colors[f]);
      
      leg1D->AddEntry(hprojR, names[f]);
      
      cproj->cd(1);
      if(f==0) hprojR->Draw();
      else hprojR->Draw("sames");
      leg1D->Draw();
      cproj->cd(2);
      gPad->SetLogy();
      if(f==0) hprojpTj->Draw();
      else hprojpTj->Draw("sames");
      
      
   }
   SaveCv(cproj);

}

//________________________________________________________________________________________________

void CompareRandomConeAndSingleTrackEmbeddingFixedInput(){
   //
   // draws comparison plots between the single track embedding and the random cone method
   // modify input strings as needed
   //
   
   TString filename = "AnalysisResults.root";
   TString dirnameRC= "BackFlucRandomCone";
   TString dirnameConst = //"JetShapeConst_JetEmb_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC";
   //"JetShapeConst_Jet_AKTChargedR040_ThrmTracksEmb_pT0150_pt_scheme_TCMCMatch";
   "JetShapeConst_Jet_AKTChargedR040_ThrmTracksEmbSingle_pT0150_E_scheme_TCMCMatch";
   CompareRandomConeAndSingleTrackEmbedding(filename, dirnameRC, dirnameConst);
}

void CompareRandomConeAndSingleTrackEmbedding(TString filename, TString dirnameRC, TString dirnameConst){
   
   SetStyle(0);
   TList *lRC=ReadFile(filename, dirnameRC);
   TList *lConst=ReadFile(filename, dirnameConst);
   if(!lRC && !lConst) {
      Printf("Input not found");
      return;  
   }

   //histograms in the random cone output to be drawn
   Int_t nRC = 1;
   TString h2namesRC[nRC] = {"fhEtaPhiConeCentre"};
   for(Int_t i=0; i<nRC; i++){
      TH2F* h = (TH2F*) lRC->FindObject(h2namesRC[i]);
      TCanvas *c = new TCanvas(Form("c%s",h2namesRC[i].Data()), h2namesRC[i], 800, 800);
      c->cd();
      h->Draw("colz");
      
      TPaveText *pvT = new TPaveText(0.4, 0.3, 0.6, 0.5, "NDC");
      pvT->SetFillStyle(0);
      pvT->SetBorderSize(0);
      pvT->AddText(Form("N = %.0f", h->GetEntries()));
      c->cd(); pvT->DrawClone();
      SaveCv(c);
   }
   
   //histograms in the single track embedding output to be drawn
   Int_t nConst = 1;
   TString h1namesConst[nConst] = {"fhAreaJet"};
   for(Int_t i=0; i<nConst; i++){
      TH2F* h = (TH2F*) lConst->FindObject(h1namesConst[i]);
      h->SetLineColor(kBlue);
      h->SetLineWidth(3);
      TCanvas *c = new TCanvas(Form("c%s",h1namesConst[i].Data()), h1namesConst[i], 800, 800);
      c->cd();
      h->Draw("colz");
      
      TPaveText *pvT = new TPaveText(0.3, 0.3, 0.7, 0.5, "NDC");
      pvT->SetFillStyle(0);
      pvT->SetBorderSize(0);
      pvT->AddText(Form("N = %.0f, #mu = %.3f, #sigma = %.3f", h->GetEntries(), h->GetMean(), h->GetRMS()));
      c->cd(); pvT->DrawClone();
      SaveCv(c);
   }
   
   //list of histograms to be compared
   Int_t nhists = 2;
   TString hnamesRCcmp[nhists] = {"fhMasspTInCone", "fhNConstituents"};
   TString hnamesConstcmp[nhists] = {"fhptjetSMinusSingleTrack", "fhNconstit"};
   TString title[nhists] = {"#it{p}_{T} distribution", "N constituents"};
   
   for(Int_t i=0; i<nhists; i++){
      
      Printf("Reading %s and %s", hnamesRCcmp[i].Data(), hnamesConstcmp[i].Data());
      TLegend *leg = new TLegend(0.35, 0.73, 0.9, 0.88, title[i]);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      
      TH1F* hRC = 0x0;
      if(i == 0){
      	 TH2F *h2tmp = (TH2F*)lRC->FindObject(hnamesRCcmp[i]);
      	 hRC = (TH1F*) h2tmp->ProjectionY("fhpTInCone");
      } else {
      	 hRC = (TH1F*) lRC->FindObject(hnamesRCcmp[i]);
      }
      if(!hRC) {
      	 lRC->ls();
      	 continue;
      }
      hRC->SetLineColor(kGreen+4);
      hRC->SetLineWidth(3);
      leg->AddEntry(hRC, Form("Random Cone (N = %.0f, #mu = %.3f, #sigma = %.3f)", hRC->GetEntries(), hRC->GetMean(), hRC->GetRMS()), "L");
      
      Double_t maximum = hRC->GetMaximum();
      
      TH1F* hConst = 0x0;
      hConst = (TH1F*) lConst->FindObject(hnamesConstcmp[i]);
      if(!hConst) {
      	 continue;
      }
      
      hConst->SetLineColor(kBlue);
      hConst->SetLineWidth(3);
      
      if(hConst->GetMaximum() > maximum) maximum = hConst->GetMaximum();
      
      leg->AddEntry(hConst, Form("Single-trk Emb (N = %.0f, #mu = %.3f, #sigma = %.3f)", hConst->GetEntries(), hConst->GetMean(), hConst->GetRMS()), "L");
      
      Printf("Drawing");
      TCanvas *c = new TCanvas(Form("cOvl%d",i), Form("Overlap %s", title[i].Data()), 800, 800);
      c->cd();
      if(i == 0 ) gPad->SetLogy();
      hRC->SetMaximum(maximum);
      
      hRC->Draw();
      hConst->Draw("sames");
      leg->Draw();
      SaveCv(c);
   }
   
}

/*
//----------------------------------------------------------------------------------------
void PerformConvolution2D(TString inputFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root",  TString listname = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged", TString correctionFile = "MassFromBkgFluctuations.root"){
   SetStyle(0);
   TList *l=ReadFile(inputFile, listname);

   TString hname = "fhnMassResponse";
   
   THnSparseF *hsp = (THnSparseF*) l->FindObject(hname);
   if(!hsp){
      Printf("THnSparse not found");
      return;
   }
   TFile *f = new TFile(correctionFile);
   TString namefout = "ModifiedMassBgkFluct.root";
   TFile *fout = new TFile(namefout, "recreate");
   Printf("This macro writes an output in a ROOT file (%s).", namefout.Data());
   Double_t pTlims[4] = {40., 60., 80., 100.};
   
   
   for(Int_t ipt = 0; ipt < 3; ipt++){
      TH1F *h = (TH1F*)f->Get(Form("hMassSubpT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]));
      if(!h){
      	 Printf("%d) not found", ipt);
      	 continue;
      }
      Int_t binpt[2] = {hsp->GetAxis(3)->FindBin(pTlims[ipt]) , hsp->GetAxis(3)->FindBin(pTlims[ipt+1])-1 }; //3 = pt particle level
      hsp->GetAxis(3)->SetRange(binpt[0], binpt[1]);
      
      TH1F* hMass = (TH1F*)hsp->Projection(1); //1 = mass particle level
      hMass->SetName(Form("hMasspT%.0f-%.0f", pTlims[ipt], pTlims[ipt+1]));
      Printf("Projection %d", ipt);
      hMass->SetMarkerStyle(20);
      
      TH1F* hMassNew = Convolution2D(hMass, h);
      fout->cd();
      hMassNew->Write();
   }
   
}
*/
//--------------------------------------------------------------------------------------------------

void NumberOfExcludedJets(TString filename = "AnalysisResults.root", TString listname = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TCnoOvl"){
   
   // draw the histogram of the distrnace between the single particle embed jet and a signal jet and integrate the number of jets excluded or accepted with respect to the total
   // compare to the number of jets filling the THnSparse (should be the same as the integral above R)
   
   // Conclusions: the number of entries in the jetO histograms is much smaller because in many case there is no signal jet at all
   // Still to be understood why the number of events with zero signal jets + the number of events with no overlap with a signal jet is larger than the number of entried in the THnSparse. They're similar though...
   
   TList * list = ReadFile(filename, listname);
   if(!list){
      Printf("Exit");
      return;
   }
   
   TString spname = "fhnDeltaMass_0";
   THnSparse *hsp = (THnSparse*)list->FindObject(spname);
   if(!hsp){
      Printf("%s  not found", spname.Data());
      return;
   }
   TH1D* hpT = (TH1D*)hsp->Projection(3); //pT true (embedded track)
   Float_t entriesProj0 = hpT->GetEntries();
   
   
   TString hname = "fhRjetTrvspTj";
   
   TH2F *h2deltaRvsptSjet = (TH2F*) list->FindObject(hname);
   
   if(!h2deltaRvsptSjet) {
      Printf("%s not found", hname.Data());
      return;
   }
   hname = "fh2MSubMatch_0";
   
   TH2F *hMatch = (TH2F*) list->FindObject(hname);
   
   if(!hMatch) {
      Printf("%s not found", hname.Data());
      return;
   }
   
   hname = "fhNJetsSelEv";
   TH1F *hNSigJetEv = (TH1F*) list->FindObject(hname);
   
   if(!hNSigJetEv) {
      Printf("%s not found", hname.Data());
      return;
   }
   
   Double_t R = 0.4;
   Int_t binSel = h2deltaRvsptSjet->FindBin(R, 0, 0), binSelX, y, z;
   h2deltaRvsptSjet->GetBinXYZ(binSel, binSelX, y, z);
   Printf("Overlap < bin %d ", binSelX);
   TH1D *hpTSjOverlap   = (TH1D*)h2deltaRvsptSjet->ProjectionY("hpTSjOverlap", 1, binSelX);
   TH1D *hpTSjNoOverlap = (TH1D*)h2deltaRvsptSjet->ProjectionY("hpTSjNoOverlap", binSelX, h2deltaRvsptSjet->GetXaxis()->GetNbins());
   hpTSjOverlap->  SetLineColor(kRed-3);
   hpTSjNoOverlap->SetLineColor(kGreen+3);
   hpTSjOverlap->  SetLineWidth(3);
   hpTSjNoOverlap->SetLineWidth(3);
   
   TPaveText *pvtxt1 = new TPaveText(0.1, 0.8, 0.8, 0.9, "NDC");
   pvtxt1->SetFillStyle(0);
   pvtxt1->SetBorderSize(0);
   pvtxt1->SetTextColor(kRed-3);
   pvtxt1->AddText(Form("Integral Overlap single Tr - signal jet = %.4e", hpTSjOverlap->Integral()));
   
   TPaveText *pvtxt2 = new TPaveText(0.1, 0.7, 0.8, 0.8, "NDC");
   pvtxt2->SetFillStyle(0);
   pvtxt2->SetBorderSize(0);
   pvtxt2->SetTextColor(kGreen+3);
   pvtxt2->AddText(Form("Integral NO Overlap single Tr - signal jet = %.4e", hpTSjNoOverlap->Integral()));
   
   TCanvas *cpTSigJ = new TCanvas("cpTSigJ", "pT signal jets", 600, 600);
   cpTSigJ->cd();
   gPad->SetLogy();
   hpTSjOverlap->Draw();  
   hpTSjNoOverlap->Draw("sames");
   pvtxt1->Draw();
   pvtxt2->Draw();
   
   SaveCv(cpTSigJ);
   
   // Compare Integral
   
   TH1F *hIntegrals = new TH1F("hIntegrals", "THnSparse (0), Matched jets (1), No Sig jets (2), No Overlap (3), Overlap (4), Must match 0&1 (5) ; ;Number of entries", 6, -0.5, 5.5);
   hIntegrals->SetLineWidth(3);
   hIntegrals->SetLineColor(kBlue+1);
   hIntegrals->Fill(3., hpTSjNoOverlap->Integral());
   hIntegrals->Fill(4., hpTSjOverlap->Integral());
   
   
   hIntegrals->Fill(0., entriesProj0);
   
   TH1F* hMatch1D = (TH1F*)hMatch->ProjectionY();
   
   hIntegrals->Fill(1., hMatch1D->GetBinContent(2));
   
   hIntegrals->Fill(2., hNSigJetEv->GetBinContent(1));
   
   hIntegrals->Fill(5., hNSigJetEv->GetBinContent(1) + hpTSjNoOverlap->Integral());

   TCanvas *cIntrgrals = new TCanvas("cIntrgrals", "Number of entries No Overlap", 600, 600);
   
   cIntrgrals->cd();
   hIntegrals->Draw("htext0");
   
}

//--------------------------------------------------------------------------------------------------

void DrawOtherHistos(TString filename = "AnalysisResults.root", TString listname = "JetShapeDeriv_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_TC"){
   
   TList * list = ReadFile(filename, listname);
   Int_t nentries = list ->GetEntries();
   
   for(Int_t i=0; i <nentries; i++){
      TClass* objtype=list->At(i)->IsA();
      TString tpname=objtype->GetName();
      
      if(tpname=="TH2F"){
      	 TH2F* h=(TH2F*)list->At(i);
      	 
      	 TCanvas *c = new TCanvas(Form("c%s", h->GetName()), Form("c%s", h->GetName()));
      	 h->Draw("colz");
      }
   }
}

//--------------------------------------------------------------------------------------------------

void WriteMPYTHIAinpTRangeFormEmbeddingInput(Double_t ppartmin = 10, Double_t ppartmax = 140, TString inputFile = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train576/AnalysisResultsWeighted.root", TString inputList = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_schemeMerged", TString thnspname = "fhnMassResponse", Int_t axproj = 1, Int_t axsel = 3, Bool_t drawF = kTRUE, Bool_t twoD = kTRUE, Int_t axproj2 = 3){
   
   TList * list = ReadFile(inputFile, inputList);
   if(!list) {
      return;
   }
   THnSparseF *hsp = (THnSparseF*) list->FindObject(thnspname);
   if(!hsp) {
      Printf("%s not found, exit", thnspname.Data());
      list->ls();
      return;
      
   }
   
   TH2D *h2CorRec = hsp->Projection(0, 2);
   TH2D *h2CorGen = hsp->Projection(1, 3);
   TH2D *h2ResM   = hsp->Projection(0, 1);
   TH2D *h2RespT  = hsp->Projection(2, 3);
   TCanvas *cRespCorr = new TCanvas("cRespCorr", "Response and pT/M correlations", 1000, 1000);
   cRespCorr->Divide(2,2);
   cRespCorr->cd(1);
   h2CorRec->Draw("colz");
   cRespCorr->cd(2);
   h2CorGen->Draw("colz");
   cRespCorr->cd(3);
   h2ResM->Draw("colz");
   cRespCorr->cd(4);
   h2RespT->Draw("colz");
   
   hsp->GetAxis(axsel)->SetRangeUser(ppartmin, ppartmax);
   TH1D *houtput = hsp->Projection(axproj);
   houtput->SetName(Form("hDProj%dRange%d_%.0f_%.0f", axproj, axsel, ppartmin, ppartmax));
   TH1F *houtputf = ConvertTH1DinF(houtput);
   houtputf->SetName(Form("hProj%dRange%d_%.0f_%.0f", axproj, axsel, ppartmin, ppartmax));
   //hsp->GetAxis(axsel)->SetRangeUser(0, -1);
   TH1D *houtselax = hsp->Projection(axsel);
   houtselax->SetName(Form("hDProj%dRange%d_%.0f_%.0f", axsel, axsel, ppartmin, ppartmax));
   TH1F *houtselaxf = ConvertTH1DinF(houtselax);
   houtselaxf->SetName(Form("hProj%dRange%d_%.0f_%.0f", axsel, axsel, ppartmin, ppartmax));
   TH2D *h2doutp = hsp->Projection(axproj, axproj2);
   h2doutp->SetName(Form("hDProj%dvs%dRange%d_%.0f_%.0f", axproj, axproj2, axsel, ppartmin, ppartmax));
   TH2F *h2doutpf = ConvertTH2DinF(h2doutp);
   h2doutpf->SetName(Form("hProj%dvs%dRange%d_%.0f_%.0f", axproj, axproj2, axsel, ppartmin, ppartmax));
   
   TCanvas *cMass = new TCanvas(Form("c%s_proj%d_Sel%d_%.0f%.0f", thnspname.Data(), axproj, axsel, ppartmin, ppartmax), "Mass");
   TPaveText *pvt = new TPaveText(0.2, 0.75, 0.5, 0.9,"NDC");
   pvt->SetBorderSize(0);
   pvt->SetFillStyle(0);
   pvt->AddText(Form("%.0f < %s < %.0f", ppartmin, hsp->GetAxis(axsel)->GetTitle(), ppartmax));

   cMass->cd();
   gPad->SetLogy();
   if(drawF) houtputf->Draw();
   else houtput->Draw();
   pvt->Draw();
   
   TCanvas *cpt = new TCanvas(Form("c%s_proj%d", thnspname.Data(), axsel), "pT");
   cpt->cd();
   gPad->SetLogy();
   if(drawF) houtselaxf->Draw();
   else houtselax->Draw();
   
   TCanvas *cmasspt = new TCanvas(Form("c%s_proj%dvs%d", thnspname.Data(), axproj, axproj2));
   cmasspt->cd();
   //gPad->SetLogz();
   if(drawF) h2doutpf->Draw("colz");
   else h2doutp->Draw("colz");
  
   TString filename = Form("MassEmbInput_%s_proj%d_Sel%d_%.0f%.0f.root", thnspname.Data(), axproj, axsel, ppartmin, ppartmax);
   if(twoD) filename = Form("MassEmbInput_%s_proj%dvs%d_Sel%d_%.0f%.0f.root", thnspname.Data(), axproj, axproj2, axsel, ppartmin, ppartmax);
   TFile *fout = new TFile(filename, "recreate");
   houtputf->Write();
   houtselaxf->Write();
   h2doutpf->Write();
   
}


void DrawUnfoldingMatrix(Int_t type = 1, Float_t minPtDet = 8, Float_t maxPtDet = 98, TString fileName = "UnfoldingThnSparse.root", TString suff = ""){
   // type 1 = read from UnfoldingThnSparse%s.root (coming from MatrixDeltaMDeltapT)
   // type 2 = read from UnfoldingMatrix.root
   TFile *f = new TFile(fileName);
   if(!f->IsOpen()){
      Printf("File %s not found", fileName.Data());
      return;
   
   }
   TH2D *h2dM  = 0x0;
   TH2D *h2dpt = 0x0;
   TH2D *hptMD = 0x0;
   TH2D *hptMP = 0x0;
   
   TH1D *hMD   = 0x0;
   TH1D *hMP   = 0x0;
   TH1D *hpTD  = 0x0;
   TH1D *hpTP  = 0x0;
   
   if(type == 1){
      TString hnamesp =  "fhnDeltaMass_0_proj_2_3_4_5";
      THnSparse *hsp = dynamic_cast<THnSparse*>(f->Get(hnamesp));
      if(!hsp) {
      	 Printf("THnSparse %s not found", hnamesp.Data());
      	 f->ls();
      	 return;
      }
      
      for(Int_t ix = 0; ix < hsp->GetNdimensions(); ix++){
      	 Printf("Axis %d -> title %s", ix, hsp->GetAxis(ix)->GetTitle());
      }
      
      Printf("Applying the following pT det cut: %.1f - %.1f", minPtDet, maxPtDet);
      hsp->GetAxis(2)->SetRange( hsp->GetAxis(2)->FindBin(minPtDet), hsp->GetAxis(2)->FindBin(maxPtDet));
      
      h2dM  = hsp->Projection(0,1);
      h2dpt = hsp->Projection(2,3);
      hptMD = hsp->Projection(0,2);
      hptMP = hsp->Projection(1,3);
      
      hMD   = hsp->Projection(0);
      hMP   = hsp->Projection(1);
      hpTD  = hsp->Projection(2);
      hpTP  = hsp->Projection(3);
      
      
   }
   
   if(type == 2){
      TString hnamesp =  "fhnMassResponse";
      THnSparse *hsp = dynamic_cast<THnSparse*>(f->Get(hnamesp));
      if(!hsp) {
      	 Printf("THnSparse %s not found", hnamesp.Data());
      	 f->ls();
      	 return;
      }
      for(Int_t ix = 0; ix < hsp->GetNdimensions(); ix++){
      	 Printf("Axis %d -> title %s", ix, hsp->GetAxis(ix)->GetTitle());
      }
      h2dM  = hsp->Projection(0,1);
      h2dpt = hsp->Projection(2,3);
      hptMD = hsp->Projection(0,2);
      hptMP = hsp->Projection(1,3);
      
      hMD   = hsp->Projection(0);
      hMP   = hsp->Projection(1);
      hpTD  = hsp->Projection(2);
      hpTP  = hsp->Projection(3);
      
   }
   hMD ->Sumw2();
   hMP ->Sumw2();
   hpTD->Sumw2();
   hpTP->Sumw2();
   hMD ->SetName(Form("hMD%s", suff.Data()));
   hMP ->SetName(Form("hMP%s", suff.Data()));
   hpTD->SetName(Form("hpTD%s", suff.Data()));
   hpTP->SetName(Form("hpTP%s", suff.Data()));
   hMD ->SetLineWidth(2);
   hMP ->SetLineWidth(2);
   hpTD->SetLineWidth(2);
   hpTP->SetLineWidth(2);
   hMD ->SetLineColor(colors[0]);
   hMP ->SetLineColor(colors[1]);
   hpTD->SetLineColor(colors[0]);
   hpTP->SetLineColor(colors[1]);
   hMD ->SetMarkerColor(colors[0]);
   hMP ->SetMarkerColor(colors[1]);
   hpTD->SetMarkerColor(colors[0]);
   hpTP->SetMarkerColor(colors[1]);
   hMD ->SetMarkerStyle(20);
   hMP ->SetMarkerStyle(24);
   hpTD->SetMarkerStyle(20);
   hpTP->SetMarkerStyle(24);
   
   h2dM ->SetName(Form("h2dM%s", suff.Data()));
   h2dpt->SetName(Form("h2dpt%s", suff.Data()));
   hptMD->SetName(Form("hptMD%s", suff.Data()));
   hptMP->SetName(Form("hptMP%s", suff.Data()));
   
   // 2D projections
   TCanvas *cResponse = new TCanvas(Form("cResponse%.0f-%.0f", minPtDet, maxPtDet), "Response", 900, 400);
   cResponse->Divide(2,1);
   cResponse->cd(1);
   h2dM->Draw("colz");
   cResponse->cd(2);
   h2dpt->Draw("colz");
   
   TCanvas *cCorrelation = new TCanvas(Form("cCorrelation%.0f-%.0f", minPtDet, maxPtDet), "M vs pT correlations", 900, 400);
   cCorrelation->Divide(2,1);
   cCorrelation->cd(1);
   
   hptMD->Draw("colz");
   cCorrelation->cd(2);
   hptMP->Draw("colz");
   
   // 1D projections
   
   TCanvas *cpT = new TCanvas(Form("cpT%.0f-%.0f", minPtDet, maxPtDet), "pT", 400, 400);
   cpT->cd();
   gPad->SetLogy();
   hpTP->Draw();
   hpTD->Draw("sames");
   
   TCanvas *cM = new TCanvas(Form("cM%.0f-%.0f", minPtDet, maxPtDet), "pM", 400, 400);
   cM->cd();
   gPad->SetLogy();
   hMD->Draw();
   hMP->Draw("sames");
   
   SaveCv(cResponse, suff);
   SaveCv(cCorrelation, suff);
   SaveCv(cpT, suff);
   SaveCv(cM, suff);
   
   TFile *fout = new TFile(Form("SaveProjections%.0f-%.0f%s.root", minPtDet, maxPtDet, suff.Data()), "recreate");
   hMD ->Write();
   hMP ->Write();
   hpTD->Write();
   hpTP->Write();
   h2dM ->Write();
   h2dpt->Write();
   hptMD->Write();
   hptMP->Write();
   
}

void CompareFinalWithEmbeddedDistr(TString saveprojfilesuff = "",  TString suff = "", TString fembname = "AnalysisResults.root", TString listname = "SingleTrackEmbedding", Bool_t normalize = kFALSE){
   
   //read file created in DrawUnfoldingMatrix and compare to the embedded distribution from the analysis output
   TString fprojname = Form("SaveProjections%s%s.root", saveprojfilesuff.Data(), suff.Data());
   TFile* fpr = new TFile(fprojname);
   if(!fpr->IsOpen()){
      Printf("File %s not found, did you run DrawUnfoldingMatrix?", fprojname.Data());
      return;
   
   }
   fpr->ls();
   
   const Int_t n = 4;
   TString namesFinal[n] = {Form("hMD%s", suff.Data()), Form("hMP%s", suff.Data()), Form("hpTD%s", suff.Data()), Form("hpTP%s", suff.Data())};
   
   //read general output file with list from embedding task
   TList * list = ReadFile(fembname, listname);
   if(!list){
      Printf("%s not found, please check", listname.Data());
      return;
   }
   list->ls();
   const Int_t ne = 4;
   TString namesEmbed[ne] = {"fhMEmb", "fhpTEmb", "fhEtaEmb", "fhPhiEmb"};
   
   TLegend *leg = new TLegend(0.5, 0.4, 0.8, 0.7);
   leg->SetFillStyle(0);
   
   //read histograms
   TH1D *hproj[n];
   for(Int_t i = 0; i<n ; i++){
      hproj[i] = 0x0;
      hproj[i] = dynamic_cast<TH1D*>(fpr->Get(namesFinal[i]));
      if(!hproj[i]) Printf("%s not found", namesFinal[i].Data());
      else {
      	 if((i < 2)){
      	    if((i % 2)) leg->AddEntry(hproj[i], "Particle", "PL");
      	    else  leg->AddEntry(hproj[i], "Detector", "PL");
      	 }
      }
   }
   
   TH1F* hEmb[ne];
   for(Int_t i = 0; i<ne ; i++){
      hEmb[i] = 0x0;
      hEmb[i] = dynamic_cast<TH1F*>(list->FindObject(namesEmbed[i]));
      if(hEmb[i]) hEmb[i]->SetMarkerStyle(25);
      else Printf("%s not found", namesEmbed[i].Data());
   }
   Int_t id = 0;
   leg->AddEntry(hEmb[id], Form("Embedded%s", normalize ? " N" : ""), "PL");
   // 1D projections
   
   TCanvas *cpT = new TCanvas(Form("cpT%s", saveprojfilesuff.Data()), "pT", 600, 600);
   cpT->cd();
   gPad->SetLogy();
   if(hproj[2]) {
      hproj[2]->GetXaxis()->SetRangeUser(-10, 150);
      hproj[2]->Draw();
   }
   if(hproj[3]) hproj[3]->Draw("sames");
   if(hEmb[1] ) hEmb[1] ->Draw("sames");
   leg->Draw();
   
   if(normalize){
      //try to normalize the embedded distribution to the others
      for(Int_t i = 0; i<ne ; i++){
      	 Double_t contentEmb = hEmb[i]->GetBinContent(4);
      	 Double_t contentDet = hproj[0]->GetBinContent(4);
      	 Printf(" %.3f/ %.3f", contentEmb, contentDet);
      	 
      	 Double_t scale =  hproj[0]->Integral(); //contentDet/contentEmb *
      	 hEmb[i]->Scale(scale / hEmb[i]->Integral());
      	 
      }
   }
   TCanvas *cM = new TCanvas(Form("cM%s", saveprojfilesuff.Data()), "pM", 600, 600);
   cM->cd();
   gPad->SetLogy();
   if(hproj[0]) {
      hproj[0]->GetXaxis()->SetRangeUser(-1, 16);
      hproj[0]->Draw();
      
   }
   if(hproj[1]) hproj[1]->Draw("sames");
   if(hEmb[0] ) hEmb[0] ->Draw("sames");
   leg->Draw();
   
   TCanvas *cEta = new TCanvas(Form("cEta%s", saveprojfilesuff.Data()), "Eta distribution embedded tracks", 600, 600);
   cEta->cd();
   if(hEmb[2] ) hEmb[2] ->Draw();
   TCanvas *cPhi = new TCanvas(Form("cPhi%s", saveprojfilesuff.Data()), "Phi distribution embedded tracks", 600, 600);
   cPhi->cd();
   if(hEmb[3] ) hEmb[3] ->Draw();
   
   SaveCv(cpT, suff);
   SaveCv(cM, suff);
   SaveCv(cEta, suff);
   SaveCv(cPhi, suff);
}
