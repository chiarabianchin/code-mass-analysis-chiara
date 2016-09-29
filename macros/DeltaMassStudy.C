#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH3F.h>
#include <TF1.h>

#include </data/Work/MyCodeJetMass/utils/CommonTools.C>

void Make3Projections(TH3F* h3, TH1D*& hX, TH1D*& hY, TH1D*& hZ);

void DeltaMassStudy(Int_t nevents = 1000, Bool_t zeroM = kTRUE){

   TRandom3 *rdm = new TRandom3(1234);
   
   Double_t ptlims[2]  = {40, 100.}; //use uniform
   Double_t etalims[2] = {-0.8, 0.8}; //use uniform
   Double_t philims[2] = {0, TMath::TwoPi()}; //use uniform
   Double_t maslims[2] = {0, 7}; //use uniform
   
   Double_t dptlims[2]  = {-40, 40.}; //use uniform
   Double_t detalims[2] = {-0.3, 0.3}; //use uniform
   Double_t dphilims[2] = {-1, 1}; //use uniform
   Double_t dmaslims[2] = {-5, 5}; //use uniform

   //TF1
   TF1 *fbkgd = new TF1("fbkgd", "[0]*TMath::Power([1], 2)*x*TMath::Exp(-[1]*x)", 0.15, 200.);
   fbkgd->SetParNames("Amplitude", "b (GeV/c)^{-1}");
   fbkgd->SetParameters(1, 2./20.5);
   
   TCanvas *cfunc = new TCanvas("cfunc", "Func");
   cfunc->cd();
   fbkgd->Draw();
   
   TF1 *fgaus = new TF1("fgaus","gaus(0)", 30, 120);
   //fgaus->SetParameter(1, 40);
   TCanvas *cgaus = new TCanvas("cgaus", "Gaussian");
   cgaus->cd();
   fgaus->Draw();
   
   //histograms
   TH3F *hLorentzGen = new TH3F("hLorentzGen", "; #it{p}_{T} (GeV/#it{c}); #eta; #phi", 100, ptlims[0], ptlims[1], 100, etalims[0], etalims[1], 100, philims[0], philims[1]);
   TH1F *hMGen = new TH1F("hMGen", "Mass", 100, maslims[0], maslims[1]);
   
   TH3F *hLorentzRec = (TH3F*)hLorentzGen->Clone("hLorentzRec");
   TH1F *hMRec = (TH1F*)hMGen->Clone("hMRec");
   
   TH3F *hDeltaLorentz = new TH3F("hDeltaLorentz", "; #it{p}_{T}^{rec} - #it{p}_{T}^{gen} (GeV/#it{c}); #eta^{rec} - #eta^{gen}; #phi^{rec} - #phi^{gen}", 100, dptlims[0], dptlims[1], 100, detalims[0], detalims[1], 100, dphilims[0], dphilims[1]);
   TH1F *hDeltaM = new TH1F("hDeltaM", "Delta Mass", 100, dmaslims[0], dmaslims[1]);
   Printf("N events = %d", nevents);
   for(Int_t ie = 0; ie < nevents ; ie++){
      TLorentzVector vecGen;
      Double_t pt, eta, phi, m = 0;
      pt  = 0; //fbkgd->GetRandom(ptlims[0], ptlims[1]); //rdm->Uniform(ptlims[0], ptlims[1]);
      while (pt == 0) pt  = fbkgd->GetRandom(ptlims[0], ptlims[1]);
      //Printf("pT = %f", pt);
      eta = rdm->Uniform(etalims[0], etalims[1]);
      phi = rdm->Uniform(philims[0], philims[1]);
      if(!zeroM) m = rdm->Gaus((maslims[1]-maslims[0])*0.5);
      vecGen.SetPtEtaPhiM(pt, eta, phi, m);
      hLorentzGen->Fill(pt, eta, phi);
      hMGen->Fill(m);
      
      //gaussian smearing to mock reconstruction
      TLorentzVector vecRec;
      Double_t ptR, etaR, phiR, mR;
      ptR  = rdm->Gaus(pt, 0.1);
      etaR = rdm->Gaus(eta,0.01);
      Printf("Eta %f gives etaR %f", eta, etaR);
      phiR = rdm->Gaus(phi, 0.1);
      mR   = rdm->Gaus(m, 0.1);
      vecRec.SetPtEtaPhiM(ptR, etaR, phiR, mR);
      hLorentzRec->Fill(ptR, etaR, phiR);
      hMRec->Fill(mR);
      
      TLorentzVector vecDelta = vecRec - vecGen; //(TLorentzVector*)vecRec->Clone();
      //Printf("Delta phi V1 = %f, V2 = %f", phiR-phi, vecDelta.Phi());
      hDeltaLorentz->Fill(vecDelta.Pt(), vecDelta.Eta(), vecDelta.Phi());
      hDeltaM->Fill(vecDelta.M());
      
      
   }

   TCanvas* cGeneratedpj = new TCanvas("cGeneratedpj", "Projection Generated", 900, 900);
   cGeneratedpj->Divide(2,2);
   
   TH1D *hpTGen = 0x0, *hEtaGen = 0x0, *hPhiGen = 0x0;
   
   Make3Projections(hLorentzGen, hpTGen, hEtaGen, hPhiGen);
   
   TH1D *hpTRec = 0x0, *hEtaRec = 0x0, *hPhiRec = 0x0;
   
   Make3Projections(hLorentzRec, hpTRec, hEtaRec, hPhiRec);
   hpTRec ->SetLineColor(kRed);
   hEtaRec->SetLineColor(kRed);
   hPhiRec->SetLineColor(kRed);
   hMRec  ->SetLineColor(kRed);
   
   cGeneratedpj->cd(1);
   hpTGen->Draw();
   hpTRec->Draw("sames");
   
   cGeneratedpj->cd(2);
   hEtaGen->Draw();
   hEtaRec->Draw("sames");
   
   cGeneratedpj->cd(3);
   hPhiGen->Draw();
   hPhiRec->Draw("sames");
   
   cGeneratedpj->cd(4);
   hMGen->Draw();
   hMRec->Draw("sames");
   
   TCanvas* cDeltapj = new TCanvas("cDeltapj", "Projection Difference", 900, 900);
   cDeltapj->Divide(2,2);
   
   TH1D *hDpT = 0x0, *hDEta = 0x0, *hDPhi = 0x0;
   
   Make3Projections(hDeltaLorentz, hDpT, hDEta, hDPhi);
   
   cDeltapj->cd(1);
   hDpT->Draw();
   
   cDeltapj->cd(2);
   hDEta->Draw();
   
   cDeltapj->cd(3);
   hDPhi->Draw();
   
   cDeltapj->cd(4);
   hDeltaM->Draw();
   
   
}

void Make3Projections(TH3F* h3, TH1D*& hX, TH1D*& hY, TH1D*& hZ){
   hX = h3->ProjectionX();
   hX->SetName(Form("%s_X",h3->GetName()));
   hY = h3->ProjectionY();
   hY->SetName(Form("%s_Y",h3->GetName()));
   hZ = h3->ProjectionZ();
   hZ->SetName(Form("%s_Z",h3->GetName()));
   
}
