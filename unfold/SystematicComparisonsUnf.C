#ifndef SystematicComparisonsUnf_C
#define SystematicComparisonsUnf_C
#include <TH3.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TClass.h>
#include <TMath.h>
#include </data/macros/LoadALICEFigures.C>
#include </data/Work/MyCodeJetMass/utils/CommonTools.C>

//global
//binning
const Int_t nptbins = 3;//4;
Double_t ptlims[nptbins+1] = {60., 80., 100., 120};//40., 
Double_t maxRangeMass[nptbins] = {26., 26., 26.};//20., 
Double_t maxRangeMassFinal[nptbins] = {18., 22., 24.};//20., 

TH1D** CompareResults(const Int_t ninputs, TString files[], TString hnamebase[], TString legs[], Int_t offset[], Bool_t changeColor = kFALSE, Bool_t writeRatios = kFALSE, Bool_t noUniform = kFALSE, Int_t base = 1, TString name = "setit", Bool_t logscale = kFALSE);
TH1D* SystematicFromRatio(TH1D* hRatio, TString systname);
TH1D* SetMassValueInSystematic(TH1D* hSyst, TH1D* hMass);
TList* ReadHistogramsInFile(TString infilename, TString basehname, Bool_t usebinlims, Int_t n);
TList* AddSystematicstoMassFromFile(TString inMfilename, TString hbasenamemass, Bool_t usebinlimsForM,  TString inSysfilename, TString hbasenamesys, Bool_t usebinlimsForSys);
TH1D* AddInQuadrature(TH1D** haddq, Int_t nh, Int_t index = -1, TString name = "");
TH1D* MaxSpread(TH1D** hsystratio, Int_t nh, Int_t index, TString name);
void MaxSpread(const Int_t nval, Double_t values[], Int_t idbase, Double_t& up, Double_t& down);
TH1D* SmoothUncertaintyByAverage(const Int_t nh, TH1D** hinputUnc, const char* newname);
TH1D** SmoothUnceraintyByFit(const Int_t nh, TH1D** hinputUnc, Int_t fittype, const char* newname, TH1D** hSysOnError, Double_t minx= -1, Double_t maxx = -1);
TH1D** GetMeanSyst(const Int_t ninputs, TString filenames[], TString hnameinput[]);
TH1D* SumMeanSyst(const Int_t ninputs, TH1D** hsysPart);

//_________________________________________________________________________________________________
//method implementations
TH1D** CompareResults(const Int_t ninputs, TString files[], TString hnamebase[], TString legs[], Int_t offset[], Bool_t changeColor, Bool_t writeRatios, Bool_t noUniform, Int_t base, TString name, Bool_t logscale){
   // the pt bin number will be appended to the hnamebase
   
   if(logscale) Printf("Draw in log scale");
   TH1 *histogram[ninputs][nptbins];
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   TCanvas *cMass = new TCanvas(Form("cMass"), Form("Mass"), dx, dy);
   cMass->Divide(nx, ny);
   
   TCanvas *cMasR = new TCanvas("cMasR", "Mass Ratios", dx, dy);
   cMasR->Divide(nx, ny);
   
   TLegend *leg = new TLegend(0.1, 0.7, 0.6, 0.9);
   leg->SetFillStyle(0);
   
   TH1D **houtputPerPtBin = new TH1D*[nptbins];
   TH1D **hSystFromRatio  = new TH1D*[ninputs];
   Double_t mean[nptbins][ninputs];
   TH1D *hSystMean = new TH1D(Form("hSysMean%s",name.Data()), "Systematic on the mean", nptbins, ptlims);
   
   for(Int_t ifile = 0; ifile<ninputs; ifile++){
      
      //need to do it only once per file
      TString hnamesuff = "", hnamepref = "";
      Bool_t isdifferentnamingscheme = kFALSE;
      if(hnamebase[ifile].Contains("Iter")){
      	 isdifferentnamingscheme = kTRUE;
      	 Int_t nchar = hnamebase[ifile].Sizeof();
      	 hnamesuff = hnamebase[ifile](6,nchar);
      	 hnamepref = hnamebase[ifile](0, 6);
      	 Printf("Pref %s Suff %s -> %d", hnamepref.Data(), hnamesuff.Data(), isdifferentnamingscheme);
      }
      cMass->SetName(Form("%s%s%s", cMass->GetName(), legs[ifile].Data(), logscale ? "logy" : ""));
      cMasR->SetName(Form("%s%s", cMasR->GetName(), legs[ifile].Data()));
      //read input file
      TFile *fin = new TFile(files[ifile]);
      if(!fin->IsOpen()) {
      	 Printf("%s not found, skip", files[ifile].Data());
      	 continue;
      }
      Printf("Reading file %s", files[ifile].Data());
      
      //read the histograms from the files
      //for the moment one histogram type at a time
      for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	 histogram[ifile][ipt] = 0x0;
      	 houtputPerPtBin[ipt] = 0x0;
      	 
      	 TString hnamethisloop = "", hnamedoublethisloop = "";
      	 if(isdifferentnamingscheme) {
      	    Printf("isdifferentnamingscheme");
      	    hnamethisloop = Form("%s%d%s", hnamepref.Data(),  ipt + offset[ifile], hnamesuff.Data());
      	    
      	     // for the comparison of the unfolded new
      	    
      	 }
      	 else {
      	 	 hnamethisloop = Form("%s%d", hnamebase[ifile].Data() ,  ipt + offset[ifile]);
      	 	 hnamedoublethisloop = Form("%s%.0f%.0fb", hnamebase[ifile].Data() ,  ptlims[ipt + offset[ifile]],  ptlims[ipt + 1 + offset[ifile]]);
      	 }
      	 histogram[ifile][ipt] = (TH1*)fin->Get(hnamethisloop);
      	 if(!histogram[ifile][ipt]) {
      	 	 Printf("Histogram %s not found, try %s ", hnamethisloop.Data(), hnamedoublethisloop.Data());
      	 	 histogram[ifile][ipt] = (TH1*)fin->Get(hnamedoublethisloop);
      	 	 if(!histogram[ifile][ipt]) {
      	 	 	 Printf("Histogram %s and %s not found ", hnamethisloop.Data(), hnamedoublethisloop.Data());
      	 	 	 fin->ls();
      	 	 	 continue;
      	 	 }
      	 }
      	 histogram[ifile][ipt]->SetName(Form("%s_%s_id%d", histogram[ifile][ipt]->GetName(), name.Data(), ifile));
      	 histogram[ifile][ipt]->SetMarkerSize(1);
      	 Printf("Read hisrogram %s", histogram[ifile][ipt]->GetName());
      }
   }

   TFile *fout = 0x0;
   
   if(writeRatios) {
      fout = new TFile(Form("Ratio%sOver%s.root", legs[0].Data(), legs[1].Data()), "recreate");
   }
   Int_t color = 0;
   for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      Bool_t firsttimeBase = kTRUE;
      Bool_t firsttime = kTRUE;
      Int_t ihindex = -1;
      for(Int_t ih = 0; ih < ninputs; ih++){
      	 if(!histogram[ih][ipt]) {
      	 	 Printf("histogram %d %d not read", ih, ipt);
      	 	 continue;
      	 }
      	 Printf("Debug ih = %d", ih);
      	 
      	 TH1* h1 = 0x0;
      	 TH1* hbase = 0x0;
      	 TH1* hR = 0x0;
      	 
      	 if(ih == base) continue;
      	 ihindex++;
      	 if(changeColor) {
      	    
      	    histogram[ih][ipt]->SetLineColor(colors[color+ih]);
      	    histogram[ih][ipt]->SetMarkerColor(colors[color+ih]);
      	    
      	    histogram[base][ipt]->SetLineColor(colors[color+base]);
      	    histogram[base][ipt]->SetMarkerColor(colors[color+base]);
      	 }
      	 mean[ipt][ih] = histogram[ih][ipt]->GetMean();
      	 mean[ipt][base] = histogram[base][ipt]->GetMean();
      	 if(!noUniform) UniformTH1FForDivide(histogram[ih][ipt], histogram[base][ipt], h1, hbase, "TH1D", noUniform);
      	 else {
      	    Printf("The histogram were not manipulated");
      	    h1 = histogram[ih][ipt];
      	    hbase = histogram[base][ipt];
      	    
      	 }
      	 
      	 if(ipt == 1) {
      	    leg->AddEntry(histogram[ih][ipt], legs[ih], "LP");
      	    if(firsttimeBase) leg->AddEntry(histogram[base][ipt], legs[base], "LP");
      	 }
      	 cMass->cd(ipt+1);
      	 if(logscale) gPad->SetLogy();
      	 if(h1) {
      	    if(!noUniform) h1->Scale(1./h1->Integral());
      	    h1->GetXaxis()->SetRangeUser(0, maxRangeMass[ipt]);
      	    if(logscale) h1->GetYaxis()->SetRangeUser(1.e-3, 1.5*h1->GetMaximum());
      	    else h1->GetYaxis()->SetRangeUser(0., 1.5*h1->GetMaximum());
      	    if(firsttime)  {
      	       h1->Draw();
      	       firsttime = kFALSE;
      	    }
      	    else h1->Draw("sames");
      	    hR = (TH1*)h1->Clone(Form("hR_%d_pt%d", ih, ipt));
      	    hR->GetYaxis()->SetRangeUser(-0., 2.);
      	    hR->GetYaxis()->SetTitle(Form("%s/%s", legs[ih].Data(), legs[base].Data()));
      	 }
      	 if(hbase){
      	    if(!noUniform) hbase->Scale(1./hbase->Integral());
      	    hbase->GetXaxis()->SetRangeUser(0, maxRangeMass[ipt]);
      	    if(firsttimeBase) hbase->Draw("sames");
      	    if(hR && hbase) hR->Divide(hbase);
      	    if(writeRatios) hbase->Write();
      	    firsttimeBase = kFALSE;
      	 }
      	 
      	 
      	 if(!noUniform){
      	    cMasR->cd(ipt+1);
      	    if(ninputs <= 2) hR->Draw();
      	    else hR->Draw("sames");
      	 }
      	 if(writeRatios) {
      	    h1->Write();
      	    if(!noUniform) hR->Write();
      	 }
      	 
      	 hSystFromRatio[ihindex] = 0x0;
      	 hSystFromRatio[ihindex] = SystematicFromRatio((TH1D*)hR, Form("hSyst%d%s_Pt%.0f_%.0f",ih, name.Data(), ptlims[ipt], ptlims[ipt+1]));
      	 //Printf("Index %d, Sys from ratio pointer %p (%s)",ihindex, hSystFromRatio[ihindex], hSystFromRatio[ihindex]->GetName());
      }
      if(ninputs > 2) {
      	  Printf("Systematics as standard deviation of the ratios");
      	  houtputPerPtBin[ipt] = MaxSpread(hSystFromRatio, ninputs -1, ipt, Form("hSystMx%s_Pt%.0f_%.0f",name.Data(), ptlims[ipt], ptlims[ipt+1]));
      
      	  cMasR->cd(ipt+1);
      	  
      	  houtputPerPtBin[ipt]->Draw("E2sames");
      }
      else {
      	      	  
      	  houtputPerPtBin[ipt] = hSystFromRatio[0];
      	  //Printf("Pointer is %p (%s)", houtputPerPtBin[ipt], houtputPerPtBin[ipt]->GetName());
      	  cMasR->cd(ipt+1);
      	  houtputPerPtBin[ipt]->Draw("E2sames");
      }
      
      cMass->cd(ipt+1);
      leg->Draw();
      
      // calculate the systematics on the mean
      Double_t erryu, erryd;
      MaxSpread(ninputs, mean[ipt], base, erryu, erryd);
      Double_t err = (erryu + TMath::Abs(erryd))/2.;
      hSystMean->SetBinContent(ipt+1, mean[ipt][base]);
      hSystMean->SetBinError(ipt+1, err);
   }
   if(writeRatios) {
      fout->cd();
      hSystMean->Write();
      Printf("The output was saved in a ROOT file");
      //fout->Close();
   }
  
   SaveCv(cMass);
   SaveCv(cMasR);
   
   return houtputPerPtBin;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TH1D* SystematicFromRatio(TH1D* hRatio, TString systname){
   
   TH1D* hsyst = (TH1D*)hRatio->Clone(systname);
   hsyst->Reset();
   hsyst->SetLineWidth(2);
   hsyst->SetFillStyle(0);
   hsyst->GetYaxis()->SetTitle("Syst = Ratio - 1");
   Printf("%s", systname.Data());
   for(Int_t ib = 0 ; ib < hsyst->GetNbinsX(); ib++){
      hsyst->SetBinContent(ib+1, 1);
      Double_t ratioVal = hRatio->GetBinContent(ib+1), error = TMath::Abs(hRatio->GetBinContent(ib+1) - 1);
      
      if(ratioVal < 1e-5) error = 0;
      hsyst->SetBinError(ib+1, error);
   }
   return hsyst;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TH1D* MaxSpread(TH1D** hsystratio, Int_t nh, Int_t index, TString name){
	
	TH1D* hmaxspread = (TH1D*)hsystratio[0]->Clone(Form("%s_Pt%.0f_%.0f", name.Data(), ptlims[index], ptlims[index+1]));
	hmaxspread->SetLineColor(1);
	hmaxspread->Reset();
   
   
	for(Int_t ib = 0 ; ib < hmaxspread->GetNbinsX(); ib++){
		Double_t min = 1, max = 1;
		//Printf("------------Start Bin %d", ib);
		Double_t content[nh];
		Double_t average = 0;
		for(Int_t ihs = 0; ihs<nh; ihs++){
			if(!hsystratio[ihs]) continue;
			//Printf("%p", hsystratio[ihs]);
			Double_t ratiobin = hsystratio[ihs]->GetBinError(ib+1);
			//Printf("Content = %f, error %f", ratiobin, hsystratio[ihs]->GetBinError(ib+1));
			content[ihs] = ratiobin;
			average += ratiobin;
			if(ratiobin < min) min = ratiobin ;
			if(ratiobin > max) max = ratiobin ;
      }
      Double_t dnh = (Double_t)nh;
      average /= dnh;
      Double_t devst = 0;
      for(Int_t ihs = 0; ihs<nh; ihs++) devst += (content[ihs]-average)*(content[ihs]-average);
      devst /= dnh; devst = TMath::Sqrt(devst);
      hmaxspread->SetBinContent(ib+1, 1);
      
      hmaxspread->SetBinError(ib+1, average);
      //Printf("Total = %f", average);
      //Printf("------------End Bin %d", ib);
   }
   
   return hmaxspread;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void MaxSpread(const Int_t nval, Double_t values[], Int_t idbase, Double_t& up, Double_t& down){
	
	up = 0; down = 0;
	for(Int_t i = 0; i< nval; i++){
		if(i == idbase) continue;
		Double_t signeddiff = values[i]-values[idbase];
		if(signeddiff > 0 && signeddiff>up) up = signeddiff;
		if(signeddiff < 0 && signeddiff<down) down = signeddiff;
	}

	Printf("%f + %f - %f", values[idbase], up, down);
	
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TH1D* AddInQuadrature(TH1D** haddq, Int_t nh, Int_t index, TString name){
   
	TH1D* hsyst = 0x0;
	for(Int_t ihs = 0; ihs<nh; ihs++){
		if(!haddq[ihs]) continue;
		if(!hsyst) {
			TString namesystot = name;
			if(index > -1) namesystot = Form("%s_Pt%.0f_%.0f", name.Data(), ptlims[index], ptlims[index+1]);
			hsyst = (TH1D*)haddq[ihs]->Clone(namesystot);
			hsyst->Reset();
			break;
		}
	}
	for(Int_t ib = 0 ; ib < hsyst->GetNbinsX(); ib++){
		Double_t addq = 0;
		//Printf("------------Start Bin %d", ib);
		for(Int_t ihs = 0; ihs<nh; ihs++){
			if(!haddq[ihs]) continue;
			
			Double_t errbin = haddq[ihs]->GetBinError(ib+1);
			addq += (errbin * errbin) ;
			//Printf("%s Add %f  = %f ", haddq[ihs]->GetName(), errbin*errbin, addq);
		}
		
		hsyst->SetBinContent(ib+1, 1);
		
		hsyst->SetBinError(ib+1, TMath::Sqrt(addq));
		//Printf("Total = %f", TMath::Sqrt(addq));
		//Printf("------------End Bin %d", ib);
	}
	
	hsyst->SetLineColor(1);
	
	return hsyst;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TH1D* SetMassValueInSystematic(TH1D* hSyst, TH1D* hMass){
   
   /// read the mass spectrum and the absolute systematic uncertainty. Modify the bin content of the latter with the value from the former and return an histogram with the relative systematic uncertainties and value 1 (or 0, tbd)
   
   if(!hSyst || !hMass){
      Printf("Inalid input (syst = %p, mass = %p)", hSyst, hMass);
      return 0x0;
   }
   
   TString hrelsysname = hSyst->GetName();
   hrelsysname.Insert(1, "Rel");
   TH1 *hh;
   TH1 *gg;
   UniformTH1FForDivide(hSyst, hMass, hh, gg, "TH1D", kTRUE);
   delete hh;
   delete gg;
   hMass->SetMarkerSize(1);
   TH1D* hRelSystUnc = (TH1D*) hSyst->Clone(hrelsysname);
   for(Int_t ib = 0; ib < hRelSystUnc->GetNbinsX(); ib++){
      Double_t absEr = hSyst->GetBinError(ib+1);
      Double_t mass = hMass->GetBinContent(ib+1);
      Double_t relEr = 0;
      if(TMath::Abs(mass) > 1e-5) relEr = absEr/mass;
      
      hRelSystUnc->SetBinError(ib+1, relEr);
      hRelSystUnc->SetBinContent(ib+1, 0);
      
      hSyst->SetBinContent(ib+1, mass);
      hSyst->SetBinError(ib+1, mass*absEr);
      
      //For debug
      //Printf("M  = %.3f +- %.2f ", mass, mass*absEr);
   }
   
   return hRelSystUnc;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TList* AddSystematicstoMassFromFile(TString inMfilename, TString hbasenamemass, Bool_t usebinlimsForM,  TString inSysfilename, TString hbasenamesys, Bool_t usebinlimsForSys){
   
   /// read input file inMfilename and histogram hbasenamemass with absolute systematic uncertainty
   /// read total systematic uncertainties in file inSysfilename with name hbasenamesys
   /// returns  a TList with the mass histogram with systematic uncertainties as error, relative systematic uncertainties, mass spectrum in the input file with stat unc only
   
  TList *listM = ReadHistogramsInFile(inMfilename, hbasenamemass, usebinlimsForM, nptbins);
  if(!listM) return 0x0;
  TList *listSys = ReadHistogramsInFile(inSysfilename, hbasenamesys, usebinlimsForSys, nptbins);
  if(!listSys) return 0x0;

  //
  TH1** hMStat = new TH1*[nptbins];
  TH1D** hMSyst = new TH1D*[nptbins];
  TH1** hAbsSyst = new TH1*[nptbins];
  TH1D** hRelSyst = new TH1D*[nptbins];
  TH1D** hRelStat = new TH1D*[nptbins];
  
  TList *listResults = new TList();
  listResults->SetOwner();
  
  for(Int_t i = 0; i<nptbins; i++){
     
     //hMStat[i]    = (TH1D*) listM->At(i);
     //hAbsSyst[i] = (TH1D*) listSys->At(i);
     
     TH1D *hMStattmp   = (TH1D*) listM->At(i);
     TH1D* hAbsSysttmp = (TH1D*) listSys->At(i);
     Printf("Mass before adding syst %s, mean value %f", hMStattmp->GetName(), hMStattmp->GetMean());
     UniformTH1FForDivide(hMStattmp, hAbsSysttmp, hMStat[i], hAbsSyst[i]);
     hRelSyst[i] = SetMassValueInSystematic((TH1D*)hAbsSyst[i], (TH1D*)hMStat[i]);
     hRelStat[i] = SetMassValueInSystematic((TH1D*)hMStat[i]->Clone(Form("hMassCopy%d", i)), (TH1D*)hMStat[i]);
     
     Printf("Mass AFTER adding syst %s, mean value %f", hMStat[i]->GetName(), hMStat[i]->GetMean());
     //number of plots added to the list = 4
     listResults->Add(hMStat[i]);
     listResults->Add(hAbsSyst[i]); // this is the systematic uncertainty at the value of the mass
     listResults->Add(hRelSyst[i]);
     listResults->Add(hRelStat[i]);
  }

   return listResults;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TList* ReadHistogramsInFile(TString infilename, TString basehname, Bool_t usebinlims, Int_t n){
   /// open the file, read the n TH1D with name basehname%d, put them in a TList and return it
   /// this method could go to the CommonTools.....
   
   TFile *fin = new TFile(infilename);
   if(!fin->IsOpen()){
      Printf("File %s not found", infilename.Data());
      return 0x0;
   }

   //need to do it only once per file
   TString hnamesuff = "", hnamepref = "";
   Bool_t isdifferentnamingscheme = kFALSE;
   if(basehname.Contains("Iter")){
   	   isdifferentnamingscheme = kTRUE;
   	   hnamesuff = basehname(6,11);
   	   hnamepref = basehname(0, 6);
   	   Printf("Pref %s Suff %s -> %d", hnamepref.Data(), hnamesuff.Data(), isdifferentnamingscheme);
   }

   TList *listout = new TList();
   listout->SetOwner(kTRUE);
   
   for(Int_t i = 0; i<n; i++){
   	   
   	   TString name = "";
   	   if(isdifferentnamingscheme) {
   	   	   Printf("isdifferentnamingscheme");
   	   	   name = Form("%s%d%s", hnamepref.Data(),  i, hnamesuff.Data());
   	   	   
   	   	   // for the comparison of the unfolded new
   	   	   
   	   }
   	   if(usebinlims) name = Form("%s%.0f_%.0f", basehname.Data(), ptlims[i], ptlims[i+1]);
   	   else  Form("%s%d", basehname.Data(), i);
   	   TH1D *h = (TH1D*)fin->Get(name);
   	   if(!h){
      	 Printf("%s not found", name.Data());
      	 fin->ls();
      	 continue;
      }
      
      listout->Add(h);
   
   }
   
   return listout;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TH1D* SmoothUncertaintyByAverage(const Int_t nh, TH1D** hinputUnc, const char* newname){
	
	//check that the input histograms have the same bins and range
	Int_t nxbins[nh];
	Int_t nsysbins = 0;
	TH1D *hsmoothsys = 0x0;
	
	for(Int_t ih = 0; ih < nh; ih++){
		if(!hinputUnc[ih]) {
			Printf("Histogram number %d not found", ih);
			continue;
		}
		nxbins[ih] = hinputUnc[ih]->GetNbinsX();
		if(ih == 0) nsysbins = nxbins[ih];
		else if(nsysbins > nxbins[ih]) {
			Printf("Different number of bins Current = %d vs %d", nxbins[ih], nsysbins);
			nsysbins = nxbins[ih];
			
		}
		
		if(!hsmoothsys) hsmoothsys = (TH1D*)hinputUnc[ih]->Clone(newname);
	}
	
	Double_t avgSys[nsysbins];
		
	for(Int_t isyb = 0; isyb < nsysbins; isyb++){
		avgSys[isyb] = 0;
		Int_t countfail = 0;
		for(Int_t ih = 0; ih < nh; ih++){
			if(!hinputUnc[ih]) {
				Printf("Histogram number %d not found", ih);
				countfail++;
				continue;
			}
			avgSys[isyb] += hinputUnc[ih] ->GetBinError(isyb+1);
			//Printf("isyb = %d, ih = %d -> %f", isyb, ih, avgSys[isyb]);
		}
		Double_t den = (Double_t)(nh-countfail);
		avgSys[isyb] /= den;
		hsmoothsys->SetBinError(isyb+1, avgSys[isyb]);
	}
	
	// uncomment for debugging
	//TCanvas *c = new TCanvas(Form("c%s", hsmoothsys->GetName()), Form("c%s", hsmoothsys->GetName()));
	//c->cd();
	//hsmoothsys->Draw();
	
	return hsmoothsys;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

TH1D** SmoothUnceraintyByFit(const Int_t nh, TH1D** hinputUnc, Int_t fittype, const char* newname, TH1D** hSysOnError, Double_t xmin, Double_t xmax){
	
	if(!hinputUnc[0]) return 0x0;
	if(xmin == -1) xmin = hinputUnc[0]->GetXaxis()->GetBinLowEdge(1);
	if(xmax == -1) xmax = hinputUnc[0]->GetXaxis()->GetBinLowEdge(hinputUnc[0]->GetNbinsX()+1);
	
	TH1D** hnewUnc = new TH1D*[nh];
	
	Int_t nx, ny, dx, dy;
	CalculatePads(nh, nx, ny, dx, dy);
	
	TCanvas *cSmooth = new TCanvas(Form("cSmooth%s", newname), Form("cSmooth%s", newname), dx, dy);
	
	cSmooth->Divide(nx, ny);
	
	for(Int_t ih = 0; ih<nh; ih++){
		
		if(!hinputUnc[ih]) {
			Printf("Histo %d not found", ih);
			continue;
		}
		hnewUnc[ih] = (TH1D*)hinputUnc[ih]->Clone(Form("%s%d", newname, ih));
		hnewUnc[ih]->Reset();
		hSysOnError[ih] = (TH1D*)hinputUnc[ih]->Clone(Form("%sSysOnErr%d", newname, ih));
		hSysOnError[ih]->Reset();
		
		for(Int_t ib = 0; ib<hnewUnc[ih]->GetNbinsX(); ib++){
			hnewUnc[ih]->SetBinContent(ib+1, hinputUnc[ih]->GetBinError(ib+1));
			
		}
		TF1 *fpol1 = new TF1(Form("fpol1%s_id%d", newname, ih), "[0]+[1]*x", xmin, xmax);
		fpol1->SetLineWidth(2);
		fpol1->SetLineStyle(2);
		fpol1->SetLineColor(kBlue);
		fpol1->SetParameters(5., 0.1);
		
		TF1 *fpol2 = new TF1(Form("fpol2%s_id%d", newname, ih), "[0]+[1]*x + [2]*x*x", xmin, xmax);
		fpol2->SetLineWidth(2);
		fpol2->SetLineStyle(3);
		fpol2->SetLineColor(kMagenta);
		fpol2->SetParameters(5., 0.1, 0.1);
		
		Int_t fitres1 = hnewUnc[ih]->Fit(fpol1, "RL0+S", "", xmin, xmax);
		Int_t fitres2 = hnewUnc[ih]->Fit(fpol2, "RL0+S", "", xmin, xmax);
		//Printf("Fit result 1 = %d, fit result 2  = %d", fit1res, fit2res);
		
		cSmooth->cd(ih+1);
		hnewUnc[ih]->GetYaxis()->SetRangeUser(0., 1);
		hnewUnc[ih]->DrawClone();
		fpol1->Draw("sames");
		fpol2->Draw("sames");
		
		Int_t nbinsx = hnewUnc[ih]->GetNbinsX();
		//Bool_t exclufromfit[nbinsx];
		//TH1D *hexclu = new TH1D(Form("hexclu%s%d", newname, ih), Form("hexclu%s%d", newname, ih), nbinsx, xmin, xmax);
		
		for(Int_t ib = 0; ib<nbinsx; ib++){
			//if(fitres2 == 0){
				if (fittype == 2){
								
				hnewUnc[ih]->SetBinContent(ib+1, fpol2->Eval(hnewUnc[ih]->GetBinCenter(ib+1)));
				
				hSysOnError[ih]->SetBinContent(ib+1, 1); 
				hSysOnError[ih]->SetBinError(ib+1, fpol2->Eval(hnewUnc[ih]->GetBinCenter(ib+1)));
			} 
			//else {
			if (fittype == 1){
				hnewUnc[ih]->SetBinContent(ib+1, fpol1->Eval(hnewUnc[ih]->GetBinCenter(ib+1)));
				
				hSysOnError[ih]->SetBinContent(ib+1, 1); 
				hSysOnError[ih]->SetBinError(ib+1, fpol1->Eval(hnewUnc[ih]->GetBinCenter(ib+1)));
			}
			
			//Double_t reldiff = TMath::Abs((hnewUnc[ih]->GetBinContent(ib+1) - hinputUnc[ih]->GetBinContent(ib+1))/hinputUnc[ih]->GetBinContent(ib+1));
			//
			//if(reldiff > 1.) {
			//	Printf("Rel diff = %f", reldiff);
			//	exclufromfit[ib] = kTRUE;
			//	hexclu->SetBinContent(ib+1, hnewUnc[ih]->GetBinContent(ib+1));	
			//}
			//else exclufromfit[ib] = kFALSE;
		}
		
		
		
		
		cSmooth->cd(ih+1);
		hnewUnc[ih]->Draw("sames");
		//hexclu->Draw("sames");
	}
	
	SaveCv(cSmooth);
	return hnewUnc;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void CompareThreeResultFixedInput(){
   
   const Int_t ninputs = 3;
   
   TString files[ninputs] = 
   {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgDeriv.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/OneRhoBinNoWeight/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/BinByBinCorrDetPtParPt/BinbyBinCorrection_111BinT0.root"};
   // compare standard unfolding deriv with one rho bin Deriv and pythia

   TString hnamebase[ninputs] =
   {"hMUnf__Iter3", "hMUnf__Iter3", "hMPar_pt"};
   
   TString legs[ninputs] = 
   { "UnfDetFlDeriv", "UnfDetFlDerivRho15_30", "PYTHIAPar"};
   
   Int_t firstMatchingBin[ninputs] = 
   {0, 0, 0};
   
   Bool_t changeColor = kFALSE;
   Bool_t writeout = kTRUE;
   Bool_t nouniform = kFALSE;
   CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform);
}

void CompareBkgNoBkgResultsFixedInput(){
	const Int_t ninputs = 6;
   
   TString files[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo/00resp_pt20_ptT10/UnfoldedDistributionsPrior0DetFlBkgConst.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MB/ConstSub/UnfoldedDistributionsPrior0.root",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetPlusFluc/Embedding20160330pth2_PPtC10/BinbyBinCorrectedMass_3Constrb2.root", "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/BinByBinCorrDetPtParPt/BinbyBinCorrection_111BinT0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetFlucNoBkgSub/PtDetPtPar/BinbyBinCorrectedMass_3NoSub.root"};
   
   TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3", "hMassCorr_pt", "hMPar_pt", "hMassCorr_pt"};
   TString legs[ninputs] = {"UnfEmbNoBkgSubMB", "UnfEmbConstSubMB",  "NoFlucDataConst", "BbBEmbConstSub","PYTHIA", "BbBEmbNoBkgSub"};
   
   Int_t firstMatchingBin[ninputs] = {0, 0, 0, 0, 0, 0};
   
   Bool_t changeColor = kFALSE;
   Bool_t writeout = kTRUE;
   Bool_t nouniform = kFALSE;
   Int_t base = 4;
   Printf("Base %d is %s", base , legs[base].Data());
   CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);

   
   //plot only cmp ConstSub All + Pythia
   const Int_t ninputsCS = 3;
   TString filesCS[ninputsCS] = {files[1], files[3], files[4]};
   TString hnamebaseCS[ninputsCS] = {hnamebase[1], hnamebase[3], hnamebase[4]};
   TString legsCS[ninputsCS] = {legs[1], legs[3], legs[4]};
   Int_t firstMatchingBinCS[ninputsCS] = {firstMatchingBin[1], firstMatchingBin[3], firstMatchingBin[4]};
   base = ninputsCS -1;
   CompareResults(ninputsCS, filesCS, hnamebaseCS, legsCS, firstMatchingBinCS, changeColor, writeout, nouniform, base);
   
   //plot only cmp NoBkg + Const Sub unf + Pythia
   const Int_t ninputsNB = 4;
   TString filesNB[ninputsNB] = {files[1], files[0], files[5], files[4]};
   TString hnamebaseNB[ninputsNB] = {hnamebase[1], hnamebase[0], hnamebase[5], hnamebase[4]};
   TString legsNB[ninputsNB] = {legs[1], legs[0], legs[5], legs[4]};
   Int_t firstMatchingBinNB[ninputsNB] = {firstMatchingBin[1], firstMatchingBin[0], firstMatchingBin[5], firstMatchingBin[4]};
   base = ninputsNB -1;
   CompareResults(ninputsNB, filesNB, hnamebaseNB, legsNB, firstMatchingBinNB, changeColor, writeout, nouniform, base);
   
   //plot only cmp NoBkg + OldNo Fl + Pythia
   const Int_t ninputsFl = 3;
   TString filesFl[ninputsFl] = {files[2], files[0], files[4]};
   TString hnamebaseFl[ninputsFl] = {hnamebase[2], hnamebase[0], hnamebase[4]};
   TString legsFl[ninputsFl] = {legs[2], legs[0], legs[4]};
   Int_t firstMatchingBinFl[ninputsFl] = {firstMatchingBin[2], firstMatchingBin[0], firstMatchingBin[4]};
   base = ninputsFl -1;
   CompareResults(ninputsFl, filesFl, hnamebaseFl, legsFl, firstMatchingBinFl, changeColor, writeout, nouniform, base);
}

void CompareOvlNoOvl_NoBkgSub(){
	const Int_t ninputs = 2;
	TString files[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/UnfoldedDistributionsPrior0.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/NoOvlEx0/analysis/NoSub/unfoldFirstLook/UnfoldedDistributionsPrior0.root"};
	TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMUnf__Iter3"};
	TString legs[ninputs] = {"UnfEmbNoBkgSubMB", "UnfEmbNoBkgSubOvlExcluMB"};
	Int_t firstMatchingBin[ninputs] = {0, 0};
	Bool_t changeColor = kFALSE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);
}

void CompareFinalSysNoSubConstSub(){
	const Int_t ninputs = 2;
	TString files[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MBEJE/00resp_pt20_80or70_110_ptT10_140or50_140_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/UnfoldedDistributionsPrior0.root"};
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MBEJE/redo2/00resp_pt20_80or70_110_ptT10_140or50_140_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root"};
	TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMUnf__Iter3"};
	TString legs[ninputs] = {"UnfEmbNoBkgSub", "CompNoSub"};//"UnfEmbConstSub"};
	Int_t firstMatchingBin[ninputs] = {0, 0};
	Bool_t changeColor = kFALSE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);

}

void CompareNoSubUnfAndBbB(){
	const Int_t ninputs = 2;
	TString files[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MBEJE/00resp_pt20_80or70_110_ptT10_140or50_140_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/Bin-By-BinCorr/DetFlucNoBkgSub/PtDetPtPar/BinbyBinCorrectedMass_3NoSub.root"};
	TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMassCorr_pt"};
	TString legs[ninputs] = {"UnfEmbNoBkgSub", "BbBEmbNoBkgSub"};
	Int_t firstMatchingBin[ninputs] = {0, 0};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);

}

void CompareMissNoMissMB(){
const Int_t ninputs = 3;
	TString files[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/UnfoldedDistributionsPrior0Miss2D.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/UnfoldedDistributionsPrior0noMiss.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/CorrectedUnfoldedMassnoMissKine.root"
		//, "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/BinByBinCorrDetPtParPt/BinbyBinCorrection_111BinT0.root"
	};
	TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMUnf__Iter3", "hMassUnfCorrPt"
	//, "hMPar_pt"
	};
	TString legs[ninputs] = {"UnfEmbNoBkgSub", "UnfEmbNoBkgSubNoMiss", "UnfEmbNoBkgSubNoMissKineEff"
	//, "PYTHIAMatched"
	};
	Int_t firstMatchingBin[ninputs] = {0, 0, 0};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);
}

void ComparisonUnfNoBkgSubVsNoFluc(){
	const Int_t ninputs = 4;
	TString files[ninputs] = {
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/UnfoldedDistributionsPrior0DetFlBkgNoSub.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/UnfoldedDistributionsPrior0Miss2D.root",
		"/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/analysis/BinByBinCorrDetPtParPt/BinbyBinCorrection_111BinT0.root"
	};
	TString hnamebase[ninputs] = {
		"hMUnf__Iter3", "hUnfMass_PtBin", "hMUnf__Iter3", "hMPar_pt"
	};
	TString legs[ninputs] = {"UnfEmbNoBkgSubDef", "pPbMarta", "UnfEmbNoBkgSubrunalone", "PYTHIAMatched"
	};
	
	Int_t firstMatchingBin[ninputs] = {0, 2, 1, 1};
	
	
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);
}

void ComparisonDebug(){
	const Int_t ninputs = 5;
	TString files[ninputs] = {
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/testUnfEJE/UnfoldedDistributionsPrior0NoBkgSub.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/testUnfEJE/UnfoldedDistributionsPrior0NoBkgSubtoday.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/testUnfEJE/UnfoldedDistributionsPrior01134Analarge.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreAppr/CmpEmbBkgSub/RatioNoSubOverDetOnly.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreAppr/CmpEmbBkgSub/RatioNoSubOverDetOnly.root"
	};
	
	TString hnamebase[ninputs] = {
		"hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3_BkgSub_id0b", "hMUnf__Iter3_BkgSub_id1b"
	};
	
	TString legs[ninputs] = {"UnfEmbNoBkgSubTodayAutResp","UnfEmbNoBkgSubToday", "DetOnlyToday", "UnfEmbNoBkgSubPreApp", "UnfEmbDetOnlyPreApp"
	};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kFALSE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Int_t firstMatchingBin[ninputs] = {0, 0, 0, 0, 0};
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);
	
}

void ComparisonDebugSyst(){
	const Int_t ninputs = 5;
	TString files[ninputs] = {
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160905/CmpEmbBkgSub/RatioNoSubOverDetOnly.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160905/CmpEmbBkgSub/RatioNoSubOverDetOnly.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/testUnfEJE/UnfoldedDistributionsPrior01134Analarge.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreAppr/CmpEmbBkgSub/RatioNoSubOverDetOnly.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreAppr/CmpEmbBkgSub/RatioNoSubOverDetOnly.root"
	};
	
	TString hnamebase[ninputs] = {
		"hMUnf__Iter3_BkgSub_id0b", "hMUnf__Iter3_BkgSub_id1b", "hMUnf__Iter3", "hMUnf__Iter3_BkgSub_id0b", "hMUnf__Iter3_BkgSub_id1b"
	};
	
	TString legs[ninputs] = {"UnfEmbNoBkgSubSysToday", "UnfEmbDetOnlySysToday", "DetOnlyToday", "UnfEmbNoBkgSubPreApp", "UnfEmbDetOnlyPreApp"
	};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kFALSE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Int_t firstMatchingBin[ninputs] = {0, 0, 0, 0, 0};
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);
	
}

void ComparisonDebugMBEJE(){
	const Int_t ninputs = 5;
	TString files[ninputs] = {
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MBEJE/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/ConstSub/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/testUnfEJE/UnfoldedDistributionsPrior01134Analarge.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreAppr/CmpEmbBkgSub/RatioNoSubOverDetOnly.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/SystPreAppr/CmpEmbBkgSub/RatioNoSubOverDetOnly.root"
	};
	
	TString hnamebase[ninputs] = {
		"hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3_BkgSub_id0b", "hMUnf__Iter3_BkgSub_id1b"
	};
	
	TString legs[ninputs] = {"UnfEmbNoBkgSubMBEJE","UnfEmbDetOnlyMBEJE", "DetOnlyToday", "UnfEmbNoBkgSubPreApp", "UnfEmbDetOnlyPreApp"
	};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kFALSE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Int_t firstMatchingBin[ninputs] = {0, 0, 0, 0, 0};
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);
	
}

void CompareFixedVariableBinW(Bool_t logy = kFALSE){
	const Int_t ninputs = 5;
	TString files[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/MB/UnfoldedDistributionsPrior0VarBinWidth.root",	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/EJE/UnfoldedDistributionsPrior0VarBinWidth.root",
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/UnfoldedDistributionsPrior0DetFlBkgNoSub.root",
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/EJE/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/UnfoldedDistributionsPrior0DetFlBkgNoSub.root",
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MBEJE/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root"
	};
	//these below are the test that match best
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/EJE/tests/UnfoldedDistributionsPrior010GeVpar.root",
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/EJE/tests/UnfoldedDistributionsPrior0.root","/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/EJE/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/UnfoldedDistributionsPrior0DetFlBkgNoSub.root"

	TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3"};
	TString legs[ninputs] = {"MBVar","EJEVar", "MBFix","EJEFix", "MBEJEFix"};
	Int_t firstMatchingBin[ninputs] = {0, 0, 0, 0, 0};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 2;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base, "", logy);

}

void CompareFixedVariableBinWCentralValueMBEJE(Bool_t logy = kFALSE){
	const Int_t ninputs = 2;
	TString files[ninputs] = {
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MBEJE/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/MassUnfSumMBEJE.root",
		"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/MBEJE/MassUnfSumMBEJE.root"	

	};
	//these below are the test that match best
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/EJE/tests/UnfoldedDistributionsPrior010GeVpar.root",
	//"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/EJE/tests/UnfoldedDistributionsPrior0.root","/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/EJE/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/UnfoldedDistributionsPrior0DetFlBkgNoSub.root"

	TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMUnf__Iter3"};
	TString legs[ninputs] = {"MBEJEFix", "MBEJEVar"};
	Int_t firstMatchingBin[ninputs] = {0, 0};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base, "", logy);

}
//______________________________________________________________________________

void ComparePriorVar(){
	const Int_t ninputs = 4;
	TString files[ninputs] = {
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB/00resp_pt20_80or70_120_ptT10_150or50_150_m0_12or0_14_mT0_40/UnfoldedDistributionsPrior0DetFlBkgNoSub.root",
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/analysis/NoSub/testPriorVariation/UnfoldedDistributionsPrior12percv2.root",
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/analysis/NoSub/testPriorVariation/UnfoldedDistributionsPrior15percv2.root",
	"/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160726/Const/analysis/NoSub/testPriorVariation/UnfoldedDistributionsPrior110percv2.root"
	};
	TString hnamebase[ninputs] = {"hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3", "hMUnf__Iter3"};
	TString legs[ninputs] = {"Std", "Varied2perc", "Varied5perc", "Varied10perc"};
	Int_t firstMatchingBin[ninputs] = {0, 0, 0, 0};
	Bool_t changeColor = kTRUE;
	Bool_t writeout = kTRUE;
	Bool_t nouniform = kFALSE;
	Int_t base = 0;
	Printf("Base %d is %s", base , legs[base].Data());
	CompareResults(ninputs, files, hnamebase, legs, firstMatchingBin, changeColor, writeout, nouniform, base);

}

//______________________________________________________________________________

TH1D** GetMeanSyst(const Int_t ninputs, TString filenames[], TString hnameinput[]){
	// TString filenames[ninputs] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Syst20160915/CmpEmbBkgSub/RatioNoSubOverDetOnly.root", ..}
	// TString hnameinput[ninputs] = {"hSysMeanBkgSub", ..}
	
	TH1D** hsysPart = new TH1D*[ninputs];
	for(Int_t isys = 0; isys<ninputs; isys++){
		hsysPart[isys] = 0x0;
		TFile *fin = new TFile(filenames[isys]);
		if(!fin->IsOpen()){
			Printf("Error, file %s not found", filenames[isys].Data());
			
			continue;
		}
		hsysPart[isys] = (TH1D*)fin->Get(hnameinput[isys]);
		
	}
	return hsysPart;
}
//______________________________________________________________________________

TH1D* SumMeanSyst(const Int_t ninputs, TH1D** hsysPart){
	
	
	TH1D *hMeanSystTot = AddInQuadrature(hsysPart, ninputs, -1, "hMeanSystTot");
	
	return hMeanSystTot;
	
}


#endif

