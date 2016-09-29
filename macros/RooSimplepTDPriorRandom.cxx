//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "TFile.h"
#include "TVectorD.h"

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TProfile.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldTestHarness2D.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================





TH2D* CorrelationHist (const TMatrixD& cov,const char* name, const char* title,
	Double_t lo, Double_t hi,Double_t lon,Double_t hin)
{
	Int_t nb= cov.GetNrows();
	Int_t na= cov.GetNcols();
	cout<<nb<<" "<<na<<endl;
	TH2D* h= new TH2D (name, title, nb, 0, nb, na, 0, na);
	h->SetAxisRange (-1.0, 1.0, "Z");
	for(int i=0; i < na; i++)
	for(int j=0; j < nb; j++) {
		Double_t Viijj= cov(i,i)*cov(j,j);
		if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
    }
    return h;
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooSimplepTDPriorRandom(Int_t rad,Int_t eff,Int_t frac,Int_t gaus)
{
	#ifdef __CINT__
	gSystem->Load("libRooUnfold");
	#endif
	Int_t difference=1;
	Int_t Ppol=0;
	cout << "==================================== pick up the response matrix for background==========================" << endl;
	///////////////////parameter setting
	RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
	
	//Get the tree for MC
	if(eff==0) TString fname = Form("/Users/leticia/2016WorkDirQGTag/Responsepp/DetRespPythia.root",rad);
	if(eff==1) TString fname = Form("/Users/leticia/2016WorkDirQGTag/ResponseppEffLess/DetRespPythia.root",rad);
	TList *list1;
	TFile *input = TFile::Open( fname );
    TTree *mc=(TTree*)input->Get("newtree");
    
    //Get the tree for data
    TString fname2 = "/Users/leticia/2016WorkDirQGTag/AnalysisResults_Datapp_NoSub.root";
    
    TList *list2;
    TFile *input2 = TFile::Open( fname2 );
    list2=(TList*) input2->Get(Form("JetQGTaggings_Jet_AKTChargedR0%i0_tracks_pT0150_E_scheme_TCRaw_Data_NoSub_Incl",rad)); 
    TTree *data=(TTree*)list2->FindObject("fTreeJetShape");
    
    fRandom=new TRandom3(0);
    
    
    //DEFINE SOME HISTOGRAMS ************************************************************
    
    //the raw correlation
    TH2F *h2raw(0);
    h2raw=new TH2F("raw","raw",7,0.3,1,6,20,80);
    //detector measure level
    TH2F *h2smeared(0);
    h2smeared=new TH2F("smeared","smeared",7,0.3,1,6,20,80);
    TH2F *h2smeartest(0);
    h2smeartest=new TH2F("smearedtest","smearedtest",7,0.3,1,6,20,80); 
    //detector measured level but no cuts
    TH2F *h2smearednocuts(0);
    h2smearednocuts=new TH2F("smearednocuts","smearednocuts",10,0.,1,8,0,160);
    //true correlations with measured cuts
    TH2F *h2true(0);
    h2true=new TH2F("true","true",10,0,1.,8,0,160);
    
    TH2F *h2fulleff(0);
    h2fulleff=new TH2F("truef","truef",10,0.,1.,8,0,160);
    TH2F *h2fulleffnocuts(0);
    h2fulleffnocuts=new TH2F("fullnocuts","fullnocuts",10,0.,1,8,0,160);
    
    TH2F *h2fulleffFast(0);
    h2fulleffFast=new TH2F("truefFast","truefFast",10,0.,1.,8,0,160);
    
    TH2F *hcovariance(0);
    hcovariance=new TH2F("covariance","covariance",10,0.,1.,10,0,1.);
    
    TH2F *histomult[11][9];
    for(Int_t i=0;i<11;i++){
    	for(Int_t j=0;j<9;j++){
    histomult[i][j]=new TH2F(Form("aver%i%i",i,j),"aver",7,0.3,1,6,20,80);}}
    
    TH2F *effnum=(TH2F*)h2fulleff->Clone("effnum");
    TH2F *effdenom=(TH2F*)h2fulleff->Clone("effdenom");
    TH2F *h2fulleffsm=(TH2F*)h2fulleff->Clone("h2fulleffsm");
    TH2F *h2truesm=(TH2F*)h2true->Clone("h2truesm");
    effnum->Sumw2();
    effdenom->Sumw2();
    h2smeared->Sumw2();
    h2true->Sumw2();
    h2raw->Sumw2();
    h2fulleff->Sumw2();
    h2smearednocuts->Sumw2();
    h2fulleffnocuts->Sumw2();
    h2fulleff->Sumw2();
    h2fulleffsm->Sumw2();
    h2truesm->Sumw2();
    Float_t partonCode=0., ptJet=0., ptDJet=0., mJet=0., nbOfConst=0., angularity=0., circularity=0., lesub=0., sigma2=0.;
    Float_t ptJetMatch=0., ptDJetMatch=0., mJetMatch=0., nbOfConstMatch=0., angularityMatch=0., circularityMatch=0., lesubMatch=0., sigma2Match=0.;
    Double_t weightPythiaFromPtHard=0.;
    Float_t weightPythia=0; 
    Float_t partonCodem=0., ptJetm=0., ptDJetm=0., mJetm=0., nbOfConstm=0., angularitym=0., circularitym=0., lesubm=0., sigma2m=0.;
    Float_t ptJetMatchm=0., ptDJetMatchm=0., mJetMatchm=0., nbOfConstMatchm=0., angularityMatchm=0., circularityMatchm=0., lesubMatchm=0., sigma2Matchm=0.,weightPythiam=0.; 
    Float_t partonCode2=0., ptJet2=0., ptDJet2=0., mJet2=0., nbOfConst2=0., angularity2=0., circularity2=0., lesub2=0., sigma22=0.;
    Float_t ptJetMatch2=0., ptDJetMatch2=0., mJetMatch2=0., nbOfConstMatch2=0., angularityMatch2=0., circularityMatch2=0., lesubMatch2=0., sigma2Match2=0.,weightPythia2=0.; 
    
    
    
    Int_t nEv=0;; 
    
    
    //FILL THE RAW DATA CORRELATION ************************************   
    cout<<"cucu"<<endl;
    nEv=data->GetEntries(); 
    
    data->SetBranchAddress("partonCode", &partonCode); 
    data->SetBranchAddress("ptJet", &ptJet); 
    data->SetBranchAddress("ptDJet", &ptDJet); 
    data->SetBranchAddress("angularity", &angularity); 
    data->SetBranchAddress("lesub", &lesub); 
    data->SetBranchAddress("ptJetMatch", &ptJetMatch); 
    data->SetBranchAddress("ptDJetMatch", &ptDJetMatch); 
    data->SetBranchAddress("angularityMatch", &angularityMatch); 
    data->SetBranchAddress("lesubMatch", &lesubMatch); 
    data->SetBranchAddress("weightPythia", &weightPythia);
    double count=0;
    for(int iEntry=0; iEntry< nEv; iEntry++){
    	data->GetEntry(iEntry); 
    	if(ptJet>80 || ptJet<20) continue;
    	if(ptDJet<0.3) continue;
    	
    h2raw->Fill(ptDJet,ptJet);}
    
    //rebin the raw distribution for statistics
    TH2F  *h2rawR=RebinHisto(h2raw);
    
    //FILL MC HISTOGRAMS  **************************************************
    nEv=mc->GetEntries(); 
    Printf ("nEv = %d", nEv);
    
    mc->SetBranchAddress("partonCode", &partonCode); 
    mc->SetBranchAddress("ptJet", &ptJet); 
    mc->SetBranchAddress("ptDJet", &ptDJet); 
    mc->SetBranchAddress("angularity", &angularity); 
    mc->SetBranchAddress("lesub", &lesub); 
    mc->SetBranchAddress("ptJetMatch", &ptJetMatch); 
    mc->SetBranchAddress("ptDJetMatch", &ptDJetMatch); 
    mc->SetBranchAddress("angularityMatch", &angularityMatch); 
    mc->SetBranchAddress("lesubMatch", &lesubMatch); 
    mc->SetBranchAddress("weightPythiaFromPtHard", &weightPythiaFromPtHard);
    
    
    
    for(int iEntry=0; iEntry< nEv; iEntry++){
    	mc->GetEntry(iEntry); 
    	
        if(ptJetMatch>160) continue;
     	h2fulleff->Fill(ptDJetMatch,ptJetMatch,weightPythiaFromPtHard);  
        h2smearednocuts->Fill(ptDJet,ptJet,weightPythiaFromPtHard);
        if(ptJet>80 || ptJet<20) continue;
        
        if(ptDJet<0.3) continue;
        
		h2smeared->Fill(ptDJet,ptJet,weightPythiaFromPtHard);
		h2true->Fill(ptDJetMatch,ptJetMatch,weightPythiaFromPtHard);
	}
	
	TH1F *htrueptd=(TH1F*) h2fulleff->ProjectionX("trueptd",1,h2fulleff->GetNbinsY());
	TH1F *htruept=(TH1F*) h2fulleff->ProjectionY( "truept",1,h2fulleff->GetNbinsX()); 
	TH1F *hpttt=(TH1F*)htrueptd->Clone("hpttt");
	//rebin the smeared distribution 
    TH1F *hsmearpt=(TH1F*) h2smeared->ProjectionY( "hsmearpt",1,h2smeared->GetNbinsX()); 
    TH2F *h2smearedR=RebinHisto(h2smeared);
    TH2F *h2trueR=RebinHisto1D(h2true);
    
    
    
    
    RooUnfoldResponse response;
    RooUnfoldResponse responsenotrunc;
    RooUnfoldResponse responsenew;
    
    response->Setup(h2smearedR,h2trueR);
    responsenotrunc->Setup(h2smearednocuts,h2fulleff);
    responsenew->Setup(h2smearedR,h2trueR);
    
    //FILL THE RESPONSE  **************************************************************
    for(int iEntry=0; iEntry< nEv; iEntry++){
    	mc->GetEntry(iEntry); 
    	if(ptJetMatch>160 ) continue;
    	int bin1=htrueptd->FindFixBin(ptDJetMatch);
    	int bin2=htruept->FindFixBin(ptJetMatch);
    	
    	
    	responsenotrunc->Fill(ptDJet,ptJet,ptDJetMatch,ptJetMatch,weightPythiaFromPtHard);
    	
    	if(ptJet>80 || ptJet<20) continue;
    	if(ptDJet<0.3) continue;
    	
    	histomult[bin1-1][bin2-1]->Fill(ptDJet,ptJet,weightPythiaFromPtHard);
    	
    	
    	response->Fill(ptDJet,ptJet,ptDJetMatch,ptJetMatch,weightPythiaFromPtHard);
    	
    	
    }
    
    
    //SMEAR THE PRIOR CORRELATION AND FILL A DIFFERENT RESPONSE  ****************************************
    //we start with (shape_true,pT_true,shape_det,pT_det)
    //we smear shape_true according to a gaussian of given width to obtain shape_true'
    //Then, I consider the 1D histogram for detector level shapes corresponding to (shape_true',pT_true,pT_det) and
    //generate shape_det' randomly. 
    //Last, I fill the response with (shape_true',pT_true,shape_det',pT_det)
    
    
    
    double xx,yy;
    for(int iEntry=0; iEntry< nEv; iEntry++){
    	mc->GetEntry(iEntry); 
    	if(ptJetMatch>160 ) continue;
    	double smear=0;
    	//ptDJetMatch is the particle level shape, I smear it with a gaussian
    	smear = fRandom->Gaus(ptDJetMatch,0.01*frac*ptDJetMatch);
    	
    	h2fulleffsm->Fill(smear,ptJetMatch,weightPythiaFromPtHard);
    	//then I check the bins of  (shape_True',pt_True,pt_Det)
    	int bin1=htrueptd->FindFixBin(TMath::Abs(smear));
    	int bin2=htruept->FindFixBin(ptJetMatch);
    	int bin3=hsmearpt->FindFixBin(ptJet);
    	
        if(ptJet>80 || ptJet<20) continue;
        
        //then I generated randomly    shape_Det according to the 1D histogram 
        TH1F *histomult1d=(TH1F*) histomult[bin1-1][bin2-1]->ProjectionX("histomult1d",bin3,bin3); 
        
        
        double xx=histomult1d->GetRandom(); 
        
        if(xx<0.3) continue;  
        
        h2truesm->Fill(smear,ptJetMatch,weightPythiaFromPtHard);  
        //and fill the new response
        responsenew->Fill(xx,ptJet,smear,ptJetMatch,weightPythiaFromPtHard);
        
        
    }
    
    
    
    
    
    
    
    TH2D* hfold=(TH2D*)h2rawR->Clone("hfold");
    hfold->Sumw2();
    //compute kinematic efficiencies////
    TH1D * effok=(TH1D *)h2truesm->ProjectionX("effok",2,2);
    TH1D * effok1=(TH1D *)h2fulleffsm->ProjectionX("effok2",2,2);
    effok->Divide(effok1);
    effok->SetName("correff20-40");
    
    TH1D * effok3=(TH1D *)h2truesm->ProjectionX("effok3",3,3);
    TH1D * effok4=(TH1D *)h2fulleffsm->ProjectionX("effok4",3,3);
    effok3->Divide(effok4);
    effok3->SetName("correff40-60"); 
    
    TH1D * effok5=(TH1D *)h2truesm->ProjectionX("effok5",4,4);
    TH1D * effok6=(TH1D *)h2fulleffsm->ProjectionX("effok6",4,4);
    effok5->Divide(effok6);
    effok5->SetName("correff60-80"); 
    
    
    
    
    TH1D * trueptd=(TH1D *)h2fulleff->ProjectionX("trueptd",2,3);
    
    TFile *fout=new TFile (Form("shapesunfoldpTDR%iEff%iPrior%iRandomType%i.root",rad,eff,frac,gaus),"RECREATE");
    
    
    for(int jar=1;jar<10;jar++){
    	Int_t iter=jar;
        for(Int_t nn=0;nn<1000;nn++){
        	
        	if(nn>0)  TH2F *h2rawRSmeared=SmearPoints(h2rawR);
        	if(nn==0)  TH2F *h2rawRSmeared=(TH2F*)h2rawR->Clone("h2rawRSmeared");
        	cout<<"iteration"<<iter<<endl;
        	cout<<"==============Unfold h1====================="<<endl;
        	
        	RooUnfoldBayes   unfold(&responsenew, h2rawRSmeared, iter);    // OR
        	TH2D* hunf= (TH2D*) unfold.Hreco(errorTreatment);
        	//FOLD BACK
        	
        	
        	
        	
        	for(int i=0;i<6;i++){
        		for(int j=0;j<5;j++){
        			double effects=0;
        			
        			for(int k=0;k<9;k++){
        				for(int l=0;l<8;l++){
        					
        					int indexm=i+6*j;
        					int indext=k+9*l;
        					
        					
        					
        					effects=effects+hunf->GetBinContent(k+1,l+1)*responsenew(indexm,indext);
        					
        				}
        			}
        			hfold->SetBinContent(i+1,j+1,effects);
        		}
        	}
        	
        	
        	
        	
        	TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
        	htempUnf->SetName(Form("Bayesian_Unfolded_Cent0020_TT2050-1520_R0.4_Cut150MeV_iter%dRandom%i.root",iter,nn));
        	
        	TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
        	htempFold->SetName(Form("Bayesian_Folded_Cent0020_TT2050-1520_R0.4_Cut150MeV_iter%dRandom%i.root",iter,nn));        
        	
        	
        	
        	
        	
        	
        	fout->cd();
        	htempUnf->Write();
        	htempFold->Write();
        	
        	
        	
        	
        	
        }
    }
    h2truepostrandom->Write();
    h2trueprerandom->Write();
    h2shapepostrandom->Write();
    effok->Write();  
    effok3->Write();
    effok5->Write();
    h2rawR->SetName("raw");
    h2rawR->Write();
    h2smearedR->SetName("smeared");
    h2smearedR->Write();
    trueptd->Write();
    
    
    
}

















void Normalize2D(TH2* h)
{
	Int_t nbinsYtmp = h->GetNbinsY();
	const Int_t nbinsY = nbinsYtmp;
	Double_t norm[nbinsY];
	for(Int_t biny=1; biny<=nbinsY; biny++)
	{
		norm[biny-1] = 0;
		for(Int_t binx=1; binx<=h->GetNbinsX(); binx++)
		{
			norm[biny-1] += h->GetBinContent(binx,biny);
		}
	}
	
	for(Int_t biny=1; biny<=nbinsY; biny++)
	{
		for(Int_t binx=1; binx<=h->GetNbinsX(); binx++)
		{
			if(norm[biny-1]==0)  continue;
			else
			{
				h->SetBinContent(binx,biny,h->GetBinContent(binx,biny)/norm[biny-1]);
				h->SetBinError(binx,biny,h->GetBinError(binx,biny)/norm[biny-1]);
			}
		}
	}
}

TH1D *TruncateHisto(TH1D *gr,Int_t nbinsold,Int_t lowold,Int_t highold,Int_t nbinsnew,Int_t lownew,Int_t highnew,Int_t lim)
{
	TH1D *hTruncate=new TH1D("hTruncate","",nbinsnew,lownew,highnew);
	hTruncate->Sumw2();
	for(Int_t binx = 1; binx <= hTruncate->GetNbinsX(); binx++){
		hTruncate->SetBinContent(binx,gr->GetBinContent(lim+binx));    
	hTruncate->SetBinError(binx,gr->GetBinError(lim+binx));}   
	
	
	
	
	return hTruncate;
}




TH2F *RebinHisto(TH2F *gr)
{
	
	
	Double_t xbins[6];
	xbins[0]=20;
	xbins[1]=30;
	xbins[2]=40;
	xbins[3]=50;
	xbins[4]=60;
	xbins[5]=80;
	
	
    
	
	Double_t xbinsb[7];
	
	xbinsb[0]=0.3;
	xbinsb[1]=0.4;
	xbinsb[2]=0.5;
	xbinsb[3]=0.6;
	xbinsb[4]=0.7;
	xbinsb[5]=0.8;
	xbinsb[6]=1;
	
	
	
	
	
    TH2F *cucu=new TH2F("cucu","",6,xbinsb,6,20,80);    
    cucu->Sumw2();
    TH1D * cucuy= cucu->ProjectionY("cucuy",1,cucu->GetNbinsX(),"");
    TH1D * cucux= cucu->ProjectionX("cucux",1,cucu->GetNbinsY(),"");
    
    
    for(Int_t j=1;j<=6;j++){
    	Double_t content=0;
    	
    	for(Int_t m=1;m<=5;m++){
    		content=gr->GetBinContent(m,j);
    		cucu->SetBinContent(m,j,content);
    	}    
    	content=0;
    	for(Int_t m=6;m<=7;m++){
    	content=content+gr->GetBinContent(m,j);}
    	cucu->SetBinContent(6,j,content);
    	content=0;
    	
    }
    
    TH2F *cucu2=new TH2F("cucu2","",6,xbinsb,5,xbins);
    
    for(Int_t ll=1;ll<=6;ll++){       		  
    	Double_t content=0;
    	for(Int_t nn=1;nn<=4;nn++){
    		content=cucu->GetBinContent(ll,nn);
    	cucu2->SetBinContent(ll,nn,content);}
    	
    	content=0;
    	for(Int_t nn=5;nn<=6;nn++){
    	content=content+cucu->GetBinContent(ll,nn);}  
    	cucu2->SetBinContent(ll,5,content);
    }
    
    return cucu2;
}

TH2F *RebinHisto1D(TH2F *gr)
{
	
	
	Double_t xbins[10];
	xbins[0]=0;
	xbins[1]=0.1;
	xbins[2]=0.2;
	xbins[3]=0.3;
	xbins[4]=0.4;
	xbins[5]=0.5;
	xbins[6]=0.6;
	xbins[7]=0.7;
	xbins[8]=0.8;
	xbins[9]=1.;
	
	
	TH2F *cucu2=new TH2F("cucu2","",9,xbins,8,0,160);
	
	for(Int_t ll=1;ll<=8;ll++){       		  
		Double_t content=0;
		for(Int_t nn=1;nn<=8;nn++){
			content=gr->GetBinContent(nn,ll);
		cucu2->SetBinContent(nn,ll,content);}
		
		content=0;
		for(Int_t nn=9;nn<=10;nn++){
		content=content+gr->GetBinContent(nn,ll);}  
		cucu2->SetBinContent(9,ll,content);
	}
	
	return cucu2;
}









TH1D *RebinHisto3(TH1 *gr)
{
	
	Double_t xbinsb[8];
	xbinsb[0]=0;
	xbinsb[1]=20;
	xbinsb[2]=36;
	xbinsb[3]=52;
	xbinsb[4]=68;
	xbinsb[5]=84;
	xbinsb[6]=100;
	xbinsb[7]=200;
    
	
	
	
	TH1D *hRawnew3=(TH1D*)gr->Rebin(7,"hRawnew3",xbinsb);
	Double_t content=0;
	Double_t scale1=0;
	Double_t erre=0;
	for(Int_t cu=1;cu<=hRawnew3->GetNbinsX();cu++){
		scale1=hRawnew3->GetBinWidth(cu);
		content=hRawnew3->GetBinContent(cu);
		erre=hRawnew3->GetBinError(cu);
		hRawnew3->SetBinContent(cu,content*(2./scale1));
		hRawnew3->SetBinError(cu,erre*(2./scale1));
	}
	
	
	
	
	return hRawnew3;
}
















TH1D *RebinHistoMany(TH1D *gr, Int_t radius,Int_t xtrunc){
	
    
	if(xtrunc==1) {Double_t xbins[14]; //12
	Int_t nn=14;}
	if(xtrunc==2){Double_t xbins[13]; //14
	Int_t nn=13;}
	if(xtrunc==3){Double_t xbins[12]; //16
	Int_t nn=12;}
	if(xtrunc==4){Double_t xbins[11]; //18
	Int_t nn=11;}
	if(xtrunc==5) {Double_t xbins[9]; //24
	Int_t nn=9;}
	if(xtrunc==6){Double_t xbins[8]; //30
	Int_t nn=8;}
	if(xtrunc==7){Double_t xbins[7]; //36
	Int_t nn=7;}
	
	Double_t xbinsa[10];
	
	xbinsa[0]=20;
	xbinsa[1]=24;
	xbinsa[2]=28;
	xbinsa[3]=34;
	xbinsa[4]=42;
	xbinsa[5]=54;
	xbinsa[6]=64;
	xbinsa[7]=74;
	xbinsa[8]=84;
	xbinsa[9]=100; 
	
	
	
	
	
	Int_t re=10-nn;
	if(nn<10){
		for(Int_t la=0;la<nn;la++){
			xbins[la]=xbinsa[la+re];
	cout<<la<<" "<<xbins[la]<<" "<<xbinsa[la+re]<<endl;}}
	
	
	
	if(nn>10){    
		Int_t re=-1*re;
		for(Int_t la=0;la<re;la++){
			xbins[la]=20-2*(re-la);
		cout<<xbins[la]<<endl;}
		for(Int_t la=0;la<10;la++){
			xbins[la+re]=xbinsa[la];
	cout<<la<<" "<<xbins[la+re]<<" "<<xbinsa[la]<<endl;}}
    
	TH1D *hRawnew2=(TH1D*)gr->Rebin(nn-1,"hRawnew2",xbins);
	return hRawnew2;
	
}



TH2F *RebinResponseMany(TH2D *gr, Int_t radius,Int_t xtrunc)
{
	
	if(xtrunc==1) {Double_t xbins[14]; //12
	Int_t nn=14;}
	if(xtrunc==2){Double_t xbins[13]; //14
	Int_t nn=13;}
	if(xtrunc==3){Double_t xbins[12]; //16
	Int_t nn=12;}
	if(xtrunc==4){Double_t xbins[11]; //18
	Int_t nn=11;}
	if(xtrunc==5) {Double_t xbins[9]; //24
	Int_t nn=9;}
	if(xtrunc==6){Double_t xbins[8]; //28
	Int_t nn=8;}
	if(xtrunc==7){Double_t xbins[7]; //34
	Int_t nn=7;}
	
	Double_t xbinsb[8];
	xbinsb[0]=0;
	xbinsb[1]=20;
	xbinsb[2]=36;
	xbinsb[3]=52;
	xbinsb[4]=68;
	xbinsb[5]=84;
	xbinsb[6]=100;
	xbinsb[7]=200;
    
	
	
	
	
	
	
	Double_t xbinsa[10];
	
	xbinsa[0]=20;
	xbinsa[1]=24;
	xbinsa[2]=28;
	xbinsa[3]=34;
	xbinsa[4]=42;
	xbinsa[5]=54;
	xbinsa[6]=64;
	xbinsa[7]=74;
	xbinsa[8]=84;
	xbinsa[9]=100; 
	
	
	
	
	Int_t re=10-nn;
	if(nn<10){
		for(Int_t la=0;la<nn;la++){
			xbins[la]=xbinsa[la+re];
	cout<<la<<" "<<xbins[la]<<" "<<xbinsa[la+re]<<endl;}}
	
	
	
	if(nn>10){    
		Int_t re=-1*re;
		for(Int_t la=0;la<re;la++){
			xbins[la]=20-2*(re-la);
		cout<<xbins[la]<<endl;}
		for(Int_t la=0;la<10;la++){
			xbins[la+re]=xbinsa[la];
	cout<<la<<" "<<xbins[la+re]<<" "<<xbinsa[la]<<endl;}}
	
	
	
	
	TH2F *cucu=new TH2F("cucu","",nn-1,xbins,100,0.,200.);    
	TH1D * cucuy= cucu->ProjectionY("cucuy",1,cucu->GetNbinsX(),"");
	TH1D * cucux= cucu->ProjectionX("cucux",1,cucu->GetNbinsY(),"");
	
    
	
	for(Int_t j=1;j<=100;j++){
		Double_t content=0;
		
		if(nn>10){
			for(Int_t la=1;la<=re;la++){
			cucu->SetBinContent(la,j,gr->GetBinContent(la,j));}
			for(Int_t m=1;m<=2;m++){
			content=content+gr->GetBinContent(m+re,j);}     
			cucu->SetBinContent(1+re,j,content);
			content=0; 
			for(Int_t m=3;m<=4;m++){
			content=content+gr->GetBinContent(m+re,j);}     
			cucu->SetBinContent(2+re,j,content);
			content=0; 
			for(Int_t m=5;m<=7;m++){
			content=content+gr->GetBinContent(m+re,j);}     
			cucu->SetBinContent(3+re,j,content);
			content=0; 
			
			for(Int_t m=8;m<=11;m++){
			content=content+gr->GetBinContent(m+re,j);}     
			cucu->SetBinContent(4+re,j,content);
			content=0; 
			
			for(Int_t m=12;m<=17;m++){
			content=content+gr->GetBinContent(m+re,j);}     
			cucu->SetBinContent(5+re,j,content);
			content=0; 
			
			for(Int_t m=0;m<3;m++){
				content=0;
				for(Int_t n=1;n<=5;n++){
				content=content+gr->GetBinContent(re+17+m*5+n,j);}
			cucu->SetBinContent(re+m+6,j,content);}
			content=0;
			
			for(Int_t n=33;n<=40;n++){
			content=content+gr->GetBinContent(n+re,j);}
		cucu->SetBinContent(9+re,j,content);}
		
		
		
		
		
		
		
		
		if(nn==9){
			content=0;
			
			
			for(Int_t m=3;m<=4;m++){
			content=content+gr->GetBinContent(m-2,j);}     
			cucu->SetBinContent(1,j,content);
			content=0; 
			for(Int_t m=5;m<=7;m++){
			content=content+gr->GetBinContent(m-2,j);}     
			cucu->SetBinContent(2,j,content);
			content=0; 
			
			for(Int_t m=8;m<=11;m++){
			content=content+gr->GetBinContent(m-2,j);}     
			cucu->SetBinContent(3,j,content);
			content=0; 
			
			for(Int_t m=12;m<=17;m++){
			content=content+gr->GetBinContent(m-2,j);}     
			cucu->SetBinContent(4,j,content);
			content=0; 
			
			for(Int_t m=0;m<3;m++){
				content=0;
				for(Int_t n=1;n<=5;n++){
				content=content+gr->GetBinContent(15+m*5+n,j);}
			cucu->SetBinContent(m+5,j,content);}
			content=0;
			
			for(Int_t n=33;n<=40;n++){
			content=content+gr->GetBinContent(n-2,j);}
		cucu->SetBinContent(8,j,content);}
		
		
		
		
		
		
		if(nn==8){
			content=0;
			
			for(Int_t m=5;m<=7;m++){
			content=content+gr->GetBinContent(m-4,j);}     
			cucu->SetBinContent(1,j,content);
			content=0; 
			
			for(Int_t m=8;m<=11;m++){
			content=content+gr->GetBinContent(m-4,j);}     
			cucu->SetBinContent(2,j,content);
			content=0; 
			
			for(Int_t m=12;m<=17;m++){
			content=content+gr->GetBinContent(m-4,j);}     
			cucu->SetBinContent(3,j,content);
			content=0; 
			
			for(Int_t m=0;m<3;m++){
				content=0;
				for(Int_t n=1;n<=5;n++){
				content=content+gr->GetBinContent(13+m*5+n,j);}
			cucu->SetBinContent(m+4,j,content);}
			content=0;
			
			for(Int_t n=33;n<=40;n++){
			content=content+gr->GetBinContent(n-4,j);}
		cucu->SetBinContent(7,j,content);}
		
		
		
		
		
		
		if(nn==7){
			content=0;
			
			
			for(Int_t m=8;m<=11;m++){
			content=content+gr->GetBinContent(m-7,j);}     
			cucu->SetBinContent(1,j,content);
			content=0; 
			
			for(Int_t m=12;m<=17;m++){
			content=content+gr->GetBinContent(m-7,j);}     
			cucu->SetBinContent(2,j,content);
			content=0; 
			
			for(Int_t m=0;m<3;m++){
				content=0;
				for(Int_t n=1;n<=5;n++){
				content=content+gr->GetBinContent(10+m*5+n,j);}
			cucu->SetBinContent(m+3,j,content);}
			content=0;
			
			for(Int_t n=33;n<=40;n++){
			content=content+gr->GetBinContent(n-7,j);}
		cucu->SetBinContent(6,j,content);}
		
		
		
		
		
		
		
	}
	
	
	
	TH2F *cucu2=new TH2F("cucu2","",nn-1,xbins,7,xbinsb);
	
	for(Int_t ll=1;ll<=nn-1;ll++){
		
		
		content=0;
		for(Int_t k=1;k<11;k++){
		content=content+cucu->GetBinContent(ll,k);}
		cucu2->SetBinContent(ll,1,content);
		
		
		
		for(Int_t mm=0;mm<4;mm++){
		    content=0;
		    for(Int_t k=1;k<=8;k++){
		    content=content+cucu->GetBinContent(ll,10+8*mm+k);}  
		cucu2->SetBinContent(ll,mm+2,content);}
		content=0;
		for(Int_t k=43;k<=50;k++){
		content=content+cucu->GetBinContent(ll,k);}
		cucu2->SetBinContent(ll,6,content);
		content=0;
		
		for(Int_t k=51;k<=100;k++){
		content=content+cucu->GetBinContent(ll,k);}
	cucu2->SetBinContent(ll,7,content);}		 
	
	
	
	
	
	
	return cucu2;
	
	
}




TH1D *RebinHistoOther(TH1D *gr, Int_t radius)
{
	
	Double_t xbins[9];
	xbins[0]=20;
	xbins[1]=28;
	xbins[2]=36; 
	for(Int_t n=1;n<=5;n++){
	xbins[n+2]=36+n*8;}
	xbins[8]=110;
	
	
	TH1D *hRawnewo1=(TH1D*)gr->Rebin(8,"hRawnewo1",xbins);
	return hRawnewo1;
	
}

TH2F *RebinResponseOther(TH2D *gr, Int_t radius)
{
	
	Double_t xbinsb[8];
	xbinsb[0]=0;
	xbinsb[1]=20;
	xbinsb[2]=36;
	xbinsb[3]=52;
	xbinsb[4]=68;
	xbinsb[5]=84;
	xbinsb[6]=100;
	xbinsb[7]=200;
	
	
	
	Double_t xbins[9];
	xbins[0]=20;
	xbins[1]=28;
	xbins[2]=36; 
	for(Int_t n=1;n<=5;n++){
	xbins[n+2]=36+n*8;}
	xbins[8]=110;
	
	
	
	TH2F *cucu=new TH2F("cucu","",8,xbins,100,0,200.);    
	TH1D * cucuy= cucu->ProjectionY("cucuy",1,cucu->GetNbinsX(),"");
	TH1D * cucux= cucu->ProjectionX("cucux",1,cucu->GetNbinsY(),"");
	Double_t content=0;
	
	for(Int_t j=1;j<=100;j++){
		
		for(Int_t m=0;m<7;m++){
			content=0;
			for(Int_t n=1;n<=4;n++){
			content=content+gr->GetBinContent(4*m+n,j);}
		cucu->SetBinContent(m+1,j,content);}
		content=0;
		
		for(Int_t n=29;n<=45;n++){
		content=content+gr->GetBinContent(n,j);}
	cucu->SetBinContent(8,j,content); }
	
	
  	TH2F *cucu2=new TH2F("cucu2","",8,xbins,7,xbinsb);
  	
  	for(Int_t ll=1;ll<=8;ll++){
  		
  		content=0;
  		for(Int_t nn=1;nn<11;nn++){
  		content=content+cucu->GetBinContent(ll,nn);}
  		cucu2->SetBinContent(ll,1,content);
  		
  		
  		for(Int_t mm=0;mm<4;mm++){
		    content=0;
		    for(Int_t nn=1;nn<=8;nn++){
		    content=content+cucu->GetBinContent(ll,10+8*mm+nn);}  
		cucu2->SetBinContent(ll,mm+2,content);}
		content=0;
		for(Int_t nn=43;nn<=50;nn++){
		content=content+cucu->GetBinContent(ll,nn);}
		cucu2->SetBinContent(ll,6,content);
		content=0;
		
		for(Int_t nn=51;nn<=100;nn++){
		content=content+cucu->GetBinContent(ll,nn);}
	cucu2->SetBinContent(ll,7,content);}  		 
return cucu2;}




TH1D *RebinHistoOther2(TH1D *gr, Int_t radius)
{
	
	Double_t xbins[7];
	
	
	xbins[0]=20;
	xbins[1]=30;
	xbins[2]=42;
	xbins[3]=52;
	xbins[4]=64;
	xbins[5]=78;
	xbins[6]=94;
	
	
	
	TH1D *hRawnew=(TH1D*)gr->Rebin(6,"hRawnew",xbins);
	
	
	return hRawnew;
}


TH2F *RebinResponseOther2(TH2D *gr, Int_t radius)
{
	
	Double_t xbinsb[8];
	xbinsb[0]=0;
	xbinsb[1]=20;
	xbinsb[2]=36;
	xbinsb[3]=52;
	xbinsb[4]=68;
	xbinsb[5]=84;
	xbinsb[6]=100;
	xbinsb[7]=200;
	
	
	Double_t xbins[7];
	
	xbins[0]=20;
	xbins[1]=30;
	xbins[2]=42;
	xbins[3]=52;
	xbins[4]=64;
	xbins[5]=78;
	xbins[6]=94;
	
	
	
	
	
	TH2F *cucu=new TH2F("cucu","",6,xbins,100,0.,200.);
	TH1D * cucuy= cucu->ProjectionY("cucuy",1,cucu->GetNbinsX(),"");
	TH1D * cucux= cucu->ProjectionX("cucux",1,cucu->GetNbinsY(),"");
	
	
	
	
	
	for(Int_t j=1;j<=100;j++){
		
		Double_t content=0;
		for(Int_t m=1;m<=5;m++){
		content=content+gr->GetBinContent(m,j);}     
		cucu->SetBinContent(1,j,content);
		content=0; 
		
		for(Int_t m=1;m<=6;m++){
		content=content+gr->GetBinContent(5+m,j);}
		cucu->SetBinContent(2,j,content);
		content=0;
		
		for(Int_t n=12;n<=16;n++){
		content=content+gr->GetBinContent(n,j);}
		cucu->SetBinContent(3,j,content);
		content=0;
		for(Int_t n=17;n<=22;n++){
		content=content+gr->GetBinContent(n,j);}
		cucu->SetBinContent(4,j,content);
		content=0;
		
		for(Int_t n=23;n<=29;n++){
		content=content+gr->GetBinContent(n,j);}
		cucu->SetBinContent(5,j,content);
		content=0;
		
		for(Int_t n=30;n<=37;n++){
		content=content+gr->GetBinContent(n,j);}
		cucu->SetBinContent(6,j,content);
	}
    
	
	
	
	
	
	TH2F *cucu2=new TH2F("cucu2","",6,xbins,7,xbinsb);
	
	for(Int_t ll=1;ll<=6;ll++){
		
		
		content=0;
		for(Int_t nn=1;nn<11;nn++){
		content=content+cucu->GetBinContent(ll,nn);}
		cucu2->SetBinContent(ll,1,content);
		
		for(Int_t mm=0;mm<4;mm++){
		    content=0;
		    for(Int_t nn=1;nn<=8;nn++){
		    content=content+cucu->GetBinContent(ll,10+8*mm+nn);}  
		cucu2->SetBinContent(ll,mm+2,content);}
		content=0;
		for(Int_t nn=43;nn<=50;nn++){
		content=content+cucu->GetBinContent(ll,nn);}
		cucu2->SetBinContent(ll,6,content);
		content=0;
		
		for(Int_t nn=51;nn<=100;nn++){
		content=content+cucu->GetBinContent(ll,nn);}
	cucu2->SetBinContent(ll,7,content);}  		 
	
	return cucu2;
	
	
	
}



TH2F *SmearPoints(TH2F *gr)
{
    
	TH2F *hnew(0);
	fRandom = new TRandom3(0);
    hnew=(TH2F*)gr->Clone("hnew");  
    TH1F *hnewx=(TH1F*) hnew->ProjectionX("hnewx",1,-1);
    TH1F *hnewy=(TH1F*) hnew->ProjectionY("hnewy",1,-1);
    for(Int_t i=1; i<=hnewx->GetNbinsX(); i++){
    	for(Int_t j=1; j<=hnewy->GetNbinsX(); j++){  
    		
    		Double_t contgr=gr->GetBinContent(i,j); 
    		Double_t errg=gr->GetBinError(i,j);
    		Double_t content=fRandom->Gaus(contgr,errg);
    		hnew->SetBinContent(i,j,content);
    		hnew->SetBinError(i,j,errg);
    }}
    return hnew;
}










#ifndef __CINT__
int main () { RooUnfoldmine(); return 0; }  // Main program when run stand-alone
#endif
