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
#include "TH4D.h"
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

void RooSimplepTD()
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif
  Int_t difference=1;
  Int_t Ppol=0;
  cout << "==================================== pick up the response matrix for background==========================" << endl;
  ///////////////////parameter setting
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  //Detector response from fast simulation
  //TString fnamea = "/Users/leticia/WorkdirQGTag/trainresultstagging0.2/AnalysisResults_TrueDet_NoSub.root";
  // TList *list2;
  // TFile *input2 = TFile::Open( fnamea );
  //list2=(TList*) input2->Get("JetQGTaggings_JetMCOnly_AKTChargedR020_DetParticles_pT0150_E_scheme_TCRawConstSub_TrueDet_NoSub_Incl"); 
  // TTree *mc=(TTree*)list2->FindObject("fTreeJetShape");
  
  //Get the tree for MC
   TString fname = "/Users/leticia/WorkdirQGTag/trainresultstagging0.2/Resultspp/DetectorResponse/DetRespPythia.root";
   TList *list1;
   TFile *input = TFile::Open( fname );
   list1=(TList*) input->Get("JetQGTaggings_Jet_AKTChargedR020_PicoTracks_pT0150_E_scheme_TCRaw_PythiaDef_NoSub_Incl"); 
   TTree *mc=(TTree*)list1->FindObject("newtree");
  
   //Get the tree for data
  TString fname2 = "/Users/leticia/WorkdirQGTag/trainresultstagging0.2/AnalysisResults_Datapp_NoSub.root";
   TList *list2;
   TFile *input2 = TFile::Open( fname2 );
   list2=(TList*) input2->Get("JetQGTaggings_Jet_AKTChargedR020_PicoTracks_pT0150_E_scheme_TCRaw_Data_NoSub_Incl"); 
   TTree *data=(TTree*)list2->FindObject("fTreeJetShape");
 
  




  TH2F *h2raw(0);
   h2raw=new TH2F("raw","raw",7,0.3,1.,9,10,100);
  TH2F *h2smeared(0);
   h2smeared=new TH2F("smeared","smeared",7,0.3,1,9,10,100);
  TH2F *h2true(0);
   h2true=new TH2F("true","true",10,0.,1.,14,10,140);
  TH2F *h2fulleff(0);
  h2fulleff=new TH2F("truef","truef",10,0.,1.,14,10,140);
   TH2F *h2fulleffFast(0);
  h2fulleffFast=new TH2F("truefFast","truefFast",10,0.,1.,14,10,140);
   RooUnfoldResponse response;
   response->Setup(h2smeared,h2true);
   h2smeared->Sumw2();
   h2true->Sumw2();
   h2raw->Sumw2();
  Float_t partonCode=0., ptJet=0., ptDJet=0., mJet=0., nbOfConst=0., angularity=0., circularity=0., lesub=0., sigma2=0.;
  Float_t ptJetMatch=0., ptDJetMatch=0., mJetMatch=0., nbOfConstMatch=0., angularityMatch=0., circularityMatch=0., lesubMatch=0., sigma2Match=0.,weightPythia=0.; 
 
  Float_t partonCode2=0., ptJet2=0., ptDJet2=0., mJet2=0., nbOfConst2=0., angularity2=0., circularity2=0., lesub2=0., sigma22=0.;
  Float_t ptJetMatch2=0., ptDJetMatch2=0., mJetMatch2=0., nbOfConstMatch2=0., angularityMatch2=0., circularityMatch2=0., lesubMatch2=0., sigma2Match2=0.,weightPythia2=0.; 
 
 Int_t nEv=mc->GetEntries(); 
  Printf ("nEv = %d", nEv);
  mc->SetBranchAddress("partonCode", &partonCode); 
  mc->SetBranchAddress("ptJet", &ptJet); 
  mc->SetBranchAddress("ptDJet", &ptDJet); 
  mc->SetBranchAddress("mJet", &mJet);
  mc->SetBranchAddress("nbOfConst", &nbOfConst);
  mc->SetBranchAddress("angularity", &angularity); 
  mc->SetBranchAddress("circularity", &circularity);
  mc->SetBranchAddress("lesub", &lesub); 
  mc->SetBranchAddress("sigma2", &sigma2);
  mc->SetBranchAddress("ptJetMatch", &ptJetMatch); 
  mc->SetBranchAddress("ptDJetMatch", &ptDJetMatch); 
  mc->SetBranchAddress("mJetMatch", &mJetMatch);
  mc->SetBranchAddress("nbOfConstMatch", &nbOfConstMatch);
  mc->SetBranchAddress("angularityMatch", &angularityMatch); 
  mc->SetBranchAddress("circularityMatch", &circularityMatch);
  mc->SetBranchAddress("lesubMatch", &lesubMatch); 
  mc->SetBranchAddress("sigma2Match", &sigma2Match);
  mc->SetBranchAddress("weightPythia", &weightPythia);


// nt_t nEv2=mc2->GetEntries(); 
 //  Printf ("nEv2 = %d", nEv2);
 //  mc2->SetBranchAddress("partonCode", &partonCode2); 
 //  mc2->SetBranchAddress("ptJet", &ptJet2); 
 //  mc2->SetBranchAddress("ptDJet", &ptDJet2); 
 //  mc2->SetBranchAddress("mJet", &mJet2);
 //  mc2->SetBranchAddress("nbOfConst", &nbOfConst2);
 //  mc2->SetBranchAddress("angularity", &angularity2); 
 //  mc2->SetBranchAddress("circularity", &circularity2);
 //  mc2->SetBranchAddress("lesub", &lesub2); 
 //  mc2->SetBranchAddress("sigma2", &sigma22);
 //  mc2->SetBranchAddress("ptJetMatch", &ptJetMatch2); 
 //  mc2->SetBranchAddress("ptDJetMatch", &ptDJetMatch2); 
 //  mc2->SetBranchAddress("mJetMatch", &mJetMatch2);
 //  mc2->SetBranchAddress("nbOfConstMatch", &nbOfConstMatch2);
 //  mc2->SetBranchAddress("angularityMatch", &angularityMatch2); 
 //  mc2->SetBranchAddress("circularityMatch", &circularityMatch2);
 //  mc2->SetBranchAddress("lesubMatch", &lesubMatch2); 
 //  mc2->SetBranchAddress("sigma2Match", &sigma2Match2);
 //  mc2->SetBranchAddress("weightPythia", &weightPythia2);



  TH1F *htrueptd=(TH1F*) h2fulleff->ProjectionX("trueptd",1,-1);
  TH1F *htruept=(TH1F*) h2fulleff->ProjectionY( "truept",1,-1);

 
 
      for(int iEntry=0; iEntry< nEv; iEntry++){
       mc->GetEntry(iEntry); 
        if(ptJetMatch>140) continue;
        h2fulleff->Fill(ptDJetMatch,ptJetMatch);  
        if(ptJet>100 || ptJet<10) continue;
	if(ptDJet<0.3) continue;
       h2smeared->Fill(ptDJet,ptJet);
       h2true->Fill(ptDJetMatch,ptJetMatch);
        }

      

 

       for(int iEntry=0; iEntry< nEv; iEntry++){
       mc->GetEntry(iEntry); 
       if(ptJetMatch>140) continue;
     
       if(ptJet>100 || ptJet<10) continue;
         if(ptDJet<0.3) continue;
       int bin1=htrueptd->FindFixBin(ptDJetMatch);
       int bin2=htruept->FindFixBin(ptJetMatch);
       double norm=h2fulleff->GetBinContent(bin1,bin2);
       if(norm!=0) norm=1./norm;
      
       response->Fill(ptDJet,ptJet,ptDJetMatch,ptJetMatch,norm);
       //response is automatically normalized to one, no matter what weight you give to it. 
       //and the problem is that this is a truncated response.
 
   }
 

   
    
  
    //so mcr is correctly normalized to one, not the response.       
  cout<<"cucu"<<endl;
  nEv=data->GetEntries(); 
  cout<<"entries"<<nEv<<endl;
  data->SetBranchAddress("partonCode", &partonCode); 
  data->SetBranchAddress("ptJet", &ptJet); 
  data->SetBranchAddress("ptDJet", &ptDJet); 
  data->SetBranchAddress("mJet", &mJet);
  data->SetBranchAddress("nbOfConst", &nbOfConst);
  data->SetBranchAddress("angularity", &angularity); 
  data->SetBranchAddress("circularity", &circularity);
  data->SetBranchAddress("lesub", &lesub); 
  data->SetBranchAddress("sigma2", &sigma2);
  data->SetBranchAddress("ptJetMatch", &ptJetMatch); 
  data->SetBranchAddress("ptDJetMatch", &ptDJetMatch); 
  data->SetBranchAddress("mJetMatch", &mJetMatch);
  data->SetBranchAddress("nbOfConstMatch", &nbOfConstMatch);
  data->SetBranchAddress("angularityMatch", &angularityMatch); 
  data->SetBranchAddress("circularityMatch", &circularityMatch);
  data->SetBranchAddress("lesubMatch", &lesubMatch); 
  data->SetBranchAddress("sigma2Match", &sigma2Match);
  data->SetBranchAddress("weightPythia", &weightPythia);
  
  for(int iEntry=0; iEntry< nEv; iEntry++){
   data->GetEntry(iEntry); 
         if(ptJet>100 || ptJet<10) continue;
		 if(ptDJet<0.3) continue;
   
       h2raw->Fill(ptDJet,ptJet);
     
    

       
}
  
  TH2D* hfold=(TH2D*)h2raw->Clone("hfold");
 TH1D * truthfull=(TH1D *)h2fulleff->ProjectionX("truthfull",3,6);
 TH1D * truthtrunc=(TH1D *)h2true->ProjectionX("truthtrunc",3,6);
 truthtrunc->Divide(truthfull);
 
   TH1D * trueptd=(TH1D *)h2fulleff->ProjectionX("trueptd",3,6);
   TFile *fout=new TFile ("shapesunfoldpTD.root","RECREATE");
  fout->cd();
  truthtrunc->Write();
  h2raw->Write();
  h2smeared->Write();
  trueptd->Write();
     for(int jar=1;jar<15;jar++){
      Int_t iter=jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;

      RooUnfoldBayes   unfold(&response, h2raw, iter);    // OR
      TH2D* hunf= (TH2D*) unfold.Hreco(errorTreatment);
      //FOLD BACK
   



 for(int i=0;i<7;i++){
 for(int j=0;j<9;j++){
 double effects=0;
 double effects2=0;
 for(int k=0;k<10;k++){
 for(int l=0;l<14;l++){
   double effcorr=0;
   double effcorr1=h2true->GetBinContent(k+1,l+1);
   double effcorr2= h2fulleff->GetBinContent(k+1,l+1);
   if(effcorr1!=0) effcorr=effcorr2/effcorr1;
   else effcorr=0;
   int indexm=i+7*j;
   int indext=k+10*l;
  
   //effects=effects+effcorr*hunf->GetBinContent(k+1,l+1)*mcr[i][j][k][l];
   effects=effects+hunf->GetBinContent(k+1,l+1)*response(indexm,indext);
   //cout<<effects<<" "<<effects2<<endl;
}}
     hfold->SetBinContent(i+1,j+1,effects);}}

 

  
          TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
	  htempUnf->SetName(Form("Bayesian_Unfolded_Cent0020_TT2050-1520_R0.4_Cut150MeV_iter%d.root",iter));
           
	    TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
	  htempFold->SetName(Form("Bayesian_Folded_Cent0020_TT2050-1520_R0.4_Cut150MeV_iter%d.root",iter));        


           






     
      	  htempUnf->Write();
	  htempFold->Write();
	  

     

	  
     


     }}






   





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


   TH1 *RebinHisto(TH1 *gr, Int_t radius)
{

      Double_t xbins[10];

         xbins[0]=20;
         xbins[1]=24;
         xbins[2]=28;
 	 xbins[3]=34;
         xbins[4]=42;
         xbins[5]=54;
      	 xbins[6]=64;
         xbins[7]=74;
         xbins[8]=84;
         xbins[9]=100;



       
         TH1D *hRawnewOrg=(TH1D*)gr->Rebin(9,"hRawnewOrg",xbins);

         return hRawnewOrg;     
	
}

 TH2F *RebinResponse(TH2D *gr, Int_t radius)
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
    
         

         Double_t xbins[10];
        
         xbins[0]=20;
         xbins[1]=24;
         xbins[2]=28;
 	 xbins[3]=34;
         xbins[4]=42;
         xbins[5]=54;
      	 xbins[6]=64;
         xbins[7]=74;
         xbins[8]=84;
         xbins[9]=100; 

      
 

    TH2F *cucu=new TH2F("cucu","",9,xbins,100,0,200);    
    
           TH1D * cucuy= cucu->ProjectionY("cucuy",1,cucu->GetNbinsX(),"");
           TH1D * cucux= cucu->ProjectionX("cucux",1,cucu->GetNbinsY(),"");

     
     for(Int_t j=1;j<=100;j++){
     Double_t content=0;

     for(Int_t m=1;m<=2;m++){
     content=content+gr->GetBinContent(m,j);}     
     cucu->SetBinContent(1,j,content);
     content=0; 
     for(Int_t m=3;m<=4;m++){
     content=content+gr->GetBinContent(m,j);}     
     cucu->SetBinContent(2,j,content);
     content=0; 
     for(Int_t m=5;m<=7;m++){
     content=content+gr->GetBinContent(m,j);}     
     cucu->SetBinContent(3,j,content);
     content=0; 

     for(Int_t m=8;m<=11;m++){
     content=content+gr->GetBinContent(m,j);}     
     cucu->SetBinContent(4,j,content);
     content=0; 

     for(Int_t m=12;m<=17;m++){
     content=content+gr->GetBinContent(m,j);}     
     cucu->SetBinContent(5,j,content);
     content=0; 

     for(Int_t m=0;m<3;m++){
     content=0;
     for(Int_t n=1;n<=5;n++){
     content=content+gr->GetBinContent(17+m*5+n,j);}
     cucu->SetBinContent(m+6,j,content);}
     content=0;

     for(Int_t n=33;n<=40;n++){
     content=content+gr->GetBinContent(n,j);}
     cucu->SetBinContent(9,j,content);}




		TH2F *cucu2=new TH2F("cucu2","",9,xbins,7,xbinsb);

		for(Int_t ll=1;ll<=9;ll++){       		  

		  Double_t content=0;
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
		   cucu2->SetBinContent(ll,7,content);  		 


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
     












#ifndef __CINT__
int main () { RooUnfoldmine(); return 0; }  // Main program when run stand-alone
#endif
