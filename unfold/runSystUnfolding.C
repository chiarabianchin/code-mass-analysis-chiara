#include <TParameter.h>
#include "/data/Work/MyCodeJetMass/unfold/CreateRooUnfoldResponse.C"
#include "/data/Work/MyCodeJetMass/unfold/unfoldData.C"
#include "/data/Work/MyCodeJetMass/unfold/SystematicComparisonsUnf.C"
#include "/data/Work/MyCodeJetMass/macros/CalculateKineEff.C"
#include <iostream>
#include <fstream>
using std::endl;

// unfolding systematics
void runSystUnfolding(Bool_t style = kFALSE){

	if(style){
		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		gStyle->SetTextFont(42);
		
		gStyle->SetLabelSize(0.05, "X");
		gStyle->SetTitleOffset(1.1, "X");
		gStyle->SetTitleSize(0.048, "X");
		
		gStyle->SetLabelSize(0.05, "Y");
		gStyle->SetTitleOffset(1.34, "Y");
		gStyle->SetTitleSize(0.048, "Y");
		//gStyle->SetPadBottomMargin(.13);
		//gStyle->SetPadLeftMargin(.16);
		//gStyle->SetPadRightMargin(.06);
		
	}
	TString suff = "";
   
   // switches
   Int_t dryrun       = 0;     // explore directories, no actions
   Int_t domatrix     = 1;    // activate CreateRooUnfoldResponse. A check on the existance of the output file is performed, no action if the file is present
   Int_t dounfold     = 1;    // activate unfold and CombineTriggersCorrected.  A check on the existance of the output file is performed, no action if the file is present
   Int_t forcedounfold     = 0; //redo the unfolding even if the file is alredy present 
   if(forcedounfold)  dounfold = 1;
   Int_t kineeffinrangesyst = 0; // correct the variations of the range with kine eff
   Int_t reckineff = 0;
   Int_t correctforkineeff = 0;
   Bool_t fillmiss = 1;
   if(fillmiss == 1) {
   	   kineeffinrangesyst = 0;
   	   reckineff = 0;
   	   correctforkineeff = 0;
   }
   Int_t variablebinning = 1; //enable/disable the variation of variable binning
   Int_t docompare = 1;
   Int_t plotResults = 1;

   const Int_t ncontribution = 8;
   Int_t activecmp[ncontribution] = {1, 1, 1, 1, 0, 1, 0, 1}; // do = 1, don't do= 0
   TString description[ncontribution] = {"Bkgsub", "RespRange", "NIter", "OvlExclusion", "PtBinW", "TrkEff", "TrkEffv2", "Prior"};
   
   Int_t fittype[ncontribution] = {-1, 2, 1, 1, 1, 2, 2, -1};
   Int_t nactivecontributions = 0; 
   for(Int_t i=0; i<ncontribution; i++){
   	   if(activecmp[i]) nactivecontributions++;
   }
   TString fileSysOutput[nactivecontributions];
   TString hnameSysMean[nactivecontributions];
   // working directories
   TString basewdircmp = gSystem->pwd();
   Printf("The output of the systematics will be written in the current directory: %s", basewdircmp.Data());
   
   // the response matrices (CreateRooUnfoldResponse) and the unfolding (unfold) are performed in this base directories. Internal directories with the requested variations are created and the output files stored in each subdirectory
   const Int_t nunf = 12; //4;
   // 2MB and 2EJE
   TString wdirresp[nunf]={
    "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MB", 
    "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MB/ConstSub",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbOvlExclu/MB",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmb96TrEff/MB",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingVarPrior/MB",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/MB",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/EJE", 
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/EJE/ConstSub", 
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbOvlExclu/EJE",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmb96TrEff/EJE",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingVarPrior/EJE",
   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/EJE"
   };
   //,
   //
   
   TString wdirrespLowTrEff[2] = {
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train1289-1292/UnfoldEmb96TrEff/MB",
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train1289-1292/UnfoldEmb96TrEff/EJE"
   };
   
   // settings response matrix
   
   Int_t nresp = 6;
   
   TString strIn[nresp] = {
      //"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/RespRhoWeighted/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root",
      //"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/RespRhoWeighted/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", 
      "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root", 
      "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/AnalysisResultsResp1134.root", 
   	  "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbOvlExclu/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root",
   	  "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmb96TrEff/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root",
   	  "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingVarPrior/TreeResponse.root",
   	  "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root"
   }; 
   //,
   //	  
   TString strL[nresp] = {"fhResponseFinal", "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme","fhResponseFinal", "fhResponseFinal", "fhResponseFinal", "fhResponseFinal"}; // this is the name of the thnsparse with the weighted response matrix ,"fhResponseFinal"
   // set different name if the directory do not correspond to different subtraction methods
   TString tag[nresp] = {"DetFlBkgNoSub", "DetDataConstS", "NoOvl", "96TrkEff", "Prior", "VarBinWidth"}; //
   
   //"DetFlBkgDerivptpar20"
   //"DetFlBkgConstptpar20"
   
   // VARIATIONS requested for the systematic study
   
   const Int_t variationrec = 5, variationgen = 2;
   
   // bin corresponding to the "central value" or default
   Int_t centralrec = 0, centralgen = 0, centralvar = centralrec*variationgen + centralgen;
   Int_t totalvar = variationrec*variationgen;
   Int_t skipVariationRec[variationrec] = {0, 0, 0, 0, 0};
   Int_t skipVariationGen[variationgen] = {0, 1};
   //reco
   Double_t pt_min[variationrec] = {20, 10, 30, 20, 20}; //MB
   Double_t pt_max[variationrec] = {80, 80, 90, 140, 80}; //MB
   Double_t pt_mineje[variationrec] = {70, 60, 60, 70, 70};   //EJE
   Double_t pt_maxeje[variationrec] = {120, 120, 130, 140, 120}; //EJE
   
   Double_t m_min[variationrec] = {0., 0., 0., 0, 0};
   Double_t m_max[variationrec] = {12., 10., 14., 40, 12};
   Double_t m_mineje[variationrec] = {0., 0., 0., 0, 0};
   Double_t m_maxeje[variationrec] = {14., 12., 16., 40, 16};
   //gen
   Double_t pt_minT[variationgen] = {10., 10.};   //MB ->Changing this??
   Double_t pt_maxT[variationgen] = {150., 150.};  //MB
   
   Double_t pt_minTeje[variationgen] = {50., 50};   //EJE
   Double_t pt_maxTeje[variationgen] = {150., 140};  //EJE
   
   Double_t m_minT[variationgen] = {0., 0.};
   Double_t m_maxT[variationgen] = {40., 20.};
   
   //define the name of the subdirectories
   Int_t nvariations = 0;
   for(Int_t ir = 0; ir < variationrec; ir++){
   	   if(skipVariationRec[ir]) continue;
   	   for(Int_t ig = 0; ig < variationgen; ig++){
   	   	   if(skipVariationGen[ig]) continue;
   	   	   if(ir == centralrec && ig == centralgen) centralvar = nvariations;
   	   	   nvariations++;
   	   }
   }
   Printf("Info: number of variations = %d, central value: %d, %d global %d", nvariations, centralrec, centralgen, centralvar);
   TString subdir[nvariations];
   TString subdireje[nvariations];
   Int_t index = 0;
   for(Int_t ir = 0; ir < variationrec; ir++){
   	   if(skipVariationRec[ir]) continue;
      for(Int_t ig = 0; ig < variationgen; ig++){
      	  if(skipVariationGen[ig]) continue;
      	 Int_t indexOfAll = ir*variationgen + ig;
      	 
      	 subdir[index] = Form("%02dresp_pt%.0f_%.0for%.0f_%.0f_ptT%.0f_%.0for%.0f_%.0f_m%.0f_%.0for%.0f_%.0f_mT%.0f_%.0f", indexOfAll, pt_min[ir], pt_max[ir], pt_mineje[ir], pt_maxeje[ir], pt_minT[ig], pt_maxT[ig], pt_minTeje[ig], pt_maxTeje[ig], m_min[ir], m_max[ir], m_mineje[ir], m_maxeje[ir], m_minT[ig], m_maxT[ig]);
      	 //subdireje[index] = Form("%02dresp_pt%.0f_ptT%.0f", index, pt_mineje[ir], pt_minTeje[ig]);
      	 if(dryrun) Printf("ir = %d, ig = %d (index = %d) -> %s", ir, ig, index, subdir[index].Data());
      	 
      	 index++;
      }
   }
   
   // here other inputs that do not change often
   Int_t colTypefixed = -1;   // in the response this setting is not used because the ranges are set from outside
   Double_t binWidthPt = 10.;
   Int_t skipBins = 0;
   TString respFileNameBase="response";  if(!fillmiss) respFileNameBase="responsenoMiss";// the base name of output file containing the response. Appended the corresponding tag[]
   
   // variation of the pT bin width
   const Int_t widthvar = 2;
   Double_t binWidthPtvar[widthvar] = {5., 10.};
   TString subdirWidth[widthvar];
   for(Int_t iw = 0; iw < widthvar; iw++){
   	   subdirWidth[iw] = Form("PtBinW%.0f", binWidthPtvar[iw]);
   
   }
   
      
   // settings unfolding
   
   // input data files, for the pPb analysis we have two, one for MB and another for EJE triggers. They contain the mass vs pT spectra
   const Int_t ntr = 2;
   TString dataIn[ntr] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/AnalysisResultsMB.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/AnalysisResultsEJE.root"};
   // position 0 is the MB, position 1 is EJE
   TString dataLowEffIn[ntr] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train1289-1292/AnalysisResultsMB.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train1289-1292/AnalysisResultsEJE.root"};
   
   Int_t itermin=0;
   Int_t itermax=7;
   Int_t iterdef=3;
   Int_t varyPrior = 0;
   Double_t R = 0.4;
   // BkgType: 0 = deriv, 1 = const sub, 2 = no bkg sub
   // it comes from the order the data input is read is unfoldData.C
   Int_t nbkg = 2;
   Int_t defaultBkgType = 2;
   Int_t BkgType[nresp] = {defaultBkgType, 1, defaultBkgType, defaultBkgType};
   // ColType: 1 = MB, 2 = EJE
   Int_t ColType[ntr] = {1, 2}; // position 0 is the MB, position 1 is EJE
   
   // variable bin width
   // index of the directory corresponding to the variable bin width. It's run in a specific loop
   Int_t idvarbinw = 5;
   TString fileWithUnfRanges[ntr] = {"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrMB.root", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/DefineUnfRange/VarBin/MassVsPtVarWBkgNoSubTrJ1.root"};
   // particle level ranges (to be, seen if it's better to set them othewise)
   Double_t massbw = 2, massminT = 0, massmaxT = 40;
   Double_t ptbw = 20, ptminT = 20., ptmaxT = 140.;
   
   
   // for these two, that then become the combined one below, all the variations are performed. The other paths are considered as variation themselves and the default of the other settings are used
   TString defaultRespPaths[ntr] = {wdirresp[0], wdirresp[nresp]}; //{"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResulchntspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo", "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/EJE/redo"};
   
   TString fileunfoutbase = "UnfoldedDistributionsPrior";  // the base name of output file containing the unfolded mass vs pT. Appended the Prior variation and the corresponding tag[]
      
   // combine triggers
   // working directories. Subdirectories for the variations requested are created inside these directories
   TString wdirunf[nresp] = {
   	   //"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MBEJE/RespRhoWeighted",
   	   //"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MBEJE/RespRhoWeighted",
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingNoBkgSub/MBEJE", 
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldDet1134/MBEJE/ConstSub",       
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbOvlExclu/MBEJE",
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmb96TrEff/MBEJE",
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldingVarPrior/MBEJE",
   	   "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldVarBinW/MBEJE/"
   };
      //,
   //"
   
   TString wdirunfLowTrEffData = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train1289-1292/UnfoldEmb96TrEff/MBEJE";

   // make the directory for the unfolding results   
   // - response per each trigger
   for(Int_t idir = 0; idir < nunf; idir++){
   	   gSystem->mkdir(wdirresp[idir], kTRUE);
   }
   // - unfolding low data trk eff
   for(Int_t idir = 0; idir < ntr; idir++){
   	   gSystem->mkdir(wdirrespLowTrEff[idir], kTRUE);
   }
   // - combine triggers
   for(Int_t idir = 0; idir < nresp; idir++){
   	   gSystem->mkdir(wdirunf[idir], kTRUE);
   }
   gSystem->mkdir(wdirunfLowTrEffData, kTRUE);
   
   gSystem->ChangeDirectory(basewdircmp);

   
   //check that these stay the same, mostly the number (that is the index or the array wdirunf corresponding to the defaultUnfPaths)
   Int_t defaultUnfNumb = 0;
   TString defaultUnfPaths = wdirunf[defaultUnfNumb]; // it's: "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MBEJE";
   Int_t binSwitchToEJE = -1; // it's calculated when combining the triggers in the range variation, see if needed to reset it in all.
   
   Int_t arrayiter[ntr] = {iterdef, iterdef};
   TString namescombo[ntr] = {Form("MB%s", suff.Data()), Form("EJE%s", suff.Data())};
   TString fileunfsumbase = "MassUnfSum";  // the base name of output file containing the sum of the unfolded mass vs pT for the two triggers. Appended the identifying of the summed files, given in namescombo[]
   //filename output
   TString outfileunfmergename = Form("%s", fileunfsumbase.Data());
   for(Int_t i = 0; i<ntr; i++) outfileunfmergename+=namescombo[i];
   outfileunfmergename+=".root";
   TString basehnname = "hMUnf__Iter3";
   // CompareResults can handle the name in this format by adding the iteration number in between the __ .
   // the name will be Form("%s%d%s%d", basehnnamebeforePt.Data(), ipt, basehnnameafterPt.Data(), arrayiter[itr]);
   
   // "hUnfMpj_Itr3_ptb" is the one taken from the projecton of the 2D spectrum sum of MB and EJE, so it's wrong!!
      
   // settings comparisons
   
   //"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/testSystematic";
   
   TString wdircmpBkg = "CmpEmbBkgSub"; // compare bkg subtraction, default response ranges
   TString inputcmpBkg[nbkg];
   TString inputhnamecmpBkg[nbkg] = {basehnname, basehnname};
   TString legcmpBkg[nbkg] = {"NoSub", "DetOnly"};
   Int_t offsetcmpBkg[nbkg] = {0, 0};
 
   // save the path for later
   for(Int_t ibg = 0; ibg<nbkg; ibg++){
   	   inputcmpBkg[ibg] = wdirunf[ibg]; // this takes the first and second directory, be aware!!
   	   inputcmpBkg[ibg] += "/";
   	   inputcmpBkg[ibg] += subdir[centralvar];
   	   inputcmpBkg[ibg] += "/";
   	   inputcmpBkg[ibg] += outfileunfmergename;
   	   
   	   Printf("************ %s , subdir %s", inputcmpBkg[ibg].Data(), subdir[centralvar].Data());
      	              
   }
   
   Int_t novl = 2;
   TString wdircmpOvl = "CmpEmbNoOvl"; // compare overlap single-track - jet exclusion, default
   TString inputcmpOvl[novl];
   TString inputhnamecmpOvl[novl] = {basehnname, basehnname};
   TString legcmpOvl[novl] = {"Ovl", "NoOvl"};
   Int_t offsetcmpOvl[novl] = {0, 0};
 
   // save the path for later
   inputcmpOvl[0] = defaultUnfPaths;
   inputcmpOvl[1] = wdirunf[2];
   
   for(Int_t ibg = 0; ibg<novl; ibg++){
   	   
   	   inputcmpOvl[ibg] += "/";
   	   inputcmpOvl[ibg] += subdir[centralvar];
   	   inputcmpOvl[ibg] += "/";
   	   inputcmpOvl[ibg] += outfileunfmergename;
      	              
   }
   
   TString wdirUnfSystRanges = "CmpRespMoreRanges"; // compare different unfolding ranges, use Deriv bkg sub
   Int_t   nranges = nvariations; if(variablebinning) nranges = nvariations+1; // add the variation with fixed (nvariations) and variable binning (1)
   TString inputUnfSystRanges[nranges];
   TString inputhnameUnfSystRanges[nranges] ;
   TString legUnfSystRanges[nranges] ;
   Int_t offsetUnfSystRanges[nranges] ;
   TString kineeffMBFilenameSystRanges[nranges] ;
   TString kineeffJEFilenameSystRanges[nranges] ;
   Int_t   idhistoRangeVarW = 1; // variable bin width histogram standard
   ofstream longLegendSysRange;
   longLegendSysRange.open("LongLegendSysRange.txt");
   
   Int_t indexg = 0;
   for(Int_t ir = 0; ir < variationrec; ir++){
   	   if(skipVariationRec[ir]) continue;
   	   for(Int_t ig = 0; ig < variationgen; ig++){
   	   	   if(skipVariationGen[ig]) continue;
   	  	   //Int_t index = ir*variationgen + ig;
   	  	   Printf("Variation rec %d, gen %d, index %d", ir, ig, indexg);
   	  	   inputhnameUnfSystRanges[indexg] = basehnname;
   	  	  
   	  	   const char *longtext = Form("MBrec_%.0f_%.0f_%.0f_%.0f_gen_%.0f_%.0f_%.0f_%.0fEJErec_%.0f_%.0f_%.0f_%.0f_gen_%.0f_%.0f_%.0f_%.0f", pt_min[ir], pt_max[ir], m_min[ir], m_max[ir], pt_minT[ig], pt_maxT[ig], m_minT[ig], m_maxT[ig], pt_mineje[ir], pt_maxeje[ir], m_mineje[ir], m_maxeje[ir], pt_minTeje[ig], pt_maxTeje[ig], m_minT[ig], m_maxT[ig]);
   	  	   longLegendSysRange<<"Rec "<< ir << " Gen "<< ig << " " << longtext << endl;
   	  	   
   	  	   legUnfSystRanges[indexg] = Form("Rec%dGen%d", ir, ig);
   	  	   
   	  	   offsetUnfSystRanges[indexg] = 0;
   	  	   
   	  	   //preparation for kine eff
   	  	   kineeffMBFilenameSystRanges[indexg] = Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", pt_min[ir], pt_max[ir], m_min[ir], m_max[ir]);
   	  	   
   	  	   kineeffJEFilenameSystRanges[indexg] = Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", pt_mineje[ir], pt_maxeje[ir], m_mineje[ir], m_maxeje[ir]);
   	  	   indexg++;
   	   }

   }
   
   
   // last bin is for the variable binning
   // use the index incremented in the loop before - indexg, don't move!
   if(variablebinning){
   	   inputhnameUnfSystRanges[indexg] = basehnname;
   	   legUnfSystRanges[indexg] = "Variable";
   	   offsetUnfSystRanges[indexg] = 0;
   }
   
   TString wdirUnfSystIter = "CmpNIterations";   // compare different iterations
   Int_t iterminsyst = 2,  itermaxsyst = 6;
   const Int_t nitervar =  itermaxsyst - iterminsyst;
   // directory combination of different iterations
   TString iterdir[nitervar];
   // for the comparison
   TString inputUnfSystIter[nitervar];
   TString inputhnameUnfSystIter[nitervar];
   TString legUnfSystIter[nitervar];
   Int_t offsetUnfSystIter[nitervar];
   TString dirunfIterSyst = defaultUnfPaths; //wdirunf[2];//
   for(Int_t index = 0; index < nitervar; index++){
      iterdir[index] = Form("Iter%d", iterminsyst+index);
      inputhnameUnfSystIter[index] = Form("hMUnf__Iter%d", iterminsyst+index);
      legUnfSystIter[index] = Form("Iter_%.d", iterminsyst+index);
      offsetUnfSystIter[index] = 0;
      //it's always the same path
      inputUnfSystIter[index] = dirunfIterSyst; 
      inputUnfSystIter[index] += "/";
      inputUnfSystIter[index] += iterdir[index];
      inputUnfSystIter[index] += "/";
      inputUnfSystIter[index] += outfileunfmergename;
   }

   TString wdirUnfSystBinW = "CmpBinW";
   TString inputcmpBinW[widthvar];
   TString inputhnameBinW[widthvar];
   TString legUnfSystBinW[widthvar];
   Int_t offsetUnfSystBinW[widthvar];
   // save the path for later
   for(Int_t iw = 0; iw<widthvar; iw++){
   	   inputhnameBinW[iw] = basehnname;
   	   legUnfSystBinW[iw] = Form("BinW%.0f", binWidthPtvar[iw]);
   	   offsetUnfSystBinW[iw] = 0;
   	   inputcmpBinW[iw] = defaultUnfPaths;
   	   inputcmpBinW[iw] += "/";
   	   inputcmpBinW[iw] += subdirWidth[iw];
   	   inputcmpBinW[iw] += "/";
   	   inputcmpBinW[iw] += outfileunfmergename;
      	              
   }
   
   Int_t nTrkE = 2;
   Int_t idTrE = 3;
   TString wdircmpTrkE = "CmpTrackingEff"; // compare overlap single-track - jet exclusion, default
   
   TString inputcmpTrkE[nTrkE];
   TString inputhnamecmpTrkE[nTrkE] = {basehnname, basehnname};
   TString legcmpTrkE[nTrkE] = {"Std", "96PercEff"};
   Int_t offsetcmpTrkE[nTrkE] = {0, 0};
 
   // this is for the reduced traking efficiency in data as well
   TString wdircmpTrkEv2 = "CmpTrackingEffv2";
   TString inputcmpTrkEv2[nTrkE];
   TString legcmpTrkEv2[nTrkE] = {"Std", "96PercEffv2"};
   // save the path for later
   inputcmpTrkE[0] = defaultUnfPaths;
   inputcmpTrkE[1] = wdirunf[idTrE];
   inputcmpTrkEv2[0] = defaultUnfPaths;
   inputcmpTrkEv2[1] = wdirunfLowTrEffData;
   
   Printf("Inputs Trk eff systematics:");
   for(Int_t ibg = 0; ibg<nTrkE; ibg++){
   	   
   	   inputcmpTrkE[ibg] += "/";
   	   inputcmpTrkE[ibg] += subdir[centralvar];
   	   inputcmpTrkE[ibg] += "/";
   	   inputcmpTrkE[ibg] += outfileunfmergename;
   	   Printf("%d. %s", ibg, inputcmpTrkE[ibg].Data());
   	   
   	   inputcmpTrkEv2[ibg] += "/";
   	   inputcmpTrkEv2[ibg] += subdir[centralvar];
   	   inputcmpTrkEv2[ibg] +=  "/";
   	   inputcmpTrkEv2[ibg] +=  outfileunfmergename;
   	   Printf("%d. %s", ibg, inputcmpTrkEv2[ibg].Data());
   }
   
   Int_t nvarprior = 2;
   Int_t idPrior = 4;
   Double_t smearingfactor = 2; //it's 2% in the smearing of the tree
   TString wdircmpPrior = "CmpVarPrior";
   TString inputcmpPrior[nvarprior];
   TString inputhnamecmpPrior[nvarprior] = {basehnname, basehnname};
   TString legcmpPrior[nvarprior] = {"Std", "2PercSmear"};
   Int_t offsetcmpPrior[nvarprior] = {0, 0};
   Int_t idhistoRangePrior = 4; //default range and bin width histogram
   
   inputcmpPrior[0] = defaultUnfPaths;
   inputcmpPrior[1] = wdirunf[idPrior];
   for(Int_t ibg = 0; ibg<nvarprior; ibg++){
   
   	   inputcmpPrior[ibg] += "/";
   	   inputcmpPrior[ibg] += subdir[centralvar];
   	   inputcmpPrior[ibg] += "/";
   	   inputcmpPrior[ibg] += outfileunfmergename;
   	   Printf("%d. %s", ibg, inputcmpTrkE[ibg].Data());
   }
   
   
   TString filetotsystname = "TotalSystematicUnc.root", basenamehsystot = "hSystTot";
   
   // ++++++++++++ Start of the loops
   
   // define response matrices and unfold
   TH1::AddDirectory(kFALSE);
   
   for(Int_t imatr = 0; imatr < nunf; imatr++){
      Printf("\n-----------Main loop %d", imatr);
      
      
      Int_t jinput = imatr;
      if (imatr >= nresp) {
      	 Int_t k = imatr/nresp;
      	 jinput = imatr - k*nresp;
      }
      
      Int_t jdata = 0;
      if(wdirresp[imatr].Contains("Const")) jdata = 1;
      
      Int_t jtrig = 0;
      if(wdirresp[imatr].Contains("EJE")) jtrig = 1;
      if(variablebinning){
      	  if(imatr == idvarbinw || imatr == idvarbinw + nresp) {
      	  	  
      	  	  //unfolding for variable bin width
      	  	  
      	  	  Printf("\n-----------Loop  creation respose with variable binning %d", jtrig);
      	  	  
      	  	  gSystem->ChangeDirectory(wdirresp[imatr]);
      	  	  //gSystem->mkdir(subdir[centralvar]);
      	  	  //gSystem->ChangeDirectory(subdir[centralvar]);
      	  	  Printf("Response will be saved in directory %s", gSystem->pwd());
      	  	  
      	  	  if(!dryrun) {
      	  	  	  TString outfilerespname = Form("%s%s.root", respFileNameBase.Data(), tag[idvarbinw].Data());
      	  	  	  Printf("Outfilename = %s, here %s", outfilerespname.Data(), gSystem->pwd());
      	  	  	  if(domatrix){
      	  	  	  	  TFile testfile(outfilerespname);
      	  	  	  	  if(testfile.IsOpen()){
      	  	  	  	  	  Printf("Response %s already present, don't run again", outfilerespname.Data());
      	  	  	  	  	  
      	  	  	  	  } else {
      	  	  	  	  	  Printf("file with ranges = %s, input %s", fileWithUnfRanges[jtrig].Data(), strIn[jinput].Data());
      	  	  	  	  	  CreateRooUnfoldResponseVarWidthFromFile(fileWithUnfRanges[jtrig], massbw, massminT, massmaxT, ptbw, ptminT, ptmaxT, strIn[jinput], Form("%s%s", fillmiss ? "" : "noMiss", tag[idvarbinw].Data()), idhistoRangeVarW, fillmiss);
      	  	  	  	  }
      	  	  	  }
      	  	  	  // implement this...
      	  	  	  
      	  	  	  //if(kineeffinrangesyst){
      	  	  	  // 	   TFile *feffMB = new TFile(kineeffMBFilenameSystRanges[indexg]);
      	  	  	  // 	   TFile *feffJE = new TFile(kineeffJEFilenameSystRanges[indexg]);
      	  	  	  // 	   if(reckineff || !feffMB->IsOpen()){
      	  	  	  // 	   	   // calculate Kinematic efficiencies
      	  	  	  // 	   	   CalculateKineEff(pt_min[ir], pt_max[ir], m_min[ir], m_max[ir], strIn[defaultUnfNumb], strL[defaultUnfNumb], "");
      	  	  	  // 	   }
      	  	  	  // 	   if(reckineff || !feffJE->IsOpen()){
      	  	  	  // 	   	   // calculate Kinematic efficiencies
      	  	  	  // 	   	   CalculateKineEff(pt_mineje[ir], pt_maxeje[ir], m_mineje[ir], m_maxeje[ir], strIn[defaultUnfNumb], strL[defaultUnfNumb], "");
      	  	  	  // 	   }
      	  	  	  //} //calculate kine eff if requested
      	  	  	  
      	  	  	  if (dounfold) {
      	  	  	  	  TString outfileunfname = Form("%s%d%s.root", fileunfoutbase.Data(), varyPrior, tag[idvarbinw].Data());
      	  	  	  	  TFile testfile(outfileunfname);
      	  	  	  	  if(!forcedounfold && testfile.IsOpen()){
      	  	  	  	  	  Printf("Unfolded output %s already present, don't run again", outfileunfname.Data() );
      	  	  	  	  	  
      	  	  	  	  }else {
      	  	  	  	  	  Printf("jinput = %d, jdata (bkg type) = %d, jtrig (trigger type) = %d, input data file = %s", jinput, jdata, jtrig, fileWithUnfRanges[jtrig].Data());
      	  	  	  	  	  unfold(outfilerespname, fileWithUnfRanges[jtrig], iterdef, itermin, itermax, BkgType[jinput], varyPrior, R, ColType[jtrig], tag[jinput].Data(), 1, idhistoRangeVarW);
      	  	  	  	  }
      	  	  	  }
      	  	  }
      	  	  // don't do the other variations below
      	  	  continue;
      	  }
      }
      
      //create response and unfold with prior variation
      Printf("imatr = %d, iprior = %d", imatr, idPrior);
      if(imatr == idPrior || imatr == idPrior + nresp){
      	      //prior variation
      	  	  
      	  	  Printf("\n-----------Loop  creation respose with varied priors %d", jtrig);
      	  	  
      	  	  gSystem->ChangeDirectory(wdirresp[imatr]);
      	  	  gSystem->mkdir(subdir[centralvar]);
      	  	  gSystem->ChangeDirectory(subdir[centralvar]);
      	  	  Printf("Response will be saved in directory %s", gSystem->pwd());
      	  	  
      	  	  if(!dryrun) {
      	  	  	  TString outfilerespname = Form("%s%s.root", respFileNameBase.Data(), tag[idPrior].Data());
      	  	  	  if(domatrix){
      	  	  	  	  TFile testfile(outfilerespname);
      	  	  	  	  if(testfile.IsOpen()){
      	  	  	  	  	  Printf("Response %s already present, don't run again", outfilerespname.Data());
      	  	  	  	  	  
      	  	  	  	  } else {
      	  	  	  	  	  Printf("file with ranges = %s, input %s", fileWithUnfRanges[jtrig].Data(), strIn[jinput].Data());
      	  	  	  	  	  SmearedResp(fileWithUnfRanges[jtrig], massbw, massminT, massmaxT, ptbw, ptminT, ptmaxT, strIn[defaultUnfNumb], Form("%s%s", fillmiss ? "" : "noMiss", tag[idPrior].Data()), idhistoRangePrior, fillmiss, smearingfactor, strIn[jinput]);
      	  	  	  	  }
      	  	  	  }
      	  	  	  // implement this...
      	  	  	  
      	  	  	  
      	  	  	  if (dounfold) {
      	  	  	  	  TString outfileunfname = Form("%s%d%s.root", fileunfoutbase.Data(), varyPrior, tag[idPrior].Data());
      	  	  	  	  TFile testfile(outfileunfname);
      	  	  	  	  if(!forcedounfold && testfile.IsOpen()){
      	  	  	  	  	  Printf("Unfolded output %s already present, don't run again", outfileunfname.Data() );
      	  	  	  	  	  
      	  	  	  	  }else {
      	  	  	  	  	  Printf("jinput = %d, jdata (bkg type) = %d, jtrig (trigger type) = %d, input data file = %s", jinput, jdata, jtrig, fileWithUnfRanges[jtrig].Data());
      	  	  	  	  	  unfold(outfilerespname, fileWithUnfRanges[jtrig], iterdef, itermin, itermax, BkgType[jinput], varyPrior, R, ColType[jtrig], tag[jinput].Data(), 1, idhistoRangePrior);
      	  	  	  	  }
      	  	  	  }
      	  	  }
      	  	  // don't do the other variations below
      	  	  continue;
      
      }
      
      Int_t countfalse = 0;
      //check the path: if not "default" don't run the variation or not (do it only if it's the default output)
      
      for(Int_t itr = 0; itr<ntr; itr++){
      	  if(wdirresp[imatr] != defaultRespPaths[itr]) countfalse++;
      	  //Printf("%s =? %s, %d, %d", wdirresp[imatr].Data(), defaultRespPaths[itr].Data(), ir, ig);
      }
      indexg = -1;
      for(Int_t ir = 0; ir < variationrec; ir++){
      	  if(skipVariationRec[ir]) continue;
      	 for(Int_t ig = 0; ig < variationgen; ig++){
      	    if(skipVariationGen[ig]) continue;
      	    Bool_t skip = kFALSE;
      	     
      	    if((countfalse == 2) && ((ig!=centralgen) || (ir!=centralrec))) {
      	       	  skip = kTRUE;
      	       	  Printf("skipping false %d %d, %d", countfalse, ir, ig);
      	       	  break;
      	    }
      	    
      	    
      	    if(skip) continue;
      	    else Printf("Variation range from default setting in default ");
      	    indexg++;
      	    gSystem->ChangeDirectory(wdirresp[imatr]);
      	    //Int_t index = ir*variationgen + ig;
      	    gSystem->mkdir(subdir[indexg]);
      	    gSystem->ChangeDirectory(subdir[indexg]);
      	    
      	    Printf("%d (ir %d, ig %d) ) Running create response pT min = %.0f, pT max = %.0f, M min = %.0f, M max = %.0f, pTtrue min = %.0f, pTtrue max = %.0f, Mtrue min = %.0f, Mtrue max = %.0f in directory %s", indexg, ir, ig, pt_min[ir], pt_max[ir], m_min[ir], m_max[ir],  pt_minT[ig],  pt_maxT[ig], m_minT[ig], m_maxT[ig], gSystem->pwd());
      	    
      	    if(!dryrun) {
      	       TString outfilerespname = Form("%s%s.root", respFileNameBase.Data(), tag[jinput].Data());
      	       if(domatrix){
      	       	  TFile testfile(outfilerespname);
      	       	  if(testfile.IsOpen()){
      	       	     Printf("Response %s already present, don't run again", outfilerespname.Data());
      	       	     
      	       	  } else {
      	       	  	  if(wdirresp[imatr].Contains("EJE")){
      	       	  	  	  CreateRooUnfoldResponse(strIn[jinput], strL[jinput], Form("%s%s", fillmiss ? "" : "noMiss", tag[jinput].Data()), binWidthPt, skipBins, pt_mineje[ir], m_mineje[ir],  pt_minTeje[ig], m_minT[ig], pt_maxeje[ir], m_maxeje[ir],  pt_maxTeje[ig], m_maxT[ig], colTypefixed, fillmiss);
      	       	  	  } else {
      	       	  	  	  CreateRooUnfoldResponse(strIn[jinput], strL[jinput], Form("%s%s", fillmiss ? "" : "noMiss", tag[jinput].Data()), binWidthPt, skipBins, pt_min[ir], m_min[ir],  pt_minT[ig], m_minT[ig], pt_max[ir], m_max[ir],  pt_maxT[ig], m_maxT[ig], colTypefixed, fillmiss);
      	       	      }
      	       	  }
      	       }
      	       if (dounfold) {
      	       	  TString outfileunfname = Form("%s%d%s.root", fileunfoutbase.Data(), varyPrior, tag[jinput].Data());
      	       	  TFile testfile(outfileunfname);
      	       	  if(!forcedounfold && testfile.IsOpen()){
      	       	     Printf("Unfolded output %s already present, don't run again", outfileunfname.Data() );
      	       	     
      	       	  }else {
      	       	     Printf("jinput = %d, jdata (bkg type) = %d, jtrig (trigger type) = %d, input data file = %s", jinput, jdata, jtrig, dataIn[jtrig].Data());
      	       	     unfold(outfilerespname, dataIn[jtrig], iterdef, itermin, itermax, BkgType[jinput], varyPrior, R, ColType[jtrig], tag[jinput].Data());
      	       	  }
      	       } // do unfold
      	       
      	       if(kineeffinrangesyst){
      	       	   TFile *feffMB = new TFile(kineeffMBFilenameSystRanges[indexg]);
      	       	   TFile *feffJE = new TFile(kineeffJEFilenameSystRanges[indexg]);
      	       	   if(reckineff || !feffMB->IsOpen()){
      	       	   	   // calculate Kinematic efficiencies
      	       	   	   CalculateKineEff(pt_min[ir], pt_max[ir], m_min[ir], m_max[ir], strIn[defaultUnfNumb], strL[defaultUnfNumb], "");
      	       	   }
      	       	   if(reckineff || !feffJE->IsOpen()){
      	       	   	   // calculate Kinematic efficiencies
      	       	   	   CalculateKineEff(pt_mineje[ir], pt_maxeje[ir], m_mineje[ir], m_maxeje[ir], strIn[defaultUnfNumb], strL[defaultUnfNumb], "");
      	       	   }
   	   		   } //calculate kine eff if requested
   	   		   
      	    } // dryrun
      	 } // variation gen
      } // variation rec    
      
      
   } // main loop to create response and unfold
   
   for(Int_t itr = 0 ; itr < ntr; itr++){ //loop on the deafult path for the other variations (no range)
   	   
   	   for(Int_t iw = 0; iw < widthvar; iw++){
   	   	   
   	   	   gSystem->ChangeDirectory(defaultRespPaths[itr]);
   	   	   gSystem->mkdir(subdirWidth[iw]);
   	   	   gSystem->ChangeDirectory(subdirWidth[iw]);
   	   	   Printf("\nVariation width from default setting in default path %s", gSystem->pwd());
   	   	   //Printf("%d Running create response with pT bin W = %.0f (pT min = %.0f, M min = %.0f, pTtrue min = %.0f, Mtrue min = %.0f) in directory %s", index, binWidthPtvar[iw], pt_min[centralrec], m_min[centralrec],  pt_minT[centralgen], m_minT[centralgen], gSystem->pwd());
   	   	   
   	   	   if(!dryrun) {
   	   	   	   TString outfilerespname = Form("%s%s.root", respFileNameBase.Data(), tag[0].Data());
   	   	   	   if(domatrix){
   	   	   	   	   TFile testfile(outfilerespname);
   	   	   	   	   if(testfile.IsOpen()){
   	   	   	   	   	   Printf("Response %s already present, don't run again", outfilerespname.Data());
   	   	   	   	   	   
   	   	   	   	   } else {
   	   	   	   	   	   if(wdirresp[itr].Contains("EJE")){
   	   	   	   	   	   	   CreateRooUnfoldResponse(strIn[0], strL[0], Form("%s%s", fillmiss ? "" : "noMiss", tag[0].Data()), binWidthPtvar[iw], skipBins, pt_mineje[centralrec], m_min[centralrec],  pt_minTeje[centralgen], m_minT[centralgen], pt_maxeje[centralrec], m_max[centralrec],  pt_maxTeje[centralgen], m_maxT[centralgen], colTypefixed, fillmiss);
   	   	   	   	   	   } else {
   	   	   	   	   	   	   CreateRooUnfoldResponse(strIn[0], strL[0], Form("%s%s", fillmiss ? "" : "noMiss", tag[0].Data()), binWidthPtvar[iw], skipBins, pt_min[centralrec], m_min[centralrec],  pt_minT[centralgen], m_minT[centralgen], pt_max[centralrec], m_max[centralrec],  pt_maxT[centralgen], m_maxT[centralgen], colTypefixed, fillmiss);
   	   	   	   	   	   }
   	   	   	   	   }
   	   	   	   }
   	   	   	   if (dounfold) {
   	   	   	   	   TString outfileunfname = Form("%s%d%s.root", fileunfoutbase.Data(), varyPrior, tag[0].Data());
   	   	   	   	   TFile testfile(outfileunfname);
   	   	   	   	   if(!forcedounfold && testfile.IsOpen()){
   	   	   	   	   	   Printf("Unfolded output %s already present, don't run again", outfileunfname.Data() );
   	   	   	   	   	   
   	   	   	   	   }else {
   	   	   	   	   	   Printf("itr = %d, input data file = %s", itr, dataIn[itr].Data());
   	   	   	   	   	   unfold(outfilerespname, dataIn[itr], iterdef, itermin, itermax, defaultBkgType, varyPrior, R, ColType[itr], tag[0].Data());
   	   	   	   	   }
      	       } // do unfold
      	   } // dry run
       } //width variation loop
       
       
       // tracking efficiency systematics using the data with lower tracking efficiency. The response is already ready, here do unfolding
       gSystem->ChangeDirectory(wdirrespLowTrEff[itr]);
       gSystem->mkdir(subdir[centralvar]);
       gSystem->ChangeDirectory(subdir[centralvar]);
       Printf("\nunfolding data with lower tracking efficiency in %s with input resp %s/%s%s.root and input data %s", gSystem->pwd(), wdirresp[(itr+1)*idTrE].Data(), respFileNameBase.Data(), tag[idTrE].Data(), dataLowEffIn[itr].Data());
       
       if(!dryrun){
       	   TString outfilerespname = Form("%s/%s%s.root", wdirresp[(itr+1)*idTrE].Data(), respFileNameBase.Data(), tag[idTrE].Data());
       	   if (dounfold) {
       	   	   TString outfileunfname = Form("%s%d%s.root", fileunfoutbase.Data(), varyPrior, tag[idTrE].Data());
       	   	   TFile testfile(outfileunfname);
       	   	   if(!forcedounfold && testfile.IsOpen()){
       	   	   	   Printf("Unfolded output %s already present, don't run again", outfileunfname.Data() );
       	   	   	   
       	   	   }else {
       	   	   	   Printf("itr = %d, input data file = %s", itr, dataLowEffIn[itr].Data());
       	   	   	   unfold(outfilerespname, dataLowEffIn[itr], iterdef, itermin, itermax, defaultBkgType, varyPrior, R, ColType[itr], tag[idTrE].Data());
       	   	   }
       	   }
       }
   }
   
   // sum (or better, stick together) the two triggered samples created in the previous loop
   
   for(Int_t idata = 0; idata <nresp; idata++){
      Printf("\n-----------Loop combine triggers %d", idata);
      
      
      if(variablebinning && idata == idvarbinw){
      	  Printf("-> variable bin width");
      	  
      	  gSystem->ChangeDirectory(wdirunf[idvarbinw]); // go to the base directory for this (idata) data sample
      	  
      	  TString filesvarbinw[ntr] = {Form("%s/%s%d%s.root", wdirresp[idvarbinw].Data(), fileunfoutbase.Data(), varyPrior, tag[idvarbinw].Data()),  Form("%s/%s%d%s.root", wdirresp[idvarbinw+nresp].Data(),  fileunfoutbase.Data(), varyPrior, tag[idvarbinw].Data()) };
      	  Printf("Combining files: %s and %s", filesvarbinw[0].Data(), filesvarbinw[1].Data());
      	  
      	  
      	  
      	  if (dounfold) {
      	  	  
      	  	  Printf("%d) Running combine trigger for variable bin width in directory  %s", indexg, gSystem->pwd());
      	  	  
      	  	  if(!dryrun) {
      	  	  	  
      	  	  	  TFile testfile(outfileunfmergename);
      	  	  	  if(!forcedounfold && testfile.IsOpen()){
      	  	  	  	  Printf("Unfolded output %s/%s already present, don't run again", gSystem->pwd(), outfileunfmergename.Data() );
      	  	  	  }else {
      	  	  	  	  binSwitchToEJE = CombineTriggersCorrected(ntr, filesvarbinw, arrayiter, namescombo);
      	  	  	  }
      	  	  }
      	  }
      	  inputUnfSystRanges[nvariations] =  wdirunf[idvarbinw];
      	  inputUnfSystRanges[nvariations] += "/";
      	  inputUnfSystRanges[nvariations] += outfileunfmergename;
      	  
      	  Printf("Save [%d] = %s (this is the varible bin)", nvariations, inputUnfSystRanges[nvariations].Data());
      	  
      	  
      } else {
      	  
      	  Printf("-> Default and range variations");
      	  indexg = -1;
      	  for(Int_t ir = 0; ir < variationrec; ir++){
      	  	  if(skipVariationRec[ir]) continue;
      	  	  for(Int_t ig = 0; ig < variationgen; ig++){
      	  	  	  if(skipVariationGen[ig]) continue;
      	  	  	  
      	  	  	  Bool_t skip = kFALSE;
      	  	  	  //check the path: if not "default" don't run the variation or not (do it only if it's the default output)
      	  	  	  if(defaultUnfPaths != wdirunf[idata]) {
      	  	  	  	  
      	  	  	  	  if((ig!=centralgen) || (ir!=centralrec)) {
      	  	  	  	  	  skip = kTRUE;
      	  	  	  	  	  Printf("skipping %d, %d", ir, ig);
      	  	  	  	  	  break;
      	  	  	  	  }
      	  	  	  }
      	  	  	  
      	  	  	  if(skip) continue;
      	  	  	  else Printf("Variation from default setting in default ");
      	  	  	  indexg++;
      	  	  	  //Int_t index = ir*variationgen + ig;
      	  	  	  
      	  	  	  gSystem->ChangeDirectory(wdirunf[idata]); // go to the base directory for this (idata) data sample
      	  	  	  
      	  	  	  TString files[ntr] = {Form("%s/%s/%s%d%s.root", wdirresp[idata].Data(), subdir[indexg].Data(), fileunfoutbase.Data(), varyPrior, tag[idata].Data()),  Form("%s/%s/%s%d%s.root", wdirresp[idata+nresp].Data(), subdir[indexg].Data(), fileunfoutbase.Data(), varyPrior, tag[idata].Data()) };
      	  	  	  Printf("Combining files: %s and %s", files[0].Data(), files[1].Data());
      	  	  	  
      	  	  	  
      	  	  	  
      	  	  	  if (dounfold) {
      	  	  	  	  
      	  	  	  	  // create and/or go to the specific directory of the current variation 
      	  	  	  	  gSystem->mkdir(subdir[indexg]);
      	  	  	  	  gSystem->ChangeDirectory(subdir[indexg]);
      	  	  	  	  
      	  	  	  	  Printf("%d ) ----->>>>>>>> Subdir = %s Running combine triggers pT min = %.0f, pT max = %.0f, pT min eje = %.0f, pT max eje = %.0f, M min = %.0f, M max = %.0f, M min eje = %.0f, M max eje = %.0f, pTtrue min = %.0f, pTtrue max = %.0f, pTtrue min eje = %.0f, pTtrue max eje = %.0f, Mtrue min = %.0f, Mtrue max = %.0f in directory  %s", indexg, subdir[indexg].Data(), pt_min[ir], pt_max[ir], pt_mineje[ir], pt_maxeje[ir], m_min[ir], m_max[ir], m_mineje[ir], m_maxeje[ir],  pt_minT[ig],  pt_maxT[ig], m_minT[ig], m_maxT[ig],  pt_minTeje[ig],  pt_maxTeje[ig], gSystem->pwd());
      	  	  	  	  
      	  	  	  	  
      	  	  	  	  if(!dryrun) {
      	  	  	  	  	  
      	  	  	  	  	  TFile testfile(outfileunfmergename);
      	  	  	  	  	  if(!forcedounfold && testfile.IsOpen()){
      	  	  	  	  	  	  Printf("Unfolded output %s already present, don't run again", outfileunfmergename.Data() );
      	  	  	  	  	  }else {
      	  	  	  	  	  	  binSwitchToEJE = CombineTriggersCorrected(ntr, files, arrayiter, namescombo);
      	  	  	  	  	  }
      	  	  	  	  }
      	  	  	  }
      	  	  	  
      	  	  	  if(idata == defaultUnfNumb){ // deriv 
      	  	  	  	  inputUnfSystRanges[indexg] =  wdirunf[idata];
      	  	  	  	  inputUnfSystRanges[indexg] += "/";
      	  	  	  	  inputUnfSystRanges[indexg] += subdir[indexg];
      	  	  	  	  inputUnfSystRanges[indexg] += "/";
      	  	  	  	  inputUnfSystRanges[indexg] += outfileunfmergename;
      	  	  	  	  
      	  	  	  	  Printf("Save [%d] = %s", indexg, inputUnfSystRanges[indexg].Data());
      	  	  	  }
      	  	  } // gen range variation
      	  }// reco range variation
      } // range variation fixed binning
      
      for(Int_t i = 0; i<nranges; i++){
      	  longLegendSysRange<<"Dir " << i << " = "<< inputUnfSystRanges[i] << endl;
      }
      longLegendSysRange.close();
      Printf("-> Default and number of iterations");
      for(Int_t iter = 0; iter < nitervar; iter++){
      	  
      	  Bool_t skip = kFALSE;
      	  //check the path: if not "default" don't run the variation or not (do it only if it's the default output)
      	  if(dirunfIterSyst != wdirunf[idata]) {
      	  	  skip = kTRUE;
      	  	  break;
      	  }
      	  
      	  if(skip) continue;
      	  else Printf("Variation from default setting in default ");
      	  gSystem->ChangeDirectory(wdirunf[idata]);
      	  
      	  
      	  TString files[ntr] = {Form("%s/%s/%s%d%s.root", wdirresp[idata].Data(), subdir[centralvar].Data(), fileunfoutbase.Data(), varyPrior, tag[idata].Data()),  Form("%s/%s/%s%d%s.root", wdirresp[idata+nresp].Data(), subdir[centralvar].Data(), fileunfoutbase.Data(), varyPrior, tag[idata].Data()) };
      	  Printf("Combining files: %s and %s", files[0].Data(), files[1].Data());
      	  
      	  
      	  if (dounfold) {
      	  	  
      	  	  // create and/or go to the specific directory of the current variation 
      	  	  gSystem->mkdir(iterdir[iter]);
      	  	  gSystem->ChangeDirectory(iterdir[iter]);
      	  	  
      	  	  Printf("%d ) Running combine triggers Iteration %d, pT min = %.0f, M min = %.0f, pTtrue min = %.0f, Mtrue min = %.0f in directory  %s", centralvar, iter, pt_min[centralrec], m_min[centralrec],  pt_minT[centralgen], m_minT[centralgen], gSystem->pwd());
      	  	  
      	      
      	  	  if(!dryrun) {
      	       	  
      	       	  TFile testfile(outfileunfmergename);
      	       	  if(!forcedounfold && testfile.IsOpen()){
      	       	  	  Printf("Unfolded output %s already present, don't run again", outfileunfmergename.Data() );
      	       	  }else {
      	       	  	  Int_t arrayitersvar[ntr] = {iterminsyst+iter, iterminsyst+iter};
      	       	  	  CombineTriggersCorrected(ntr, files, arrayitersvar, namescombo);
      	       	  }
      	      }
      	  }
      }
      
      Printf("-> Default and Pt bin width");
      for(Int_t iw = 0; iw < widthvar; iw++){
      	  
      	  Bool_t skip = kFALSE;
      	  //check the path: if not "default" don't run the variation or not (do it only if it's the default output)
      	  if(defaultUnfPaths != wdirunf[idata]) {
      	  	  skip = kTRUE;
      	  	  break;
      	  }
      	  
      	  if(skip) continue;
      	  else Printf("Variation bin W from default setting in default %s/%s", defaultUnfPaths.Data(), subdirWidth[iw].Data());
      	  gSystem->ChangeDirectory(wdirunf[idata]);
      	  gSystem->mkdir(subdirWidth[iw]);
      	  gSystem->ChangeDirectory(subdirWidth[iw]);
      	  	  
      	  TString files[ntr] = {Form("%s/%s/%s%d%s.root", wdirresp[idata].Data(), subdirWidth[iw].Data(), fileunfoutbase.Data(), varyPrior, tag[0].Data()),  Form("%s/%s/%s%d%s.root", wdirresp[idata+nresp].Data(), subdirWidth[iw].Data(), fileunfoutbase.Data(), varyPrior, tag[0].Data()) };
      	  Printf("Combining files: %s and %s", files[0].Data(), files[1].Data());
      	  
      	  
      	  if (dounfold) {
      	  	  
      	  	  Printf("%d ) Running combine triggers PtBin W %.0f, pT min = %.0f, M min = %.0f, pTtrue min = %.0f, Mtrue min = %.0f in directory  %s", iw, binWidthPtvar[iw], pt_min[centralrec], m_min[centralrec],  pt_minT[centralgen], m_minT[centralgen], gSystem->pwd());
      	  	  
      	      
      	  	  if(!dryrun) {
      	       	  
      	       	  TFile testfile(outfileunfmergename);
      	       	  if(!forcedounfold && testfile.IsOpen()){
      	       	  	  Printf("Unfolded output %s already present, don't run again", outfileunfmergename.Data() );
      	       	  }else {
      	       	  	  CombineTriggersCorrected(ntr, files, arrayiter, namescombo);
      	       	  }
      	      } //dry run
      	  } //do unfold
      } // bin width variation
   } //loop on response directories (main ones, need to add the subdir of the variation)
   
   Printf("\nCombine triggers for Lower Traking efficiency in data and response %s/%s", wdirunfLowTrEffData.Data(), subdir[centralvar].Data());
   gSystem->ChangeDirectory(wdirunfLowTrEffData);
   gSystem->mkdir(subdir[centralvar]);
   gSystem->ChangeDirectory(subdir[centralvar]);
      	  	  
   TString files[ntr] = {Form("%s/%s/%s%d%s.root", wdirrespLowTrEff[0].Data(), subdir[centralvar].Data(), fileunfoutbase.Data(), varyPrior, tag[idTrE].Data()),  Form("%s/%s/%s%d%s.root", wdirrespLowTrEff[1].Data(), subdir[centralvar].Data(), fileunfoutbase.Data(), varyPrior, tag[idTrE].Data()) };
   Printf("Combining files: %s and %s", files[0].Data(), files[1].Data());
      	  
   if (dounfold) {
      	  	  
   	   Printf("Running combine triggers lowerTr eff in data and response in directory  %s", gSystem->pwd());
      	  	  
   	   if(!dryrun) {
      	       	  
   	   	   TFile testfile(outfileunfmergename);
   	   	   if(!forcedounfold && testfile.IsOpen()){
   	   	   	   Printf("Unfolded output %s already present, don't run again", outfileunfmergename.Data() );
   	   	   }else {
   	   	   	   CombineTriggersCorrected(ntr, files, arrayiter, namescombo);
   	   	   }
   	   } //dry run
   } //do unfold      	  
      	  
   if(dryrun) return;
   
   Int_t nx, ny, dx, dy;
   CalculatePads(nptbins, nx, ny, dx, dy);
   Int_t massbinW = 2.; //reset below with the right bin width
   
   TH1D *hSystTotForDrawing[nptbins]; // one per pT bin
   TH1D *hStatUncForDrawing[nptbins]; // one per pT bin
   TH1D **hSystPartForDrawing[ncontribution] = {0x0}; // nptbins x contributions
   
   if(docompare){
   	   Printf("\nCalculate the systematics from the variations:");
   	   for(Int_t iv = 0; iv < ncontribution; iv++){
   	   	   if(activecmp[iv]) Printf("- Comparison %d enabled", iv);
   	   	   else Printf("- Comparison %d disabled", iv);
   	   }
   	   
   	   // compare different variations, the paths are defined in the previous loop
   	   TH1D *hSystTot[nptbins]; // one per pT bin
   	   TH1D **hSystPart[ncontribution] = {0x0}; // nptbins x contributions
   	   
   	   
   	   
   	   Int_t act = 0;
   	   if(activecmp[0]){
   	   	   Printf("\nRunning compare background subtraction types");
   	   	   for(Int_t idata = 0; idata < ntr; idata++){
   	   	   	   Printf("%d ) %s", idata, inputcmpBkg[idata].Data());
   	   	   }
   	   	   
   	   	   gSystem->ChangeDirectory(basewdircmp);
   	   	   gSystem->mkdir(wdircmpBkg);
   	   	   gSystem->ChangeDirectory(wdircmpBkg);
   	   	   Printf("in directory  %s", gSystem->pwd());
   	   	   fileSysOutput[act] = gSystem->pwd();
   	   	   fileSysOutput[act] += Form("/Ratio%sOver%s.root", legcmpBkg[0].Data(), legcmpBkg[1].Data());
   	   	   hnameSysMean[act] = "hSysMeanBkgSub";
   	   	   act++;
   	   	   if(!dryrun) {
   	   	   	   TH1D** htmpSysPart = CompareResults(nbkg, inputcmpBkg, inputhnamecmpBkg, legcmpBkg, offsetcmpBkg, kTRUE, kTRUE, kFALSE, 1, "BkgSub");
   	   	   	   
   	   	   	   //this part below smooths the systematics for this source
   	   	   	   TString namestobereused[nptbins];
   	   	   	   
   	   	   	   Int_t arrayusedptbins[nptbins] = {1, 1, 0};
   	   	   	   for(Int_t ipt = 0; ipt<nptbins; ipt++){
   	   	   	   	   namestobereused[ipt] = htmpSysPart[ipt]->GetName();
   	   	   	   	   if(!arrayusedptbins[ipt]) htmpSysPart[ipt] = 0x0;
   	   	   	   }
   	   	   	   TH1D* hnewSysPart = SmoothUncertaintyByAverage(nptbins, htmpSysPart, Form("nametobeassigned"));
   	   	   	   
   	   	   	   TCanvas *cSysBkgFinal = new TCanvas("cSysBkgFinal", "Final systematics background", dx, dy); 
   	   	   	   cSysBkgFinal->Divide(nx, ny);
   	   	   	   hSystPart[0] = new TH1D*[nptbins];
   	   	   	   hSystPartForDrawing[0] = new TH1D*[nptbins];
   	   	   	   for(Int_t ipt = 0; ipt< nptbins; ipt++){
   	   	   	   	   hSystPart[0][ipt] = new TH1D(*hnewSysPart);
   	   	   	   	   hSystPart[0][ipt]->SetName(namestobereused[ipt]);
   	   	   	   	   Printf("isys 0 ---> %p (%s)", hSystPart[0][ipt], hSystPart[0][ipt]->GetName());
   	   	   	   	   cSysBkgFinal->cd(ipt+1);
   	   	   	   	   hSystPart[0][ipt]->Draw("E2");
   	   	   	   	   //delete htmpSysPart[ipt];
   	   	   	   	   hSystPartForDrawing[0][ipt] = (TH1D*)hSystPart[0][ipt]->Clone(Form("%sdraw", hSystPart[0][ipt]->GetName()));
   	   	   	   	   hSystPartForDrawing[0][ipt]->Reset();
   	   	   	   	   
   	   	   	   }
   	   	   	   delete hnewSysPart;
   	   	   	   SaveCv(cSysBkgFinal);
   	   	   }
      }
      
      if(activecmp[1]){
      	 Printf("\nRunning compare matrix ranges");
      	 hSystPartForDrawing[1] = new TH1D*[nptbins];
      	 for(Int_t index = 0; index <nranges; index++){
      	 	 
      	    Printf("%d ) %s, ", index, inputUnfSystRanges[index].Data());
      	    
      	    if(kineeffinrangesyst){
      	 	 	 
      	 	 	 TH1D *hKineNum[nptbins];
      	 	 	 TH1D *hKineDen[nptbins];
      	 	 	 TH1D *hKineEff[nptbins];
      	 	 	 TFile *feffMB = new TFile(Form("%s/%s/%s", defaultRespPaths[0].Data(), subdir[index].Data(),  kineeffMBFilenameSystRanges[index].Data()));
      	 	 	 TFile *feffJE = new TFile(Form("%s/%s/%s", defaultRespPaths[1].Data(), subdir[index].Data(), kineeffJEFilenameSystRanges[index].Data()));
      	 	 	 if(!feffJE->IsOpen() || !feffMB->IsOpen()){
      	 	 	 	 Printf("******************************************************No Kine eff found, will be skipped");
      	 	 	 	 for(Int_t ipt = 0; ipt<nptbins; ipt++){
      	 	 	 	 	 hKineNum[ipt] = 0x0;
      	 	 	 	 	 hKineDen[ipt] = 0x0;
      	 	 	 	 }
      	 	 	 } else {
      	 	 	 	 TString basenameNum = "hNumMPtPar";
      	 	 	 	 TString basenameDen = "hDenMPtPar";
      	 	 	 	 
      	 	 	 	 TFile *fMassSysRanges = new TFile(inputUnfSystRanges[index], "update");
      	 	 	 	 inputhnameUnfSystRanges[index] = Form("%sEfKCor_Ranges_id%db_Pt", inputhnameUnfSystRanges[index].Data(), index);
      	 	 	 	 
      	 	 	 	 for(Int_t ipt = 0; ipt<nptbins; ipt++){
      	 	 	 	 	 if(ipt < binSwitchToEJE) {
      	 	 	 	 	 	 hKineNum[ipt] = (TH1D*)feffMB->Get(Form("%s%.0f%0.f", basenameNum.Data(), ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 hKineNum[ipt]->SetName(Form("%s_Idx%d_%.0f%0.f", basenameNum.Data(), index, ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 hKineDen[ipt] = (TH1D*)feffMB->Get(Form("%s%.0f%0.f", basenameDen.Data(), ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 hKineDen[ipt]->SetName(Form("%s_Idx%d_%.0f%0.f", basenameDen.Data(), index, ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 
      	 	 	 	 	 } else {
      	 	 	 	 	 	 hKineNum[ipt] = (TH1D*)feffJE->Get(Form("%s%.0f%0.f", basenameNum.Data(), ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 hKineNum[ipt]->SetName(Form("%s_Idx%d_%.0f%0.f", basenameNum.Data(), index, ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 hKineDen[ipt] = (TH1D*)feffJE->Get(Form("%s%.0f%0.f", basenameDen.Data(), ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 hKineDen[ipt]->SetName(Form("%s_Idx%d_%.0f%0.f", basenameDen.Data(), index, ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 
      	 	 	 	 	 }
      	 	 	 	 	 
      	 	 	 	 	 TH1D* hMass = 0x0;
      	 	 	 	 	 hMass = (TH1D*)fMassSysRanges->Get(Form("%s%d", basehnname.Data(), ipt));
      	 	 	 	 	 if(!hMass) {
      	 	 	 	 	 	 Printf("%s%d not found", basehnname.Data(), ipt);
      	 	 	 	 	 	 fMassSysRanges->ls();
      	 	 	 	 	 	 continue;
      	 	 	 	 	 }
      	 	 	 	 	 
      	 	 	 	 	 
     	 	 	 	 	 TH1D* hMassK= 0x0;
      	 	 	 	 	 if(hKineNum[ipt] && hKineDen[ipt]) {
      	 	 	 	 	 	 hMassK = ApplyKineEff(hMass, hKineNum[ipt], hKineDen[ipt]);
      	 	 	 	 	 	 if(!hMassK) {
      	 	 	 	 	 	 	 Printf("00000000000 Problem");
      	 	 	 	 	 	 	 continue;
      	 	 	 	 	 	 }
      	 	 	 	 	 	 hMassK->SetName(Form("%s%d",  inputhnameUnfSystRanges[index].Data(), ipt));
      	 	 	 	 	 	 hMassK->SetMarkerColor(kRed+2);
      	 	 	 	 	 	 hMassK->SetLineColor(kRed+2);
      	 	 	 	 	 	 hMassK->GetYaxis()->SetRangeUser(0, 0.25);
      	 	 	 	 	 	 // 2nd mass corrected for kine eff (rescaled to integral 1)
      	 	 	 	 	 	 
      	 	 	 	 	 	 hKineEff[ipt] = (TH1D*) hKineNum[ipt]->Clone(Form("hKineEff_Idx%d_%.0f%0.f", index, ptlims[ipt], ptlims[ipt+1]));
      	 	 	 	 	 	 hKineEff[ipt]->Divide(hKineDen[ipt]);
      	 	 	 	 	 	 
      	 	 	 	 	 	 Int_t lastBabove = hKineEff[ipt]->FindLastBinAbove(0.7);
      	 	 	 	 	 	 if(lastBabove < 1) lastBabove = 1;
      	 	 	 	 	 	 TParameter<Double_t> *param = new TParameter<Double_t>(Form("lastXAbove_Idx%d_%.0f%0.f", index, ptlims[ipt], ptlims[ipt+1]), hKineEff[ipt]->GetBinCenter(lastBabove));
      	 	 	 	 	 	 
      	 	 	 	 	 	 
      	 	 	 	 	 	 fMassSysRanges->cd();
      	 	 	 	 	 	 hMassK->Write();
      	 	 	 	 	 	 param->Write();
      	 	 	 	 	 }
      	 	 	 	 }
      	 	 	 	 //Printf("Check if %s is there", inputhnameUnfSystRanges[index].Data());
      	 	 	 	 //fMassSysRanges->ls();
      	 	 	 	 fMassSysRanges->Close();
      	 	 	 	 delete fMassSysRanges;
      	 	 	 }
      	 	 	 
      	 	 }
      	 }
      	 
      	 gSystem->ChangeDirectory(basewdircmp);
      	 gSystem->mkdir(wdirUnfSystRanges);
      	 gSystem->ChangeDirectory(wdirUnfSystRanges);
      	 Printf("in directory  %s", gSystem->pwd());
      	 fileSysOutput[act] = gSystem->pwd();
      	 fileSysOutput[act] += Form("/Ratio%sOver%s.root", legUnfSystRanges[0].Data(), legUnfSystRanges[1].Data());
      	 hnameSysMean[act] = "hSysMeanRanges";
      	 act++;
      	 if(!dryrun) {
      	 	       	 	     	 	 
      	 	 TH1D** htmpSysPart = CompareResults(nranges, inputUnfSystRanges, inputhnameUnfSystRanges, legUnfSystRanges, offsetUnfSystRanges, kTRUE, kTRUE, kFALSE, centralvar, "Ranges");
      	 	 //Printf("isys 1 ---> %p (%s)", htmpSysPart[1][0], htmpSysPart[1][0]->GetName());
      	 	 
      	 	 hSystPart[1] = new TH1D*[nptbins];
      	 	 TH1D** hnewSysPart = SmoothUnceraintyByFit(nptbins, htmpSysPart, fittype[1], Form("h%ssmoot", description[1].Data()), hSystPart[1]);
      	 	 
      	 	 for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	 	 	 hSystPartForDrawing[1][ipt] = (TH1D*)hSystPart[1][ipt]->Clone(Form("%sdraw", hSystPart[1][ipt]->GetName()));
      	 	 }
      	 }
      }
      
      if(activecmp[2]){
      	  hSystPartForDrawing[2] = new TH1D*[nptbins];
      	  Printf("\nRunning compare iterations");
      	  for(Int_t index = 0; index < nitervar; index++){
      	  	  Printf("%d ) iteration %d (File %s) vs default %d (File %s)  ", index, iterminsyst+index, inputUnfSystIter[index].Data(), iterdef, inputUnfSystIter[iterdef-iterminsyst].Data());
      	  }
      	  gSystem->ChangeDirectory(basewdircmp);
      	  gSystem->mkdir(wdirUnfSystIter);
      	  gSystem->ChangeDirectory(wdirUnfSystIter);
      	  Printf("in directory  %s", gSystem->pwd());
      	  fileSysOutput[act] = gSystem->pwd();
      	  fileSysOutput[act] += Form("/Ratio%sOver%s.root", legUnfSystIter[0].Data(), legUnfSystIter[1].Data());
      	  hnameSysMean[act] = "hSysMeanIter";
      	  act++;
      	  if(!dryrun) {
      	  	  
      	  	  TH1D** htmpSysPart = CompareResults(nitervar, inputUnfSystIter, inputhnameUnfSystIter, legUnfSystIter, offsetUnfSystIter, kTRUE, kTRUE, kFALSE, iterdef-iterminsyst, "Iter");
      	  	  
      	  	  hSystPart[2]  = new TH1D*[nptbins];
      	  	  
      	  	  
      	  	  //TH1D** hnewSysPart = SmoothUnceraintyByFit(nptbins, htmpSysPart, fittype[2], Form("h%ssmoot", description[2].Data()), hSystPart[2], 2);
      	  	  
      	  	  for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	  	  	  hSystPart[2][ipt] = htmpSysPart[ipt];
      	  	  	  
      	  	  	  hSystPartForDrawing[2][ipt] = (TH1D*)hSystPart[2][ipt]->Clone(Form("%sdraw", hSystPart[2][ipt]->GetName()));
      	  	  }
      	  }
      }
      
      if(activecmp[3]){
      	  hSystPartForDrawing[3] = new TH1D*[nptbins];
   	   	   Printf("\nRunning compare Overlap/no Overlap");
   	   	   for(Int_t idata = 0; idata <novl; idata++){
   	   	   	   Printf("%d ) %s, ", idata, inputcmpOvl[idata].Data());
   	   	   }
   	   	   
   	   	   gSystem->ChangeDirectory(basewdircmp);
   	   	   gSystem->mkdir(wdircmpOvl);
   	   	   gSystem->ChangeDirectory(wdircmpOvl);
   	   	   Printf("in directory  %s", gSystem->pwd());
   	   	   fileSysOutput[act] = gSystem->pwd();
   	   	   fileSysOutput[act] += Form("/Ratio%sOver%s.root", legcmpOvl[0].Data(), legcmpOvl[1].Data());
   	   	   hnameSysMean[act] = "hSysMeanOvlExl";
   	   	   act++;
   	   	   if(!dryrun) {
   	   	   	   TH1D** htmpSysPart = CompareResults(novl, inputcmpOvl, inputhnamecmpOvl, legcmpOvl, offsetcmpOvl, kTRUE, kTRUE, kFALSE, 0, "OvlExl");
   	   	   	   
   	   	   	   hSystPart[3] = new TH1D*[nptbins];
   	   	   	   TH1D** hnewSysPart = SmoothUnceraintyByFit(nptbins, htmpSysPart, fittype[3], Form("h%ssmoot", description[3].Data()), hSystPart[3]);
   	   	   	   
   	   	   	   for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	 	 	 hSystPartForDrawing[3][ipt] = (TH1D*)hSystPart[3][ipt]->Clone(Form("%sdraw", hSystPart[3][ipt]->GetName()));
      	 	 }
      	 }
      }
      
      if(activecmp[4]){
      	  hSystPartForDrawing[4] = new TH1D*[nptbins];
   	   	   Printf("\nRunning compare PtBinWidth");
   	   	   for(Int_t iw = 0; iw < widthvar; iw++){
   	   	   	   Printf("%d ) %s, ", iw, inputcmpBinW[iw].Data());
   	   	   }
   	   	   
   	   	   gSystem->ChangeDirectory(basewdircmp);
   	   	   gSystem->mkdir(wdirUnfSystBinW);
   	   	   gSystem->ChangeDirectory(wdirUnfSystBinW);
   	   	   Printf("in directory  %s", gSystem->pwd());
   	   	   fileSysOutput[act] = gSystem->pwd();
   	   	   fileSysOutput[act] += Form("/Ratio%sOver%s.root", legUnfSystBinW[0].Data(), legUnfSystBinW[1].Data());
   	   	   hnameSysMean[act] = "hSysMeanBinW";
   	   	   act++;
   	   	   if(!dryrun) {
   	   	   	   TH1D** htmpSysPart = CompareResults(widthvar, inputcmpBinW, inputhnameBinW, legUnfSystBinW, offsetUnfSystBinW, kTRUE, kTRUE, kFALSE, 0, "BinW");
   	   	   	   
   	   	   	   hSystPart[4] = new TH1D*[nptbins];
   	   	   	   TH1D** hnewSysPart = SmoothUnceraintyByFit(nptbins, htmpSysPart, fittype[4], Form("h%ssmoot", description[4].Data()), hSystPart[4]);
   	   	   	   Printf("isys 0 ---> %p (%s)", hSystPart[4][0], hSystPart[4][0]->GetName());
   	   	   	  
   	   	   	   for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	 	 	 hSystPartForDrawing[4][ipt] = (TH1D*)hSystPart[4][ipt]->Clone(Form("%sdraw", hSystPart[4][ipt]->GetName()));
      	 	 }
      	 }
      }
      
      
      if(activecmp[5]){
      	  hSystPartForDrawing[5] = new TH1D*[nptbins];
   	   	   Printf("\nRunning compare Reduced Tracking Efficiency in response");
   	   	   for(Int_t idata = 0; idata <nTrkE; idata++){
   	   	   	   Printf("%d ) %s, ", idata, inputcmpTrkE[idata].Data());
   	   	   }
   	   	   
   	   	   gSystem->ChangeDirectory(basewdircmp);
   	   	   gSystem->mkdir(wdircmpTrkE);
   	   	   gSystem->ChangeDirectory(wdircmpTrkE);
   	   	   Printf("in directory  %s", gSystem->pwd());
   	   	   fileSysOutput[act] = gSystem->pwd();
   	   	   fileSysOutput[act] += Form("/Ratio%sOver%s.root", legcmpTrkE[0].Data(), legcmpTrkE[1].Data());
   	   	   hnameSysMean[act] = "hSysMean96TrkEff";
   	   	   act++;
   	   	   if(!dryrun) {
   	   	   	   TH1D** htmpSysPart = CompareResults(nTrkE, inputcmpTrkE, inputhnamecmpTrkE, legcmpTrkE, offsetcmpTrkE, kTRUE, kTRUE, kFALSE, 0, "96TrkEff");
   	   	   	   hSystPart[5] = new TH1D*[nptbins];
   	   	   	   TH1D** hnewSysPart = SmoothUnceraintyByFit(nptbins, htmpSysPart, fittype[5], Form("h%ssmoot", description[5].Data()), hSystPart[5]);
   	   	   	   for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	 	 	 hSystPartForDrawing[5][ipt] = (TH1D*)hSystPart[5][ipt]->Clone(Form("%sdraw", hSystPart[5][ipt]->GetName()));
      	 	 }
      	 }
      }
      
      if(activecmp[6]){
      	  hSystPartForDrawing[6] = new TH1D*[nptbins];
      	  Printf("\nRunning compare Reduced Tracking Efficiency in data and response");
      	  for(Int_t idata = 0; idata <nTrkE; idata++){
      	  	  Printf("%d ) %s, ", idata, inputcmpTrkE[idata].Data());
      	  }
      	  
      	  gSystem->ChangeDirectory(basewdircmp);
      	  gSystem->mkdir(wdircmpTrkEv2);
      	  gSystem->ChangeDirectory(wdircmpTrkEv2);
      	  Printf("in directory  %s", gSystem->pwd());
      	  fileSysOutput[act] = gSystem->pwd();
      	  fileSysOutput[act] += Form("/Ratio%sOver%s.root", legcmpTrkEv2[0].Data(), legcmpTrkEv2[1].Data());
      	  hnameSysMean[act] = "hSysMean96TrkEffv2";
      	  act++;
      	  if(!dryrun) {
      	  	  TH1D** htmpSysPart = CompareResults(nTrkE, inputcmpTrkEv2, inputhnamecmpTrkE, legcmpTrkEv2, offsetcmpTrkE, kTRUE, kTRUE, kFALSE, 0, "96TrkEffv2");
      	  	  
      	  	  hSystPart[6] = new TH1D*[nptbins];
   	   	   	   TH1D** hnewSysPart = SmoothUnceraintyByFit(nptbins, htmpSysPart, fittype[6], Form("h%ssmoot", description[6].Data()), hSystPart[6]);
   	   	   	   
      	  	  for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	  	  	  hSystPartForDrawing[6][ipt] = (TH1D*)hSystPart[6][ipt]->Clone(Form("%sdraw", hSystPart[6][ipt]->GetName()));
      	  	  	  hSystPartForDrawing[6][ipt]->Reset();
      	  	  }
      	  }
      }
      
      if(activecmp[7]){
      	  hSystPartForDrawing[7] = new TH1D*[nptbins];
      	  Printf("\nRunning compare Reduced Tracking Efficiency in data and response");
      	  for(Int_t idata = 0; idata <nvarprior; idata++){
      	  	  Printf("%d ) %s, ", idata, inputcmpPrior[idata].Data());
      	  }
      	  
      	  gSystem->ChangeDirectory(basewdircmp);
      	  gSystem->mkdir(wdircmpPrior);
      	  gSystem->ChangeDirectory(wdircmpPrior);
      	  Printf("in directory  %s", gSystem->pwd());
      	  fileSysOutput[act] = gSystem->pwd();
      	  fileSysOutput[act] += Form("/Ratio%sOver%s.root", legcmpPrior[0].Data(), legcmpPrior[1].Data());
      	  hnameSysMean[act] = "hSysMeanPrior";
      	  act++;
      	  if(!dryrun) {
      	  	  TH1D** htmpSysPart = CompareResults(nvarprior, inputcmpPrior, inputhnamecmpPrior, legcmpPrior, offsetcmpPrior, kTRUE, kTRUE, kFALSE, 0, "Prior");
      	  	  
      	  	  hSystPart[7] = new TH1D*[nptbins];
   	   	   	   
      	  	  //TH1D** hnewSysPart = SmoothUnceraintyByFit(nptbins, htmpSysPart, fittype[7], Form("h%ssmoot", description[7].Data()), hSystPart[7]);
   	   	   	   
      	  	  for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	  	  	  hSystPart[7][ipt] = htmpSysPart[ipt];
      	  	  	  hSystPartForDrawing[7][ipt] = (TH1D*)hSystPart[7][ipt]->Clone(Form("%sdraw", hSystPart[7][ipt]->GetName()));
      	  	  	  hSystPartForDrawing[7][ipt]->Reset();
      	  	  }
      	  }
      }
      
      gSystem->ChangeDirectory(basewdircmp);

      TH1D** hMeanSys = GetMeanSyst(nactivecontributions, fileSysOutput, hnameSysMean);
      TH1D* hMeanSysTot = SumMeanSyst(nactivecontributions, hMeanSys);
      //SetMassValueInSystematic(hMeanSysTot, );
      TCanvas *cMeanSystTot = new TCanvas("cMeanSystTot", "Tot Systematics on the mean");
      cMeanSystTot->cd();
      hMeanSysTot->Draw();
      TCanvas *cMeanSyst = new TCanvas("cMeanSyst", "Systematics on the mean");
      cMeanSyst->cd();
      act = 0;
      for(Int_t isys = 0; isys < ncontribution ; isys++){
      	  if(!activecmp[isys]) continue;
      	  if(!hMeanSys[act]) continue;
      	  Printf("Print %d, %d, %p", isys, act, hMeanSys[act]);
      	  hMeanSys[act]->SetLineColor(hSystPart[isys][0]->GetLineColor());
      	  if (act == 0) hMeanSys[act]->Draw();
      	  else hMeanSys[act]->Draw("sames");
      	  act++;
      }
      
      
      TLegend *legSysTot = new TLegend(0.3, 0.15, 0.7, 0.45);
      legSysTot->SetBorderSize(0);
      legSysTot->SetFillStyle(0);
      
      TCanvas *cSyst = new TCanvas("cSyst", "Systematics", dx, dy);
      cSyst->Divide(nx, ny);
      //gStyle->SetEndErrorSize(3.);
      //gStyle->SetErrorX(1.);
      
      TFile *foutSysTot = new TFile(filetotsystname, "recreate");
      
      hMeanSysTot->Write();
      
      for(Int_t ipt = 0 ; ipt< nptbins; ipt++){
      	 hSystTot[ipt] = 0x0;
      	 TH1D *hSyst[ncontribution] = {0x0};
      	 
      	 for(Int_t isys = 0; isys < ncontribution ; isys++){
      	 	 if(!activecmp[isys]) continue;
      	 	 if(!hSystPart[isys]) continue;
      	 	 if(!hSystPart[isys][ipt]) continue;
      	 	 
      	 	 // fill the content of the systematic histogram used for drawing in the paper
      	 	 for(Int_t ib = 0; ib < hSystPartForDrawing[isys][ipt]->GetNbinsX(); ib++){
      	 	 
      	 	 	 hSystPartForDrawing[isys][ipt]->SetBinContent(ib+1, hSystPart[isys][ipt]->GetBinError(ib+1));
      	 	 	 hSystPartForDrawing[isys][ipt]->SetBinError(ib+1, 0);
      	 	 	 
      	 	 }
      	 	 hSystPart[isys][ipt]->SetMarkerSize(1);
      	    if(ipt == 0){
      	    	legSysTot->AddEntry(hSystPart[isys][ipt], description[isys], "l");
      	    }
      	    Printf("%p, isys = %d, ipt = %d", hSystPart[isys][ipt], isys, ipt );
      	    
      	    hSystPart[isys][ipt]->SetLineColor(colors[isys]);
      	    hSystPartForDrawing[isys][ipt]->SetLineColor(colors[isys]);
      	    hSystPartForDrawing[isys][ipt]->SetMarkerStyle(markers[isys]);
      	    hSystPartForDrawing[isys][ipt]->SetMarkerColor(colors[isys]);
      	    cSyst->cd(ipt+1);
      	    if(isys == 0) hSystPart[isys][ipt]->Draw("E2");
      	    else hSystPart[isys][ipt]->Draw("E2sames");
      	    
      	    
      	    hSyst[isys] = hSystPart[isys][ipt];
      	    
      	 }
      	 
      	 hSystTot[ipt] = AddInQuadrature(hSyst, ncontribution, ipt, basenamehsystot);
      	 if(!hSystTot[ipt]) {
      	 	 Printf("Error, total not found!");
      	 	 continue;
      	 }
      	 hSystTot[ipt]->SetLineWidth(2);
      	 hSystTot[ipt]->SetLineColor(kOrange+8);
      	 hSystTot[ipt]->SetMarkerColor(kOrange+8);
      	 if(ipt == 0) legSysTot->AddEntry(hSystTot[ipt], "Total", "l");
      	 
      	 hSystTotForDrawing[ipt] = (TH1D*)hSystTot[ipt]->Clone(Form("%sdraw", hSystTot[ipt]->GetName()));
      	 hSystTotForDrawing[ipt]->Reset();
      	 
      	 
      	 // fill the content of the systematic histogram used for drawing in the paper
      	 for(Int_t ib = 0; ib < hSystTotForDrawing[ipt]->GetNbinsX(); ib++){
      	 	 hSystTotForDrawing[ipt]->SetBinContent(ib+1, hSystTot[ipt]->GetBinError(ib+1));
      	 	 hSystTotForDrawing[ipt]->SetBinError(ib+1, 0);
      	 }
      	    
      	 cSyst->cd(ipt+1);
      	 hSystTot[ipt]->DrawClone("E2sames");
      	 if(ipt == 1) legSysTot->Draw(); 
      	 foutSysTot->cd();
      	 hSystTot[ipt]->Write();
      }
      
      cMeanSyst->cd();
      legSysTot->Draw();
      
      massbinW = hSystTot[0]->GetBinWidth(1);
      SaveCv(cSyst);
      foutSysTot->Close();
   }
   
   if(plotResults){
   	   
   	   // text for the final plots
   	   TPaveText *pvgeneral = new TPaveText(0.58, 0.5, 0.83, 0.6, "NDC");
   	   pvgeneral->SetFillStyle(0);
   	   pvgeneral->SetBorderSize(0);
   	   pvgeneral->AddText("Anti-#it{k}_{T}, #it{R} = 0.4");
   	   
   	   TPaveText *pvpPb = new  TPaveText(0.25, 0.7, 0.6, 0.8, "NDC");
   	   pvpPb->SetFillStyle(0);
   	   pvpPb->SetBorderSize(0);
   	   pvpPb->AddText(Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
   	   
   	   TCanvas *cSystPaper = new TCanvas("cSystPaper", "Systematics", dx, dy);
   	   cSystPaper->Divide(nx, ny);
   	   
   	   TLegend *legSysPaper = new TLegend(0.28, 0.4, 0.68, 0.8);
   	   legSysPaper->SetBorderSize(0);
   	   legSysPaper->SetFillStyle(0);
   	   
   	   TH1D *hKineNum[nptbins];
   	   TH1D *hKineDen[nptbins];
   	   
   	   if(correctforkineeff){
   	   	   TString kineeffMBFilename = Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", pt_min[centralrec], pt_max[centralrec], m_min[centralrec], m_max[centralrec]);
   	   	   
   	   	   TString kineeffJEFilename = Form("KineEffPtR%.0f_%0.fMR%.0f_%0.f.root", pt_mineje[centralrec], pt_maxeje[centralrec], m_mineje[centralrec], m_maxeje[centralrec]);
   	   	   
   	   	   
   	   	   TFile *feffMB = new TFile(kineeffMBFilename);
   	   	   TFile *feffJE = new TFile(kineeffJEFilename);
   	   	   if(!feffMB->IsOpen()){
   	   	   	   // calculate Kinematic efficiencies
   	   	   	   CalculateKineEff(pt_min[centralrec], pt_max[centralrec], m_min[centralrec], m_max[centralrec], strIn[defaultUnfNumb], strL[defaultUnfNumb], "");
   	   	   	   feffMB = new TFile(kineeffMBFilename);
   	   	   }
   	   	   if(!feffJE->IsOpen()){
   	   	   	   // calculate Kinematic efficiencies
   	   	   	   CalculateKineEff(pt_mineje[centralrec], pt_maxeje[centralrec], m_mineje[centralrec], m_maxeje[centralrec], strIn[defaultUnfNumb], strL[defaultUnfNumb], "");
   	   	   	   feffJE = new TFile(kineeffJEFilename);
   	   	   }
   	   	   
   	   	   if(!feffJE->IsOpen() || !feffMB->IsOpen()){
   	   	   	   Printf("******************************************************No Kine eff found, will be skipped");
   	   	   	   for(Int_t ipt = 0; ipt<nptbins; ipt++){
   	   	   	   	   hKineNum[ipt] = 0x0;
   	   	   	   	   hKineDen[ipt] = 0x0;
   	   	   	   }
   	   	   } else {
   	   	   	   TString basenameNum = "hNumMPtPar";
   	   	   	   TString basenameDen = "hDenMPtPar";
   	   	   	   for(Int_t ipt = 0; ipt<nptbins; ipt++){
   	   	   	   	   if(ipt < binSwitchToEJE) {
   	   	   	   	   	   hKineNum[ipt] = (TH1D*)feffMB->Get(Form("%s%.0f%0.f", basenameNum.Data(), ptlims[ipt], ptlims[ipt+1]));
   	   	   	   	   	   hKineDen[ipt] = (TH1D*)feffMB->Get(Form("%s%.0f%0.f", basenameDen.Data(), ptlims[ipt], ptlims[ipt+1]));
   	   	   	   	   } else {
   	   	   	   	   	   hKineNum[ipt] = (TH1D*)feffJE->Get(Form("%s%.0f%0.f", basenameNum.Data(), ptlims[ipt], ptlims[ipt+1]));
   	   	   	   	   	   hKineDen[ipt] = (TH1D*)feffJE->Get(Form("%s%.0f%0.f", basenameDen.Data(), ptlims[ipt], ptlims[ipt+1]));
   	   	   	   	   }
   	   	   	   }
   	   	   	   
   	   	   }
   	   } else {
   	   	   Printf("Setting: No Kine Eff");
   	   	   for(Int_t ipt = 0; ipt<nptbins; ipt++){
   	   	   	   hKineNum[ipt] = 0x0;
   	   	   	   hKineDen[ipt] = 0x0;
   	   	   }
   	   }
   	   //CorrectForKineEff(outfileunfmergename, kineeffFilename, "", pt_min[centralrec], pt_max[centralrec], m_min[centralrec], m_max[centralrec]);
   	   
   	   TCanvas *cResults = new TCanvas("cResults", "Mass spectra with systematic uncertainties", 1100, 800);
   	   cResults->Divide(nx, ny);
   	   //TCanvas *cRelUnc = new TCanvas("cRelUnc", "Relative uncertainties", nx, ny);
   	   //cRelUnc->Divide(nx, ny);
   	   
   	   // PYTHIA 
   	   TString listnamePythia = "JetMassResponseDet_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme";
   	   TString pathPythia = "/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/AnalysisResults.root";

   	   TString thnspname = "fhnMassResponse", projnamepythia = "hMPythiaPar";
   	   Int_t axrange = 3, axproj = 1;
   	   
   	   TH1D** hMassPythia = GetPythiaOrThnSpaseProjections(axrange, axproj, massbinW, nptbins, ptlims, pathPythia, listnamePythia, thnspname, projnamepythia);
   	   if(!hMassPythia) Printf("Not possible to find pythia");
   	   
   	   //method
   	   TList *results = AddSystematicstoMassFromFile(Form("%s/%s/%s", defaultUnfPaths.Data(), subdir[centralvar].Data(), outfileunfmergename.Data()), basehnname, kFALSE, Form("%s/%s", basewdircmp.Data(), filetotsystname.Data()), Form("%s_Pt", basenamehsystot.Data()), kTRUE);
   	   
   	   TFile *fMassResults = new TFile("MasspPbResults.root", "recreate");
   	   TFile *fFigPaper = new TFile("InputFiguresPaperSyst.root", "recreate");
   	   // end method
   	   
   	   //manual
   	   /*
   	   TFile *finM = new TFile(Form("%s/%s/%s", defaultUnfPaths.Data(), subdir[centralvar].Data(), outfileunfmergename.Data()));
   	   TFile *finSys = new TFile(Form("%s/%s", basewdircmp.Data(), filetotsystname.Data()));
   	   */
   	   // end manual
   	   
   	   TLegend *legRes = new TLegend(0.4, 0.2, 0.9, 0.8);
   	   legRes->SetFillStyle(0);
   	   TPaveText **pvpt = new TPaveText*[nptbins];
   	   for(Int_t ipt = 0; ipt < nptbins; ipt++){
   	   	   
   	   	   pvpt[ipt] = new TPaveText(0.33, 0.8, 0.83, 0.9, "NDC");
   	   	   pvpt[ipt]->SetFillStyle(0);
   	   	   pvpt[ipt]->SetBorderSize(0);
   	   	   pvpt[ipt]->AddText(Form("%.0f < #it{p}_{T, ch jet} (GeV/#it{c}) < %.0f", ptlims[ipt], ptlims[ipt+1]));
   	   	   //method
   	   	   
   	   	   Int_t hnumbperptbin = 4; // check this in AddSystematicstoMassFromFile
   	   	   TH1D* hMass = (TH1D*)results->At(ipt*hnumbperptbin);
   	   	   //Printf("Reading result mass histogram %s, mean %f", hMass->GetName(), hMass->GetMean());
   	   	   hMass->SetName(Form("hUnfM_Itr%d_Pt%.0f_%.0f", iterdef, ptlims[ipt],ptlims[ipt+1]));
   	   	   hMass->SetMarkerStyle(21);
   	   	   hMass->SetMarkerColor(kOrange+8);
   	   	   hMass->SetLineColor(hMass->GetMarkerColor());
   	   	   
   	   	   //copy the histogram of the mass for its stat unc to be drawn in the systematics
   	   	   
   	   	   hStatUncForDrawing[ipt] = (TH1D*)hMass->Clone(Form("%sdraw", hMass->GetName()));
   	   	   hStatUncForDrawing[ipt]->GetXaxis()->SetTitle("#it{M}_{ch jet} (GeV/#it{c}^{2})");
   	   	   hStatUncForDrawing[ipt]->GetYaxis()->SetTitle("Relative Difference");
   	   	   hStatUncForDrawing[ipt]->GetYaxis()->SetTitleOffset(1.45);
   	   	   hStatUncForDrawing[ipt]->Reset();
   	   	   //stat unc as content
   	   	   for(Int_t ib = 0; ib < hStatUncForDrawing[ipt]->GetNbinsX(); ib++){
   	   	   	   
   	   	   	   hStatUncForDrawing[ipt]->SetBinContent(ib+1, hMass->GetBinError(ib+1)/hMass->GetBinContent(ib+1));
   	   	   	   hStatUncForDrawing[ipt]->SetBinError(ib+1, 0);
   	   	   	   
   	   	   }
   	   	   hStatUncForDrawing[ipt]->SetFillColor(kAzure+5);
   	   	   hStatUncForDrawing[ipt]->SetLineColor(kAzure+5);
   	   	   hStatUncForDrawing[ipt]->SetLineWidth(1);
   	   	   hStatUncForDrawing[ipt]->SetBarWidth(0.8);
   	   	   hStatUncForDrawing[ipt]->SetBarOffset(0.2);
   	   	   	   	   
   	   	   TH1D* hMaSy = (TH1D*)results->At(ipt*hnumbperptbin + 1);
   	   	   hMaSy->SetMarkerStyle(1);
   	   	   hMaSy->SetMarkerColor(kOrange+8);
   	   	   hMaSy->SetLineColor(hMass->GetMarkerColor());
   	   	   hMaSy->SetFillStyle(0);
   	   	   TH1D* hReSy = (TH1D*)results->At(ipt*hnumbperptbin + 2);
   	   	   TH1D* hReSt = (TH1D*)results->At(ipt*hnumbperptbin + 3);
   	   	   
   	   	   //1st = mass 
   	   	   fMassResults->cd();
   	   	   hMass->Write();
   	   	   
   	   	   TH1D* hMassK= 0x0;
   	   	   if(hKineNum[ipt] && hKineDen[ipt]) {
   	   	   	   hMassK = ApplyKineEff(hMass, hKineNum[ipt], hKineDen[ipt]);
   	   	   	   hMassK->SetName(Form("hUnfEffCorM_Itr%d_Pt%.0f_%.0f", iterdef, ptlims[ipt],ptlims[ipt+1]));
   	   	   	   hMassK->SetMarkerColor(kRed+2);
   	   	   	   hMassK->SetLineColor(kRed+2);
   	   	   	   hMassK->GetYaxis()->SetRangeUser(0, 0.25);
   	   	   	   // 2nd mass corrected for kine eff (rescaled to integral 1)
   	   	   	   fMassResults->cd();
   	   	   	   hMassK->Write();
   	   	   }
   	   	   
   	   	   // 3rd systematic uncertainties
   	   	   fMassResults->cd();
   	   	   hMaSy->Write();
   	   	   
   	   	   hMass->GetYaxis()->SetRangeUser(0, 0.25);
   	   	   cResults->cd(ipt+1);
   	   	   hMass->Draw();
   	   	   if(hMassK) hMassK->Draw("Esames");
   	   	   hMaSy->Draw("E2sames");
   	   	   //cRelUnc->cd(ipt+1);
   	   	   //hReSy->Draw("E2");
   	   	   //hReSt->Draw("sames");
   	   	   // end method  	   	   
   	   	   if(hMassPythia[ipt]){
   	   	   	   cResults->cd(ipt+1);
   	   	   	   hMassPythia[ipt]->SetLineColor(kBlack);
   	   	   	   hMassPythia[ipt]->SetMarkerColor(kBlack);
   	   	   	   hMassPythia[ipt]->SetMarkerStyle(24);
   	   	   	   hMassPythia[ipt]->Scale(hMass->Integral("width")/hMassPythia[ipt]->Integral("width"));

   	   	   }
 
   	   	   TLine *lmaxM = 0x0;
   	   	   if(ipt<2) lmaxM = new TLine(m_max[centralrec], 0, m_max[centralrec], 0.25);
   	   	   else lmaxM = new TLine(m_maxeje[centralrec], 0, m_maxeje[centralrec], 0.25);
   	   	   
   	   	   lmaxM->SetLineWidth(2);
   	   	   lmaxM->SetLineStyle(2);
   	   	   lmaxM->SetLineColor(kGray);
   	   	   
   	   	   cResults->cd(ipt+1);
   	   	   lmaxM->Draw();
   	   	   
   	   	   
   	   	   if(ipt == 0) {
   	   	   	   legRes->AddEntry(hMass, "Jet Mass pPb", "LP");
   	   	   	   if(hMassK) legRes->AddEntry(hMassK, "Corrected for KineEff", "LP");
   	   	   	   legRes->AddEntry(hMaSy, "Total Syst Unc", "F");
   	   	   	   legRes->AddEntry(hMassPythia[ipt], "PYTHIA", "LP");
   	   	   	   legRes->AddEntry(lmaxM, "Max #it{M}_{rec} unfolding", "l");
   	   	   	   legSysPaper->AddEntry(hStatUncForDrawing[ipt], "Statistical", "F");
   	   	   	   legSysPaper->AddEntry(hSystTotForDrawing[ipt], "Total", "l");
   	   	   }
   	   	   
   	   	   cSystPaper->cd(ipt+1);
   	   	   gPad->SetBottomMargin(.12);
   	   	   gPad->SetLeftMargin(.13);
   	   	   hStatUncForDrawing[ipt]->GetXaxis()->SetLabelSize(0.05);
   	   	   hStatUncForDrawing[ipt]->GetXaxis()->SetTitleOffset(1.1);
   	   	   hStatUncForDrawing[ipt]->GetXaxis()->SetTitleSize(0.048);
   	   	   
   	   	   hStatUncForDrawing[ipt]->GetYaxis()->SetLabelSize(0.05);
   	   	   hStatUncForDrawing[ipt]->GetYaxis()->SetTitleOffset(1.40);
   	   	   hStatUncForDrawing[ipt]->GetYaxis()->SetTitleSize(0.048);
   	   	   hStatUncForDrawing[ipt]->GetYaxis()->SetRangeUser(0, 0.6);
   	   	   hStatUncForDrawing[ipt]->GetXaxis()->SetRangeUser(0, maxRangeMassFinal[ipt]);
   	   	   hStatUncForDrawing[ipt]->Draw("B");
   	   	   for(Int_t isys = 0; isys < ncontribution; isys++){
   	   	   	   if(!activecmp[isys]) continue;
   	   	   	   if(ipt == 0) legSysPaper->AddEntry(hSystPartForDrawing[isys][ipt], description[isys], "l");
   	   	   	   hSystPartForDrawing[isys][ipt]->GetXaxis()->SetRangeUser(0, maxRangeMassFinal[ipt]);
   	   	   	   hSystPartForDrawing[isys][ipt]->Draw("PLsames");
   	   	   	   fFigPaper->cd();
   	   	   	   hSystPartForDrawing[isys][ipt]->Write();
   	   	   }
   	   	   hSystTotForDrawing[ipt]->GetXaxis()->SetRangeUser(0, maxRangeMassFinal[ipt]);
   	   	   hSystTotForDrawing[ipt]->Draw("PLsames");
   	   	   pvpt[ipt]->Draw();
   	   	   
   	   	   if(ipt == 0) {
   	   	   	   cSystPaper->cd(ipt+1);
   	   	   	   legSysPaper->Draw();
   	   	   }
   	   	   if(ipt == 1) {
   	   	   	   pvpPb->Draw();
   	   	   	   pvgeneral->Draw();	   
   	   	   }
   	   	   fFigPaper->cd();
   	   	   hStatUncForDrawing[ipt]->Write();
   	   	   hSystTotForDrawing[ipt]->Write();
   	   	   legSysPaper->Write();
   	   	   //manual
   	   	   /*
   	   	   TString hname = Form("hMUnf_%d_Iter%d", ipt, iterdef);
   	   	   TH1D* hMass = (TH1D*)finM->Get(hname);
   	   	   if(!hMass){
   	   	   	   Printf("%s not found, continue", hname.Data());
   	   	   	   continue;
   	   	   }
   	   	   TString hsysname = Form("%s_Pt%.0f_%.0f", basenamehsystot.Data(), ptlims[ipt],ptlims[ipt+1]); 
   	   	   TH1D* hsysTot = (TH1D*)finSys->Get(hsysname);
   	   	   if(!hsysTot){
   	   	   	   Printf("%s not found, continue", hsysname.Data());
   	   	   	   continue;
   	   	   }
   	   	   hsysTot->SetMarkerStyle(1);
   	   	   for(Int_t ib = 0; ib<hsysTot->GetNbinsX(); ib++){
   	   	   	   
   	   	   	   hsysTot->SetBinContent(ib+1, hMass->GetBinContent(ib+1));
   	   	   	   
   	   	   }
   	   	   
   	   	   
   	   	   cResults->cd(ipt+1);
   	   	   hMass->Draw("P");
   	   	   hsysTot->Draw("E2sames");
   	   	   */
   	   	   // end manual
   	   }
   	   
   	   cSystPaper->cd(nptbins+1);
   	   DrawLogo(1, 0.45, 0.7, 0.85, 0.8, "", 42, "");
   	   
   	   fMassResults->Close();
   	   fFigPaper->Close();
   	   //Marta's results
   	   TString pathMarta = "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/MartapPb/TranformedIntoTH1.root";
   	   //PbPb results
   	   TString pathMartaPbPb = "/data/Work/jets/JetMass/PbPbResults/TranformedIntoTH1.root";
   	   
   	   const Int_t nhM = 3;//4;
   	   TH1D* hMassMStat[nhM];
   	   TH1D* hMassMSyst[nhM];
   	   TH1D* hMassMPbPbStat[nhM];
   	   TH1D* hMassMPbPbSyst[nhM];
   	   Int_t offset = 2; //check this: it's the bin we start from
   	   TString nameMSt = "hUnfMass_PtBin";
   	   TString nameMSy = "hUnfMassSyst_PtBin";
   	   
   	   TFile *fResM = new TFile(pathMarta);
   	   if(!fResM->IsOpen()){
   	   	   Printf("File %s not found", pathMarta.Data());
   	   	   return;
   	   }
   	   TFile *fResPbPbM = new TFile(pathMartaPbPb);
   	   if(!fResPbPbM->IsOpen()){
   	   	   Printf("File %s not found", pathMartaPbPb.Data());
   	   	   return;
   	   }
   	   
   	   for(Int_t ih = 0; ih<nhM ; ih++){
   	   	   hMassMStat[ih] = (TH1D*)fResM->Get(Form("%s%d", nameMSt.Data(), ih+offset));
   	   	   hMassMSyst[ih] = (TH1D*)fResM->Get(Form("%s%d", nameMSy.Data(),ih+offset));
   	   	   
   	   	   hMassMStat[ih]->SetName(Form("pPbM%s%d", nameMSt.Data(), ih+offset));
   	   	   hMassMStat[ih]->SetMarkerColor(kBlue+3);
   	   	   hMassMStat[ih]->SetLineColor(kBlue+3);
   	   	   hMassMSyst[ih]->SetName(Form("pPbM%s%d", nameMSy.Data(), ih+offset));
   	   	   hMassMSyst[ih]->SetFillColor(kAzure+1);
   	   	   hMassMSyst[ih]->SetMarkerColor(kBlue+3);
   	   	   
   	   	   hMassMPbPbStat[ih] = (TH1D*)fResPbPbM->Get(Form("%s%d", nameMSt.Data(), ih+offset));
   	   	   hMassMPbPbSyst[ih] = (TH1D*)fResPbPbM->Get(Form("%s%d", nameMSy.Data(),ih+offset));
   	   	   hMassMPbPbStat[ih]->SetMarkerColor(kOrange+3);
   	   	   hMassMPbPbStat[ih]->SetLineColor(kOrange+3);
   	   	   hMassMPbPbSyst[ih]->SetFillColor(kYellow+1);
   	   	   hMassMPbPbSyst[ih]->SetMarkerColor(kOrange+3); 	   
   	   	   
   	   	   cResults->cd(ih+1);
   	   	   
   	   	   hMassMSyst[ih]->Draw("E2sames");
   	   	   hMassMStat[ih]->Draw("sames");
   	   	   
   	   	   hMassMPbPbSyst[ih]->Draw("E2sames");
   	   	   hMassMPbPbStat[ih]->Draw("sames");
   	   	   
   	   	   hMassPythia[ih]->Draw("sames");
   	   	   
   	   	   //debug
   	   	   //Printf("Drawn %s mean %f, %s mean %f, %s mean %f, %s mean %f, %s mean %f", hMassMSyst[ih]->GetName(), hMassMSyst[ih]->GetMean(), hMassMStat[ih]->GetName(), hMassMStat[ih]->GetMean(), hMassMPbPbSyst[ih]->GetName(), hMassMPbPbSyst[ih]->GetMean(), hMassMPbPbStat[ih]->GetName(), hMassMPbPbStat[ih]->GetMean(), hMassPythia[ih]->GetName(), hMassPythia[ih]->GetMean());
   	   	   
   	   	   if(ih == 0) {
   	   	   	   legRes->AddEntry(hMassMStat[ih], "Jet Mass no Fluc pPb", "LP");
   	   	   	   legRes->AddEntry(hMassMSyst[ih], "Tot Syst Unc no Fluc pPb", "F");
   	   	   	   legRes->AddEntry(hMassMPbPbStat[ih], "Jet Mass PbPb", "LP");
   	   	   	   legRes->AddEntry(hMassMPbPbSyst[ih], "Tot Syst Unc PbPb", "F");
   	   	   	   
   	   	   	   legRes->Draw();
   	   	   }
   	   	   
   	   	   
   	   }
   	   
   	   SaveCv(cResults);
   	   SaveCv(cSystPaper);
   }
   
   gSystem->ChangeDirectory(basewdircmp);
   
   Printf("Done");
   return;
}



