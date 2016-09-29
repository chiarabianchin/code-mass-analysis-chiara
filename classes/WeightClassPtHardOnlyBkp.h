/// \class WeightClass
/// \brief class for pT-hard bin weighting
///
/// \author Chiara Bianchin
/// \date January 2015

#ifndef WEIGHTCLASS_H
#define WEIGHTCLASS_H

class THnSparse;
class TFile;
class TList;
class TH1;
class TH2;
class TH1D;
class TString;

class WeightClass : public TNamed {
   
public:
   
   WeightClass();
   WeightClass(const char* name, const char* title, const Int_t npthbins);
   ~WeightClass() {;}
   
   // setters & getters
  
   
   void SetFirstPtHard(Int_t first)                        { fFirstpthB = first;}
   void SetResponseName(TString name)                      { fResponsethnsp = name;}
   void SetDoNormalisation(Bool_t doN = kTRUE)             { fDoNormPerInt = doN; if(fIs2D) Printf("Not implemented for 2D");}
   void SetDoNormalisationAndReadNJet(Bool_t doN = kTRUE, TString filename = "NJetsPerPtHBin.root")             { SetDoNormalisation(doN); ReadNumberOfJetPerPtHBin(filename);}
   void SetIs2D(Bool_t is2d = kTRUE, Int_t axpar = 3, Int_t axdet = 2) { fIs2D = is2d; fAxisPtp = axpar; fAxisPtd = axdet;  }
   void SetOrderAxesRespose(Int_t axMp, Int_t axMd, Int_t axptp, Int_t axptd) { fAxisMd = axMd; fAxisMp = axMp; fAxisPtd = axptd; fAxisPtp = axptp; SetArrayOrderAxesResponse(); Printf("Order of axes set to :  axMp = %d, axMd = %d, axptp = %d, axptd = %d", fAxisMp, fAxisMd, fAxisPtp, fAxisPtd); }
   void SetOrderAxesMPt3D(Int_t axMp, Int_t axMd, Int_t axptp, Int_t axptd, Int_t axMdE, Int_t axptdE) { fAxisMd = axMd; fAxisMp = axMp; fAxisPtd = axptd; fAxisPtp = axptp; fAxisMdEmb = axMdE; fAxisPtdEmb = axptdE; Printf("SetAxesWithArrayMPt3D(). Axes set to :  axMp = %d, axMd = %d, axptp = %d, axptd = %d, axMdEmb = %d, axptdEmb = %d", fAxisMp, fAxisMd, fAxisPtp, fAxisPtd, fAxisMdEmb, fAxisPtdEmb); return;}
   
   void SetCut(Int_t axcut, Double_t cutlow = -999, Double_t cuthigh = -999)  { fCutAxis = axcut; fMinCut = cutlow; fMaxCut = cuthigh; }
   void  SetFineFinning(Bool_t set = kTRUE)            { fFineBinningSlow = set;}
   
   THnSparseF*  GetResponsePtHbin(Int_t ipthb) const       { return fhResponse[ipthb]; }
   THnSparseF** GetResponse()const                         { return fhResponse; }
   THnSparseF*  GetMPt3DPtHbin(Int_t ipthb) const          { return fhMPt3D[ipthb];}  
   THnSparseF** GetMPt3D() const                           { return fhMPt3D;}
   THnSparseF*  GetResponseWPtHbin(Int_t ipthb)const       { return fhResponseW[ipthb]; }
   THnSparseF** GetResponseW() const                       { return fhResponseW; }
   THnSparseF*  GetMPt3DWPtHbin(Int_t ipthb) const         { return fhMPt3DW[ipthb];}  
   THnSparseF** GetMPt3DW() const                          { return fhMPt3DW;}
   THnSparseF** GetResponseWNorm() const                   { return fhResponseN; }
   THnSparseF*  GetResponseWNormPtHbin(Int_t ipthb)const   { return fhResponseN[ipthb]; }
   
   THnSparseF*  GetResponseWFinal() const                  { return fhResponseWFinal; }
   THnSparseF*  GetMPt3DWFinal() const                     { return fhMPt3DWFinal;}
   TH2*         Get2DPartProjPtHbin(Int_t ipthb) const     { return fhptMpar[ipthb]; }
   TH2**        Get2DPartProj() const                      { return fhptMpar; }
   TH2**        Get2DPartProjW() const                     { return fhptMparW; }
   TH2*         Get2DPartProjWPtHbin(Int_t ipthb) const    { return fhptMparW[ipthb]; }
   TH1*         GetCrossSec() const                        { return fhCrossSec ;}
   TH1*         GetNTrials()   const                       { return fhTrials   ;}
   TH1*         GetScaleF()    const                       { return fhScaleF   ;}
   TH1*         GetNormInt()   const                       { return fhNormInt  ;}
   Int_t        GetMParAxis()  const                       { return fAxisMp; }
   Int_t        GetMDetAxis()  const                       { return fAxisMd; }
   Int_t        GetPtParAxis() const                       { return fAxisPtp; }
   Int_t        GetPtDetAxis() const                       { return fAxisPtd; }
   TCanvas*     DrawHistogramTH1(TH1* histogram, Bool_t log = kFALSE) const;
   TCanvas*     DrawCrossSec() const                       { return DrawHistogramTH1(fhCrossSec, kTRUE) ;}
   TCanvas*     DrawNTrials()  const                       { return DrawHistogramTH1(fhTrials  , kTRUE) ;}
   TCanvas*     DrawScaleF()   const                       { return DrawHistogramTH1(fhScaleF  , kTRUE) ;}
   TCanvas*     DrawNormInt()  const                       { return DrawHistogramTH1(fhNormInt , kTRUE) ;}
  
   
   TH1D**       All1DResponseProj(Int_t ipthb, Int_t color = 0, TString nameh = "") const  { Printf("Projecting response"); return All1DProjections(fhResponse[ipthb], color, nameh); }
   TH1D**       All1DResponseWProj(Int_t ipthb, Int_t color = 0, TString nameh = "") const { Printf("Projecting weighted  response");return All1DProjections(fhResponseW[ipthb], color, nameh); }
   TH1D**       All1DResponseWNormProj(Int_t ipthb, Int_t color = 0, TString nameh = "") const         { Printf("Projecting weighted normalised response"); return All1DProjections(fhResponseN[ipthb], color, nameh); }
   TH1D**       All1DResponseWFinalProj(Int_t color = 0, TString nameh = "") const         { return All1DProjections(fhResponseWFinal, color, nameh); }
   
   TCanvas*     CompareResponse(THnSparseF *response2) const;
   TCanvas*     CompareResponses(THnSparseF *response1, THnSparseF *response2, Int_t color1 = 1, Int_t color2 = 2, TString leg1 = "My response", TString leg2 = "Reference resp") const;
   void         CleanPointerResponseFinal()                         { fhResponseWFinal = 0; Printf("Clean Response final -> %p", fhResponseWFinal); }
   
   // main utilities
   void   Initialise(TString path, TString list, TString inputFile = "", TString listnameInp = "");
   Bool_t ReadWeight(Int_t ipthb);
   Bool_t ReadInputDataFromList(Int_t ipthb);
   void   ReadInputData(Int_t ipthb, THnSparseF* hspresp);
   Bool_t SetResponse(Int_t ipthb, TString hnspname);
   Bool_t WeightResponse(Int_t ipthb);
   void   NormalizePerIntegral(Int_t ipthb);
   void   NormalizePerBin(Int_t ipthb);
   Bool_t Result();
   
   
   // these methods are for the 2D response (pT only). It's used for cross-check
   void    ReadInputData2D(Int_t ipthb, TH2* h2dresp);
   Bool_t  SetResponse2D(Int_t ipthb);
   Bool_t  WeightResponse2D(Int_t ipthb);
   TH2**   GetResponse2DW() const                        { return fh2dResponseW; }
   TH2*    GetResponse2DWPtHbin(Int_t ipthb) const       { return fh2dResponseW[ipthb]; }
   TH2*    GetResponse2DWFinal() const                   { return fh2DRespFinal; }
   TH1*    GetResponse2DProjPar(Int_t ipthb) const       { return GetProjPar(fh2dResponse[ipthb], "", ipthb); }
   TH1*    GetResponse2DProjDet(Int_t ipthb) const       { return GetProjDet(fh2dResponse[ipthb], "", ipthb); }
   TH1*    GetResponse2DWProjPar(Int_t ipthb) const      { return GetProjPar(fh2dResponseW[ipthb], "", ipthb); }
   TH1*    GetResponse2DWProjDet(Int_t ipthb) const      { return GetProjDet(fh2dResponseW[ipthb], "", ipthb); }
   TH1*    GetResponse2DWFinalProjPar() const            { return GetProjPar(fh2DRespFinal, "", 1); }
   TH1*    GetResponse2DWFinalProjDet() const            { return GetProjDet(fh2DRespFinal, "", 2); }
   
   TCanvas* DrawResponseWProjPar(Int_t ipthb) const      { return DrawHistogramTH1(GetResponse2DWProjPar(ipthb)); }
   TCanvas* DrawResponseWProjDet(Int_t ipthb) const      { return DrawHistogramTH1(GetResponse2DWProjDet(ipthb)); }

   TCanvas* DrawResponse2DWFinalProjPar() const          { return DrawHistogramTH1(GetResponse2DWFinalProjPar()); }
   TCanvas* DrawRespons2DeWFinalProjDet() const          { return DrawHistogramTH1(GetResponse2DWFinalProjDet()); }
   
   
   //more basic methods
   
   TH1D**       All1DProjections(THnSparseF *hnsp, Int_t color = 0, TString nameh = "") const;
   TList*       ReadFile(TString strIn, TString strLst); //to be removed when include is available
   void         UniformTH1FForDivide(TH1 *h1, TH1 *h2, TH1 *&h1bis, TH1 *&h2bis, TString type) const;
   TH1*         GetProjPar(TH2* h2, TString name = "", Int_t color = 1) const { return GetProj(h2, fAxis2DPtp, name, color); }
   TH1*         GetProjDet(TH2* h2, TString name = "", Int_t color = 1) const { return GetProj(h2, fAxis2DPtd, name, color); }
   TH1*         GetProj(TH2* h2, Int_t ax, TString name = "", Int_t color = 1) const;
   
   // save the number of jets per pT hard bin in the MC jobs
   void         SaveNumberOfJetPerPtHBin(TString filename = "NJetsPerPtHBin.root", Int_t level = 1) const;
   void         ReadNumberOfJetPerPtHBin(TString filename = "NJetsPerPtHBin.root");
   
protected:
   
   void         SetArrayOrderAxesResponse()             { fAxesResp[0] = fAxisMd; fAxesResp[1] = fAxisMp; fAxesResp[2] = fAxisPtd; fAxesResp[3] = fAxisPtp; }
   void         SetArrayOrderAxesMPt3D()                { fAxesMPt3D[0] = fAxisMd; fAxesMPt3D[1] = fAxisMp; fAxesMPt3D[2] = fAxisPtd; fAxesMPt3D[3] = fAxisPtp; fAxesMPt3D[4] = fAxisMdEmb; fAxesMPt3D[5] = fAxisPtdEmb;}
   
   void         SetAxesMPt3DWithArray()                 { fAxisMd = fAxesMPt3D[0]; fAxisMp = fAxesMPt3D[1]; fAxisPtd = fAxesMPt3D[2]; fAxisPtp = fAxesMPt3D[3]; fAxisMdEmb = fAxesMPt3D[4]; fAxisPtdEmb = fAxesMPt3D[5];}
   
   TH1D*        IntegralParPt(Int_t ipthb, Int_t level, Double_t integralptpar) const;
   
   Int_t           CheckBinning(Int_t ipthb, TH1D* h) const;
   
   TString fpathPythiaFile;                          ///< local path where to look for the file(s) with weights
   TString flistPythiaFile;                          ///< name of the list where to find the input for cross sec and trials
   TString fnamehtrials    ;                         ///< name of the histogram
   TString fnamehxsec      ;                         ///< name of the histogram
   TString fnamehevents    ;                         ///< name of the histogram
   TString fResponsethnsp  ;                         ///< name of the histogram
   TString fDMdPtthnsp  ;                         ///< name of the other histogram
   
   
   const Int_t fNpthB;                        ///< Number of pT hard bins
   Int_t fFirstpthB;                                 ///< first (of the fNpthB) pt hard bin to be considered
   Int_t fAxesResp[4];                               ///< Axis numbers for Mass_det, Mass_par, Pt_det, Pt_par
   Int_t fAxesMPt3D[4];                              ///< Axis numbers for Mass_det, Mass_par, Pt_det, Pt_par, Mass_detFluc, Pt_detFluc
   Int_t fAxisMp  ;                                  ///< Axis number for Mass_par
   Int_t fAxisMd  ;                                  ///< Axis number for Mass_det (plus fluctuation if embedding)
   Int_t fAxisPtp ;                                  ///< Axis number for Pt_par
   Int_t fAxisPtd ;                                  ///< Axis number for Pt_det (plus fluctuation if embedding)
   Int_t fAxis2DPtp ;                                ///< Axis number for Pt_par in 2D distribution
   Int_t fAxis2DPtd ;                                ///< Axis number for Pt_det in 2D distribution
   Int_t fAxisMdEmb;                                 ///< Axis number for Mass_det used for embedding (defined only if embedding)
   Int_t fAxisPtdEmb;                                ///< Axis number for Pt_det used for embedding (defined only if embedding)
   
   Float_t fxsec;                                    //!<! Current value of cross section
   Float_t ftrials;                                  //!<! Current value of n trials
   Float_t fnevts;                                   //!<! Current value of number of events
   
   TString fInputFileName;                           ///< Input file name
   TString fInputListName;                           ///< Input list name
   TList *fInputList;                                //!<! TList containing the THnSparse of the response from embedding
   
   TH1 *fhCrossSec ;                                 ///< Store the input x-sec  used
   TH1 *fhTrials   ;                                 ///< Store the input trials used
   TH1 *fhScaleF   ;                                 ///< Store the factors applied 
   TH1 *fhNormInt  ;                                 ///< Store the factors applied 

   //raw
   ///< 2D response (e.g. only pT)
   TH2          **fhptMpar   ;                        //[fNpthB]
   ///< Response per jet
   THnSparseF   **fhResponse ;                        //[fNpthB]
   ///< Another THnSparse to be weighted, e.g. variation of M and pT
   THnSparseF   **fhDMdPt;                            //[fNpthB]
   ///< Response per jet saving also the Detector level embedded particle M,pT (used only with embedding)
   THnSparseF   **fhMPt3D ;                           //[fNpthB]
  
   
   //weighted
   /// Response weighted by sigma/Ntrials
   THnSparseF   **fhResponseW;                        //[fNpthB]
   /// Response normalized per jet saving also the Detector level embedded particle M,pT (used only with embedding)   
   THnSparseF   **fhMPt3DW   ;                        //[fNpthB]

   /// Another THnSparse- e.g. variation of M and pT, weighted by sigma/Ntrials
   THnSparseF    **fhDMdPtW;                            //[fNpthB]
   
   /// projection of weighted response fhResponseWFinal on the particle level quantities
   TH2          **fhptMparW  ;                        //[fNpthB] 
   /// projection of weighted response fhResponseWFinal on the particle detector quantities
   TH2          **fhptMdetW  ;                        //[fNpthB]
   
   // normalized
   /// Response per jet weighted by sigma/Ntrials and normalized, per bin, by the fraction of total number of jets over those selected. Further inputs required
   THnSparseF   **fhResponseN;                        //[fNpthB]
   /// Another THnSparse- e.g. variation of M and pT, weighted by sigma/Ntrials and normalized, per bin, by the fraction of total number of jets over those selected. Further inputs required
   THnSparseF   **fhDMdPtN;                            //[fNpthB]
   // projections
   TH2          **fhptMparN  ;                        //[fNpthB] 
   TH2          **fhptMdetN  ;                        //[fNpthB]
   //final
   THnSparseF   *fhResponseWFinal;                    ///< Final response after weighting and summing
   THnSparseF   *fhMPt3DWFinal  ;                    ///< Weighted response containing also the Detector level embedded particle M,pT (used only with embedding)
   THnSparseF   *fhDMdPtWFinal;                    ///< other histogram after weighting and summing
      
   Bool_t        fDoNormPerInt;                       ///< Set whether to apply further normalisations
   
   // these data members are for the 2D response (pT only). It's used for cross-check
   /// The 2D weighting is not implemented for the response containing also the Detector level embedded particle M,pT (used only with embedding)
   Bool_t         fIs2D;                              ///< flag for 2D response
   ///< Response in 2D in pT-hard bins
   TH2          **fh2dResponse;                        //[fNpthB]
   ///< Wighted response in 2D in pT-hard bins        
   TH2          **fh2dResponseW;                       //[fNpthB]
   
   TH2           *fh2DRespFinal;                       ///< Final 2D response
   
   // if requested, apply a cut on one axis
   Int_t         fCutAxis;                             ///< Apply cut on the initial response to this axis (done before reshuffling the axes according to fAxesResp)
   Double_t      fMinCut;                              ///< Lower lim of the cut, applied only if != -999
   Double_t      fMaxCut;                              ///< Upper lim of the cut, applied only if != -999
   
   // weighting inputs
   TList        *fListNJets;                           ///< list of TParameter that contains the number of jets at particle level in the MC jobs and an histogram with the content of each pT,par bin, per each pT-hard bin
   Bool_t      fFineBinningSlow;               ///< Use the same binning as the original response for the mass. If kFALSE, merge 2 bins (faster!!). It is used in the method NormalizePerBin, if fDoNormPerInt = kTRUE
private:
   WeightClass(const WeightClass &weightcl);
   WeightClass& operator=(const WeightClass &weightcl);

   ClassDef(WeightClass, 1)
};
#endif
