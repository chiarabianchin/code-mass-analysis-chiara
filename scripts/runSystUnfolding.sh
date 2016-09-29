#!/bin/sh
# unfolding systematics

# switches
DOMATRIX=0
DOUNFOLD=1
DOCOMPARE=1

# working directories
wdirresp=("/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/MB/redo" "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/MB/redo" "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/EJE/redo" "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/EJE/redo")

# settings response matrix
strIn=("/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbDeriv/ResponseWJetShapeDeriv_JetRhosub_AKTChargedR040_PicoTracks.root" 
"/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/UnfoldEmbConst/ResponseWJetShapeConst_JetRhosub_AKTChargedR040_PicoTracks.root") #(response Deriv response Const)
strL="fhResponseFinal"
tag=("DetFlBkgDeriv" "DetFlBkgConst")
#"DetFlBkgDerivptpar20"
#"DetFlBkgConstptpar20"
pt_min=-40.
m_min=-20.
pt_minT=0.
m_minT=0.
# in the response this setting is not used because the ranges are set from outside
colTypefixed=-1

binWidthPt=5.
skipBins=0

# settings unfolding
dataIn=("/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/AnalysisResultsMB.root" "/data/Work/jets/JetMass/pPbJetMassAnalysis/ResultspPbJetMass/Train806-807-810-811/AnalysisResultsEJE.root")
# the name structure is: $respFileNameBase${tag[$i]}.root, the path is ${wdirresp[$i]}/
respFileNameBase="response"
itermin=1
itermax=5
iterdef=3
# 0 = deriv, 1 = const sub
BkgType=(0 1)
# 1 = MB, 2 = EJE
ColType=(1 2)

# settings comparison
