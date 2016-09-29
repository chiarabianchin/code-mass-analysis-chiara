#!/bin/bash
CURDIR=/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns
FILEEMB=/data/Work/jets/JetMass/pPbJetMassAnalysis/Embedding20160325/output/merge/AnalysisResults.root
LISTEMB=SingleTrackEmbedding
LISTRESa=JetShapeConst_JetRhosub_AKTChargedR040_PicoTracksPtH
LISTRESb=_pT0150_E_scheme_TC
FILETRE=/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train1134/outputpTHardBins/mergeRuns/
NAMETRE=outRbis
EMBPTCUT=20.
SAVETRE=1

for j in `seq 3 9`;
do
   mkdir $CURDIR/0$j
   cd $CURDIR/0$j
   pwd
   root -b -q "/data/Work/MyCodeJetMass/macros/RunDrawEntriesUsed.C+( \"$FILEEMB\", \"$LISTEMB$j\", \"$LISTRESa$j$LISTRESb\", \"$FILETRE/0$j/$NAMETRE$j.root\", $EMBPTCUT, $SAVETRE, $j)"
   
done
#bin 10
lastbin=10
mkdir $CURDIR/$lastbin
cd $CURDIR/$lastbin
pwd
root -b -q "/data/Work/MyCodeJetMass/macros/RunDrawEntriesUsed.C+( \"$FILEEMB\", \"$LISTEMB$lastbin\", \"$LISTRESa$lastbin$LISTRESb\", \"$FILETRE/$lastbin/$NAMETRE$lastbin.root\", $EMBPTCUT, $SAVETRE, $lastbin)"
   
cd $CURDIR

