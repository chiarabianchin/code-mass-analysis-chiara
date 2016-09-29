#!/bin/bash
TRAINNUM=1261
BASEDIR=/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train$TRAINNUM/outputpTHardBins/mergeRuns
MINPTDCUT=-1.
MINPTPCUT=10.
MaxEntries=1500000
BASEOUTNAME=outR
for j in `seq 1 9`;
do
    cd $BASEDIR/0$j
    pwd
    root -b -q "/data/Work/MyCodeJetMass/macros/RandomizeTreeForEmbedding.C+(\"AnalysisResults.root\", \"fTreeJet\", 0, $j, $MINPTDCUT, $MINPTPCUT, $MaxEntries, \"$BASEOUTNAME\" )"
    ls
    alien_cp "$BASEOUTNAME$j.root" alien:///alice/cern.ch/user/c/cbianchi/FileForEmbedding/$BASEOUTNAME$TRAINNUM$j.root
done
j=10
cd $BASEDIR/$j
pwd
root -b -q "/data/Work/MyCodeJetMass/macros/RandomizeTreeForEmbedding.C+(\"AnalysisResults.root\", \"fTreeJet\", 0, $j , $MINPTDCUT, $MINPTPCUT, $MaxEntries, \"$BASEOUTNAME\")"
ls
alien_cp "$BASEOUTNAME$j.root" alien:///alice/cern.ch/user/c/cbianchi/FileForEmbedding/$BASEOUTNAME$TRAINNUM$j.root

cd $BASEDIR

