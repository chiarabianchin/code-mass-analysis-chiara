#!/bin/bash
TRAINNUM=1134
BASEDIR=/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train$TRAINNUM/outputpTHardBins/mergeRuns
BASEOUTNAME=outRbis1
for j in `seq 1 9`;
do
    cd $BASEDIR/0$j
    pwd
    ls
    alien_cp "$BASEOUTNAME$j.root" alien:///alice/cern.ch/user/c/cbianchi/FileForEmbedding/$BASEOUTNAME$TRAINNUM$j.root
done
j=10
cd $BASEDIR/$j
pwd
ls
alien_cp "$BASEOUTNAME$j.root" alien:///alice/cern.ch/user/c/cbianchi/FileForEmbedding/$BASEOUTNAME$TRAINNUM$j.root

cd $BASEDIR

