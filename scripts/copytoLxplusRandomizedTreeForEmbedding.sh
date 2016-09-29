#!/bin/bash
TRAINNUM=1134
BASEDIR=/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train$TRAINNUM/outputpTHardBins/mergeRuns
BASEOUTNAME=outRbis1

cd $BASEDIR
ls

scp "*/$BASEOUTNAME*.root" cbianchi@lxplus.cern.ch:/afs/cern.ch/user/c/cbianchi/CopyTmp/


