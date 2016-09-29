#!/bin/bash
TRAINNUM=1134
BASEDIR=/afs/cern.ch/user/c/cbianchi/CopyTmp/
BASEOUTNAME=outRbis1
for j in `seq 1 10`;
do
    pwd
    ls
    alien_cp "$BASEOUTNAME$j.root" alien:///alice/cern.ch/user/c/cbianchi/FileForEmbedding/$BASEOUTNAME$TRAINNUM$j.root
done

