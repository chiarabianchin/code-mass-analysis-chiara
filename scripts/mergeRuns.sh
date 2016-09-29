#!/bin/bash

BASEDIR=`pwd`

[[ ! -d $BASEDIR/mergeRuns  ]] && mkdir $BASEDIR/mergeRuns

for bin in `seq 1 9`;
do
    [[ ! -d $BASEDIR/mergeRuns/0$bin  ]] && mkdir $BASEDIR/mergeRuns/0$bin

    cd $BASEDIR/mergeRuns/0$bin

    ls $BASEDIR/*/0$bin/AnalysisResults.root > output.list

    root -b -q "/data/macros/MergeFilesLocal.C(\"output.list\",\"AnalysisResults.root\")"

    cd $BASEDIR/mergeRuns
    ls $BASEDIR/mergeRuns/0$bin/AnalysisResults.root >> mergedoutput.txt
    cd $BASEDIR

done

[[ ! -d $BASEDIR/mergeRuns/10  ]] && mkdir $BASEDIR/mergeRuns/10

cd $BASEDIR/mergeRuns/10

ls $BASEDIR/*/10/AnalysisResults.root > output.list

root -b -q "/data/macros/MergeFilesLocal.C(\"output.list\",\"AnalysisResults.root\")"

cd $BASEDIR/mergeRuns
ls $BASEDIR/mergeRuns/10/AnalysisResults.root >> mergedoutput.txt
cd $BASEDIR
