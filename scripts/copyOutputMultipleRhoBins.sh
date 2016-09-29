#!/bin/bash

curdir=`pwd`
echo "Directory: $curdir"
nmaster=2
nruns=8
basealien="alien:///alice/cern.ch/user/c/cbianchi/work"
terminate=("Emb1134HRhoBinsOvlExcl_20160708_143453" "Emb1134HRhoBinsOvlExcl_20160708_143604")
listtxt="filestomerge.txt"
for j in `seq 0 $nmaster`;
do
	cd 0$j
	echo `pwd`
	for k in `seq 0 $nruns`
	do
		if [ "$k" -lt "10" ] 
		then 
			alien_cp $basealien/${terminate[$j]}/output/00$k/AnalysisResults.root AnalysisResults$k.root
		fi
		if [ "$k" -gt "9" ]
		then 
			alien_cp $basealien/${terminate[$j]}/output/0$k/AnalysisResults.root AnalysisResults$k.root
		fi
		echo AnalysisResults$k.root >> $listtxt
	done
	root -b -q -x /data/Work/MyCodeJetMass/macros/TFileMergeNOutputs.C+(4,\"$listtxt\",\"AnalysisResults.root\")
	cd $curdir
done

