#!/bin/bash

#path=/alice/sim/2012/LHC12a15e_fix
path=/alice/sim/2013/LHC13b4_plus
Trainnumb=1261
Time=_20160805-1053/
train=/PWGJE/Jets_EMC_pPb/$Trainnumb$Time
AOD=AOD152

#BASEDIR=/data/Work/jets/JetMass/DetectorCorrections/LHC12a15e/Train415/outputpTHardBins #`pwd`
BASEDIR=/data/Work/jets/JetMass/DetectorCorrections/LHC13b4_plus/Train$Trainnumb/outputpTHardBins
mkdir $BASEDIR
#runs=(169838 170040) #good  LHC12a15e
#runs=(170040) #semi good
runs=(195531 195567 195568 195593 195644 195673 195783 195831 196085 196310) # LHC13b4_plus
#for j in `seq 0 1`;

outputfile=$BASEDIR/filepaths.txt
for j in `seq 0 9`;
do
    cd $BASEDIR
    mkdir ${runs[$j]}

    for bin in `seq 1 9`;
    do
       echo "Run $j is ${runs[$j]}"
	[[ ! -d $BASEDIR/${runs[$j]}/0$bin  ]] && mkdir $BASEDIR/${runs[$j]}/0$bin
	cd $BASEDIR/${runs[$j]}/0$bin
	 alien_cp alien://$path/${runs[$j]}/$bin/$AOD/$train/AnalysisResults.root AnalysisResults.root
	 echo "$BASEDIR/${runs[$j]}/0$bin/AnalysisResults.root" >> $outputfile 
    done
    
    [[ ! -d $BASEDIR/${runs[$j]}/10  ]] && mkdir $BASEDIR/${runs[$j]}/10
    cd $BASEDIR/${runs[$j]}/10
     alien_cp alien://$path/${runs[$j]}/10/$AOD/$train/AnalysisResults.root AnalysisResults.root
     echo "$BASEDIR/${runs[$j]}/10/AnalysisResults.root" >> $outputfile
done

    cd $BASEDIR
