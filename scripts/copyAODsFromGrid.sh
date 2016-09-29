#!/bin/bash

path=/alice/sim/2013/LHC13b4_plus
AOD=AOD152
PTH=2

BASELOCDIR=/data/Work/ROOTFiles/LHC13b4_plus
cd $BASEDIR

runs=(195568) # LHC13b4_plus

outputfile=$BASELOCDIR/LHC13b4_plus_AODtest.txt
for k in `seq 0 0`;
do
	for j in `seq 1 9`;
	do
		cd $BASELOCDIR/${runs[$k]}/$PTH
		mkdir 000$j
		cd 000$j
		echo alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/000$j/AliAOD.root file:AnalysisResults.root
		#alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/000$j/pyxsec.root file:pyxsec.root
		alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/000$j/pyxsec_hists.root file:pyxsec_hists.root
		echo "$BASELOCDIR/${runs[$k]}/$PTH/000$j/AliAOD.root" >> $outputfile 
	done 
	for j in `seq 10 99`;
	do
		cd $BASELOCDIR/${runs[$k]}/$PTH
		mkdir 00$j
		cd 00$j
		echo alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/00$j/AliAOD.root file:AnalysisResults.root
		#alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/00$j/pyxsec.root file:pyxsec.root
		alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/00$j/pyxsec_hists.root file:pyxsec_hists.root
		echo "$BASELOCDIR/${runs[$k]}/$PTH/00$j/AliAOD.root" >> $outputfile 
	done
	
	for j in `seq 100 386`;
	do
		cd $BASELOCDIR/${runs[$k]}/$PTH
		mkdir 0$j
		cd 0$j
		echo alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/0$j/AliAOD.root file:AnalysisResults.root
		#alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/0$j/pyxsec.root file:pyxsec.root
		alien_cp alien://$path/${runs[$k]}/$PTH/$AOD/0$j/pyxsec_hists.root file:pyxsec_hists.root
		echo "$BASELOCDIR/${runs[$k]}/$PTH/0$j/AliAOD.root" >> $outputfile 
	done
done
cd $BASELOCDIR
