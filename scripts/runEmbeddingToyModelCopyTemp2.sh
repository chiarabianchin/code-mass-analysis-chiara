#!/bin/sh

#simulation and merging can be run togehter or separately setting the following switches
declare -i SIMU=1
declare -i MERGE=1

# indicate the working directory
BASEDIR=/data/Work/jets/JetMass/BkgFluctStudies
SUBDIR=CleanEmbedding
EMBTYPEDIR=SingleMassive5GeVTrackInThermal20k
#AODBASEDIR=/data/Work/JEQA/LHC15g6/outputpTHardBins/d/Train296/mergeRuns
#SLIDEDIR=/data/Documents/Presentazioni/PWGJE/QAchecks/JETQA/LHC15g6LHC10pass4/pictures/TrackOnly/LHC15g6/d/Train457

# set the properties of the simulation 
PYTHIAJET=0   # jets from PYTHIA
THERMAL=1     # underlying particle distribution with average pT and multiplicity given by the setting in the run macro
EMBCODE=0     # merge the two above (2), nothing (0), as of now it runs thermal -not so clear (1)
SINGLEPART=1  # embed a single high pT track as in macro

JETRECSCHEME=0 # recombination scheme chose in jet reconstruction (0 = E-scheme, 1 = pt-scheme)
FRACSHARE=0.2  # fraction of PYTHIA pT shared by the reconstructed jet
PTHARDMin=100  # minimum pT hard bin for PYTHIA simulation
PTHARDMax=500 # maximum pT hard bin for PYTHIA simulation
PTPARTMin=40   # minimum pT single part embedding
PTPARTMax=120  # maximum pT single part embedding
MASS=5     # single track is massless = -1
outputfile=$BASEDIR/$SUBDIR/$EMBTYPEDIR/inputfiles.txt

cd $BASEDIR
[[ ! -d $SUBDIR ]] && mkdir $SUBDIR
cd $SUBDIR
[[ ! -d $EMBTYPEDIR ]] && mkdir $EMBTYPEDIR
cd $EMBTYPEDIR

# save setting into settingsscript.dat (add by hand additional info)
FILESETTINGS=settingsscript.dat
echo "PYTHIAJET=$PYTHIAJET"  >> $FILESETTINGS
echo "THERMAL=$THERMAL"  >> $FILESETTINGS
echo "EMBCODE=$EMBCODE"  >> $FILESETTINGS
echo "SINGLEPART=$SINGLEPART"  >> $FILESETTINGS
echo "JETRECSCHEME=$JETRECSCHEME"  >> $FILESETTINGS
echo "FRACSHARE=$FRACSHARE"  >> $FILESETTINGS
echo "PTHARDMin=$PTHARDMin"  >> $FILESETTINGS
echo "PTHARDMax=$PTHARDMax"  >> $FILESETTINGS
echo "PTPARTMin=$PTPARTMin"  >> $FILESETTINGS
echo "PTPARTMax=$PTPARTMax"  >> $FILESETTINGS

if [ "$SIMU" -eq 1 ]; then
#run simulation

for bin in `seq 1 7`;
do
   
   [[ ! -d 0$bin  ]] && mkdir 0$bin
   cd 0$bin
   
   pwd
   bininput=$(( bin - 1 ))
   echo "Bin of multiplicity $bininput"
   if [ $SINGLEPART -eq 1 ]; then
      root -b -q "$BASEDIR/$SUBDIR/runTestJetMassDerivWithScript.C($bininput, $PYTHIAJET, $THERMAL, $EMBCODE, $SINGLEPART, $JETRECSCHEME, $FRACSHARE, $PTHARDMin, $PTHARDMax, $PTPARTMin, $PTPARTMax, \"ESD\", \"local\", $MASS)"
   else # need to use a different macro when running without single track embedding (temporary?)
      root -b -q "$BASEDIR/$SUBDIR/JetPYTHIAPlusThermal/runTestJetMassDerivWithScript.C($bininput, $JETRECSCHEME, $FRACSHARE, $PTHARDMin, $PTHARDMax)"
   fi
   cd $BASEDIR/$SUBDIR/$EMBTYPEDIR
   echo "$BASEDIR/$SUBDIR/$EMBTYPEDIR/0$bin/AnalysisResults.root" >> $outputfile
done

fi

#merge results
if [ $MERGE = 1 ]; then
   cd $BASEDIR/$SUBDIR/$EMBTYPEDIR
   pwd
   #create text files with list names and file names for the two steps of merging
   root -b -q "/data/Work/MyCodeJetMass/macros/ExtractTextFileList.C+()"
   
   #create files with each list merged
   
   . /data/Work/MyCodeJetMass/scripts/WeightingToyModelOutput.sh $EMBTYPEDIR
   
   #merge the lists in one file
   cd merge
   root -b -q "/data/Work/CodeJetMass/PtHardUtil/PutTListInOneFile.C(\"../listOfFilesToMerge.txt\")"
   
fi
