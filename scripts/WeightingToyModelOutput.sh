#!/bin/sh

# Create the file listOfListNamesToMerge.txt with /data/Work/MyCodeJetMass/macros/ExtractTextFileList.C
# 

#file with the lists
export fileinput=listOfListNamesToMerge.txt
#file with the paths of the files to be merged
BASEDIR=/data/Work/jets/JetMass/BkgFluctStudies/CleanEmbedding
DIR=$1 # e.g.SingleTrackInThermal, the subdirectory of the embedding exercise
FILENAME=inputfiles.txt

cd $BASEDIR/$DIR
# this is done in the root macro
#[[ ! -d merge ]] && mkdir merge
#cd merge

while read line; do

   echo $line
   root -b -q "/data/Work/MyCodeJetMass/macros/MergeEmbeddingThrmModelBins.C+(\"$line\", \"/data/Work/jets/JetMass/BkgFluctStudies/CleanEmbedding/WeightFactors.root\", \"gfractionsvsNch\", \"$BASEDIR/$DIR/$FILENAME\")"

done <$fileinput

