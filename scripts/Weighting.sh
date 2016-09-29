#!/bin/sh
export fileinput=listNamesAll.txt
#listNames.txt
#\"JetTagger_Jet_AKTChargedR040_PicoTracks_pT0150_E_scheme_Jet_AKTChargedR040_MCParticlesSelected_pT0000_E_scheme_TC\"

while read line; do

   echo $line
   root -b -q "/data/Work/CodeJetMass/PtHardUtil/MergePtHardSerial.C(\"mergedoutput.txt\",\"$line\", \"JetTagger_Jet_AKTChargedR040_tracks_pT0150_E_scheme_Jet_AKTChargedR040_mcparticles_pT0000_E_scheme_TC\")"
   #root -b -q "/data/Work/MyCodeJetMass/macros/MergePtHardSerial.C(\"mergedoutput.txt\",\"$line\")"
done <$fileinput
