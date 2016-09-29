#!/bin/sh 
#use first argument if xsec and trials are stored in a different list than the one that is being merged

export filename=listNames.txt
export pytname=$1

echo "$pytname"

while read line; do
    echo "$line"
    root -b -q "/data/Work/CodeJetMass/PtHardUtil/MergePtHardSerial.C(\"output.list\",\"$line\",\"$1\")"
done < "$filename"


#use $gitJetMass/PtHardUtil/PutTListInOneFile to get all merged lists in one file
#need to create txt file with all root files first (all those produced in do-loop above)
