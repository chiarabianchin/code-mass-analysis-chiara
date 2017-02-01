#!/bin/bash
cdir=`pwd`

echo "Working on dir $cdir"
ls
basename="JetMassALICE"
listofentries=("MasspPb" "MassPbPb" "RatioPbPbpPb") 
listofptbins=("Pt60_80" "Pt80_100" "Pt100_120")
listmean=("MeanMasspPb" "MeanMassPbPb")
filenameallyaml=listofyamlfiles.txt
first="no"

c++ -o txt2yaml /data/Work/code-mass-analysis-chiara/utils/txt2yaml.cpp

for i in "${listofentries[@]}";
do 
	echo $i
	for j in "${listofptbins[@]}";
	do
		echo $j
		filename="$basename"_"$i"_"$j"
		filetxt=$filename.txt
		fileyaml=$filename.yaml
		./txt2yaml $filetxt $fileyaml
		if [ $first == "no" ];
		then 
			echo "$fileyaml" > "$filenameallyaml"
			echo $first
			first="yes"
			echo $first
		else echo "$fileyaml" >> "$filenameallyaml"
		fi
	done
done

for i in "${listmean[@]}";
do
	echo $i
	filename="$basename"_"$i"
	filetxt=$filename.txt
	fileyaml=$filename.yaml
	./txt2yaml $filetxt $fileyaml
	echo "$fileyaml" >> "$filenameallyaml"
done
	
echo "Written file $filenameallyaml"
