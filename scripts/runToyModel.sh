#!/bin/sh
# $1=hadTrig $2=dijet corr type (0=ChCh 1=FuCh)

#loop for all centralities
for i in `seq 1 1`;#6
do

    echo "processing R=0.2 RV=0 cent=$i hadTrig=$1 dijetCorrType=$2"
#    root -b -q "/data/wrk/postdoc/kt/Analysis/PlotMacros/runToyModelCorrectionsMultiple.C(\"../../AnalysisResultsWeighted.root\",\"/data/wrk/postdoc/kt/Analysis/Data/grid/LHC13c/Train237/AnalysisResults.root\",\"DiJetResponse_Jet_AKTFullR020_PicoTracks_pT0150_CaloClustersCorr_ET0300_Jet_AKTFullR020_MCParticlesSelected_pT0000_Rho0TCMatch3HadTrig$1RV0_histosMerged\",0.2,0,$i,$1,$2)"

    echo "processing R=0.2 RV=1 cent=$i hadTrig=$1 dijetCorrType=$2"
#    root -b -q "/data/wrk/postdoc/kt/Analysis/PlotMacros/runToyModelCorrectionsMultiple.C(\"../../AnalysisResultsWeighted.root\",\"/data/wrk/postdoc/kt/Analysis/Data/grid/LHC13c/Train237/AnalysisResults.root\",\"DiJetResponse_Jet_AKTFullR020_PicoTracks_pT0150_CaloClustersCorr_ET0300_Jet_AKTFullR020_MCParticlesSelected_pT0000_Rho0TCMatch3HadTrig$1RV1_histosMerged\",0.2,1,$i,$1,$2)"

    echo "processing R=0.4 RV=0 cent=$i hadTrig=$1 dijetCorrType=$2"
    root -b -q "/data/wrk/postdoc/kt/Analysis/PlotMacros/runToyModelCorrectionsMultiple.C(\"../../AnalysisResultsWeighted.root\",\"/data/wrk/postdoc/kt/Analysis/Data/grid/LHC13c/Train291/AnalysisResults.root\",\"DiJetResponse_Jet_AKTFullR040_PicoTracks_pT0150_CaloClustersCorr_ET0300_pt_scheme_Jet_AKTFullR040_MCParticlesSelected_pT0000_pt_scheme_Rho0TCMatch3HadTrig$1RV0_histosMerged\",0.4,0,$i,$1,$2)"

    #root -b -q "/data/wrk/postdoc/kt/Analysis/PlotMacros/runToyModelCorrectionsMultiple.C(\"../../DiJetResponse_Jet_AKTFullR040_PicoTracks_pT0150_CaloClustersCorr_ET0300_Jet_AKTFullR040_MCParticlesSelected_pT0000_Rho0TCMatch3HadTrig5RV0_histos.root\",\"/data/wrk/postdoc/kt/Analysis/Data/grid/LHC13c/Train237/AnalysisResults.root\",\"DiJetResponse_Jet_AKTFullR040_PicoTracks_pT0150_CaloClustersCorr_ET0300_Jet_AKTFullR040_MCParticlesSelected_pT0000_Rho0TCMatch3HadTrig$1RV0_histosMerged\",0.4,0,$i,$1,$2)"

    echo "processing R=0.4 RV=1 cent=$i hadTrig=$1 dijetCorrType=$2"
 #   root -b -q "/data/wrk/postdoc/kt/Analysis/PlotMacros/runToyModelCorrectionsMultiple.C(\"../../AnalysisResultsWeighted.root\",\"/data/wrk/postdoc/kt/Analysis/Data/grid/LHC13c/Train237/AnalysisResults.root\",\"DiJetResponse_Jet_AKTFullR040_PicoTracks_pT0150_CaloClustersCorr_ET0300_Jet_AKTFullR040_MCParticlesSelected_pT0000_Rho0TCMatch3HadTrig$1RV1_histosMerged\",0.4,1,$i,$1,$2)"

done


#ls DiJetResponseToy*.root > DiJetResponseToy.list

#root -b -q "PutTListInOneFile(\"DiJetResponseToy.list\",\"DiJetResponseToyAll.root\")"




