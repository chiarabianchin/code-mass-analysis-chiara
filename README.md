# code-mass-analysis-chiara
Code developped for the Jet Mass anlaysis in pPb collision, originally based on the PbPb code by Marta

# Content of the project
The project contains mainly macros for analysis and plotting. Some common code is contained in utils/CommonTools.C that has to be included in most of the macros. Another commonly used set of functions is in LoadALICEFigures.C.

In classes one can find the class for pT hard bin weightig, to be used together with the macro macros/WeightResponse.C

In  the directory unfold, one can find the macros for creating the response matrix from the input THnSparse unfold/CreateRooUnfoldResponse.C and perform the unfolding with unfold/unfoldData.C. 
The directory unfold contains also the macros for the systematic uncertainties estimation. unfold/SystematicComparisonsUnf.C contains the basic methods, while unfold/runSystUnfolding.C takes care of running all the variations (check the input directory names first). The macros unfold/runSystUnfoldingVarbW.C is dedicated to the systematic variations using the variable bin width unfolding.

The directory scripts contains scripts for copying sets of files from grid (e.g. copyFromGrid.sh), randomize the order of tree entries (runRandomizeTreeForEmbedding.sh), copy the (randomized) input tree for embedding to alien (copyRandomizedTreeForEmbedding.sh), a set of scripts for running the embedding on toy model background.
The macros used inside those scripts are in general in the directory macros.

The methos macros/DrawComparisonsForPaper.C draws the figures for the paper (those related to pPb and some PbPb ones).

# Practical note
The includes of the macros have to be changed according to the local path
In the macro macros/DrawComparisonsForPaper.C the input files are hard coded in each method to avoid messing up inputs, but if needed, thet should be updated with the local file paths
