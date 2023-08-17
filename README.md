Source code to generate the data and figures for the paper:

Pulse-Adding of Temporal Dissipative Solitons: Resonant Homoclinic Points and the Orbit Flip of Case B with Delay
by Andrus Giraldo and Stefan Ruschel

The source code is organized in the following main folders:

	- AutoResonantBifurcationSolution: This folder contains data computed from Auto07p that is used to start the continuation runs in DDE-Biftool.
	- Creating Figures: This folder contains the scripts that generate the different figures in the manuscript. To run these scripts, the necessary files in the folder “DDEResults” inside the folder “SourceCode2GenerateDataDDEBifTools” should exist.
	- SourceCode2GenerateDataDDEBifTools: This folder contains the scripts that generate the data for creating the figures. 

At the moment, all the data from DDE-Biftool have been saved as “.mat” files inside the “DDEResults folder” except for the following files (due to their size):

	- SlicesResOri.mat: Run the script "mainSandResSlicesPaper.m" to generate this file.
	- SpecialPDBranchFig3: Run the script "mainExtraBifCurvesPaperBiftoolGit.m" to generate this file.

After creating these last two files, running the scripts in the folder “Creating Figures” should generate all the figures in the manuscript after downloading the source code from GitHub. If necessary, one can run the scripts in “SourceCode2GenerateDataDDEBifTools”  to regenerate the data again. Please change the scripts such that they load the correct address where your installation of DDE-Biftool is saved (you need the latest commit of DDE-Biftool on Github to run mainExtraBifCurvesPaperBiftoolGit.m).
