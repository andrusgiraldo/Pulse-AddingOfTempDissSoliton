Source code to generate the data and figures for the paper:

Pulse-Adding of Temporal Dissipative Solitons: Resonant Homoclinic Points and the Orbit Flip of Case B with Delay
by Andrus Giraldo and Stefan Ruschel

The source code is organized in the following main folders:
	•	AutoResonantBifurcationSolution: This folder contains data computed from Auto07p that is used to start the continuation runs in DDEBIfTools.
	•	Creating Figures: This folder contains the scripts that generate the different figures in the manuscript. To run these scripts, the necessary files in the folder “DDEResutls” inside the folder “SourceCode2GenerateDataDDEBifTools” should exist.
	•	SourceCode2GenerateDataDDEBifTools: This folder contains the scripts that generate the data for creating the figures. 

At the moment, all the data from DDEBifTools have been saved as “.mat” files inside the “DDEResutls folder”. Hence, running the scripts in the folder “Creating Figures” should generate all the figures in the manuscript after downloading the source code from GitHub. If necessary, one can run the scripts in “SourceCode2GenerateDataDDEBifTools”  to generate the data again. Please change the scripts such that they load the correct address were your installation of DDEBifTools is saved.
