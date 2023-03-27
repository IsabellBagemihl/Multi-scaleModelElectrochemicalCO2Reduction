This Matlab script is a multi-scale modelling framework for the techno-economic evaluation of a CO2 electrolyser as described in the manuscript
"Techno-economic assessment of CO2 electrolysis: How interdependencies propagate across scales"
DOI: XXXXX
Authors: Isabell Bagemihl and Lucas Cammann
Affiliation: Delft University of Technology, Department of Chemical Engineering, van der Maasweg 9, 2629HZ Delft, The Netherlands

The structure of the files is as following:
- The file Data.m contains all relevant parameters as given in the supplementary info of the manuscript above
- The simplistic channel model (M2) is given in the file channelmodel_sim.m
- The full channel model (M3) is given in the file channelmodel_full.m
- The electrolyser scale and process scale equations are included in the file Finances.m
- The file LinePlots_Figure2and3.m contains a function which accesses the Data.m, channelmodel_sim.m, channelmodel_full.m and Finances.m file to calculate the
  channel scale output for the different models and calculates the NPV for varying gas flowrates and current densities. The no channel model (M1) is furhter defined in this file.
  This file further contains a function for plotting the data.
- The file contourplots.m contains a function which accesses the Data.m, channelmodel_sim.m, channelmodel_full.m and Finances.m file to calculate the
  channel scale output for the different models and calculates the NPV for varying gas flowrates and current densities. The no channel model (M1) is furhter defined in this file.
  This file further contains a function for plotting the data.
- The file Optimum.m 
------------------------------------------------------------------------------------
- Datafitting.m is used to fit the Tafel equation with the experimental data by Tan et al. to obtain the kinetic constants (Section 2.4 supplementary)
- Validation_main.m calls the different channel model files and reads and plots the experimental data from Tan, Choi and Kas (Section 2.5 supplementary)
- channelmodel_full_CO.m conatins the full channel model (M3) for a two electron transfer reaction
- channelmodel_full_Ag.m conatins the full channel model (M3) with the Tafel equations from Kas et al.
- Sensetivity_Analysis.m calls the file Data_Sensetivity.m and TorPlot.m (Section 4.1 supplementary) 
