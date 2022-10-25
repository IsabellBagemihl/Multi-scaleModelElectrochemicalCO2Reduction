# Multi-scaleModelElectrochemicalCO2Reduction
This Matlab script is a multi-scale modelling framework for the techno-economic evaluation of a CO2 electrolyser as described in the manuscript
"Techno-economic assessment of CO2 electrolysis: How interdependencies propagate across scales"
DOI: XXXXX
Authors of matlab files: Lucas Cammann and Isabell Bagemihl
Affiliation: Delft University of Technology, Department of Chemical Engineering, van der Maasweg 9, 2629HZ Delft, The Netherlands
Matlab Version: R2020b

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
- The file Optimum.m uses the Matlab built-in non-linear optimiser fmincon. The objective function includes the channel model and the finance calculations. These are summarized in the file obj_func_full.m (M3) and obj_func_simple.m (M2),
  which access the files channelmodel_full.m/channelmodel_sim.m and Finances.m. Furhter a plotting function is embedded in the Optimum.m file.
