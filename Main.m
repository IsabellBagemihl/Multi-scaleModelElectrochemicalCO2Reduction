%%%%% Plotting of Figures %%%
clc;
clear variables;
close all;
[Opt_Model2, Opt_Model3] = Optimum ;

clc;
clearvars -except Opt_Model2 Opt_Model3;
[NPV_lim]=LinePlots_Figure2and3;

clc;
clearvars -except Opt_Model2 Opt_Model3 NPV_lim;
contourplots(NPV_lim, Opt_Model2, Opt_Model3);