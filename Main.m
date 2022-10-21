%%%%% Plotting of Figures %%%
clc;
clear variables;
close all;
%Optimum ;

clc;
clear variables;
[NPV_lim]=LinePlots_Figure2and3;

clc;
clearvars -except NPV_lim;
contourplots(NPV_lim);