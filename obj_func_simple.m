function [NPV] = obj_func_simple(x)
%This function calculates the NPV of the process as a function of various 
%parameters. In doing so it calls the relevant financial and physical data
%as well as the single channel model. For the optimization this function is
%itself called by optimization.m

%% Initialize the run
%Note: Upon actuating the below switches the user might want to change one
%of the corresponding variables, which is described in further detail in
%the definitions below.
Data;
SetupBest;


%% Output simplistic channel model (M2) 
v  = x(1);
Ec = x(2);

CD = j0*exp(-alpha_c*(Ec-E0_C2H4)*const.F/(const.T*const.R));  %CD as f(Ec)
CDhom = 500;     %Constant homogeneous consumption rate, set in the units of A m^-2 for comparability

[X,FE,delP] = channelmodel_simp(Ly,v, const.F,y0,CD,CDhom,L);

%% Output Process scale model
[NPV,~] = Finances(X,FE,CD,Ly,v,Ec,const,2);
NPV =-1*NPV;
end 

