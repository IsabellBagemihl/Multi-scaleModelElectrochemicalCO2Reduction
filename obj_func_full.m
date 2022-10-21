function [NPV] = obj_func_full(x)
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

%% Output full channel model (M3) 
v  = x(1);
Ec = x(2);

vL = Re/dh*(vis_H2O);
CD = j0*exp(-alpha_c*(Ec-E0_C2H4)*const.F/(const.T*const.R));  %CD as f(Ec)

[X,FE,y,delP] = channelmodel_full(CD,Ly,v,vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);

%% Output Process scale model
NPV = Finances(X,FE,CD,Ly,v,Ec,const,3);

end 

