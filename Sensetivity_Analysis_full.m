clc;
clear all;

% Base Case
Price.A     = 920;          % $ per m^2 
Price.PSA   = 1990000;      % $  
Price.Elec  = 0.03;         % $ per kWh
Price.Prod  = 1.3;          % $ per kg 
Price.CO2   = 0.04;         % $ per kg
ProdRate    = 10000;        % kg d^-1
lifetime    = 20;           % a
IR          = 0.1;          % Interest rate
CD          = 2000;         %CD base case 100 mA/cm2
v           = 10/60*10^-6/(10^-2*10^-3); %Base case 10 sccm

data = [lifetime, ProdRate, IR, Price.PSA, Price.A, Price.CO2, Price.Elec, Price.Prod, CD, v];
names = {'Electrolyser life time';'Production rate';'Interest rate';'PSA cost';'Electrolyser cost';'CO2 cost';'Electricity price';'Selling price';'Total current density'; 'Gas velocity'};

% Better Case
Price.A     = 450;          % $ per m^2 
Price.PSA   = 1990000*0.8;      % $  
Price.Elec  = 0.02;         % $ per kWh
Price.Prod  = 1.3*1.15;          % $ per kg 
Price.CO2   = 0;         % $ per kg
ProdRate    = 10000*0.8;        % kg d^-1
lifetime    = 20*0.8;           % a
IR          = 0.1*1.2;          % Interest rate
CD          = 3000;             %CD base case 100 mA/cm2
v           = 50/60*10^-6/(10^-2*10^-3); %Base case 10 sccm

data_better = [lifetime, ProdRate, IR, Price.PSA, Price.A, Price.CO2, Price.Elec, Price.Prod, CD, v];

% Worse Case
Price.A     = 1840;          % $ per m^2 
Price.PSA   = 1990000*1.2;      % $  
Price.Elec  = 0.04;         % $ per kWh
Price.Prod  = 1.3*0.85;          % $ per kg 
Price.CO2   = 0.07;         % $ per kg
ProdRate    = 10000*1.2;        % kg d^-1
lifetime    = 20*1.2;           % a
IR          = 0.1*0.8;          % Interest rate
CD          = 1000;         %CD base case 100 mA/cm2
v           = 5/60*10^-6/(10^-2*10^-3); %Base case 10 sccm

data_worse = [lifetime, ProdRate, IR, Price.PSA, Price.A, Price.CO2, Price.Elec, Price.Prod, CD, v];

TornadoPlot(data, data_worse, data_better,names, @Finances_sens);

function [NPV] = Finances_sens(y,X,FE)
%This function calculates the NPV of the process as a function of various 
%parameters. In doing so it calls the relevant financial and physical data
%as well as the single channel model. For the optimization this function is
%itself called by optimization.m

%% Initialize the run
%Note: Upon actuating the below switches the user might want to change one
%of the corresponding variables, which is described in further detail in
%the definitions below.
lifetime    = y(1);         % a
ProdRate    = y(2);         % kg d^-1
IR          = y(3);         % Interest rate
Price.PSA   = y(4);         % $ 
Price.A     = y(5);         % $ per m^2 
Price.CO2   = y(6);         % $ per kg
Price.Elec  = y(7);         % $ per kWh
Price.Prod  = y(8);         % $ per kg 
CD          = y(9);
v           = y(10);

Data_Sensetivity;

optime      = 350;          % d per a
scale       = 0.7;
Price.Stack = 250;          % $ per kW


%% Channel solution 
%Solution from optimisation
Ec = -const.T*const.R/(alpha_c*const.F)*log(2200/j0)+E0_C2H4;

%% Calculate flowrates and reactor characteristics
Fl.Out = ProdRate/0.028;         %Moles C2H4 d^-1
A = Fl.Out*12*const.F/(CD*FE*3600*24); %Electrolyzer area 
Fl.tot = Fl.Out*2/X.het;         %Total flowrate through electrolyzer
Fl.elec = Fl.tot*X.hom;          %CO2 lost in hom. RX
Fl.In = 2*Fl.Out + Fl.elec;      %Total CO2 feed to system
Fl.rec = Fl.tot - Fl.In;         %Recycle flowrate
n = A/(Lw*Ly);                   %Channel number
%Fl.tot2 = n*v*Lw*L*y0(1);        %As sanity check for Fl.tot
CO2Rate = (1 + X.hom/X.het)*44/28*2*ProdRate; %Annual CO2 consumption in kg per year
V = Fl.tot/(y0(1)*3600*24);      %Conversion from mole/d to m^3/s


%% Calculate economic for different anodic reactants 
%Calculate cell voltage 
eta.actA = const.R*const.T/(0.5*const.F)*asinh(CD/(2*1e-7));
eta.ohm = CD*(L/sigma_el+Lm/sigma_m);
%eta.tot = abs(Ec) + eta.actA + eta.ohm + 1.23;
eta.actA    = const.R*const.T/(0.5*const.F)*asinh(CD/(2*1e-7));
eta.ohm     = CD*(L/sigma_el+Lm/sigma_m);
eta.tot     = 1.23 + eta.actA + E0_C2H4 + abs(Ec) + eta.ohm;
Power = A*CD*eta.tot;   

%Economics
%Investment cost
Price.Inv  = (A*Price.A)*(1+35/65) + (Price.PSA*(V*3600/1000)^scale);       
%Revenue 
Price.Rev  =  (ProdRate*Price.Prod - CO2Rate*Price.CO2)*optime; 
%Operating cost 
Price.Op   = ((Power)*8.4 + 0.25*V*3600*24*optime)*Price.Elec*1.025; 

%% Final NPV calculation 
CF = zeros(1,lifetime-1);
for i = 1:lifetime-1
    CF(i) = (Price.Rev - Price.Op)/(1 + IR)^i; 
end 
CCF = sum(CF(:));
NPV = (- Price.Inv + CCF)*10^-6; 

end 

function TornadoPlot(data, data_worse, data_better,names,fh)
Data_Sensetivity;
[X,FE,y,delP] = channelmodel_full(data(9),Ly,data(10),vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);

for i=1:length(names)
    Objective_low = data;
    Objective_high = data;
    Objective_low(i)=data_worse(i);
    Objective_high(i)=data_better(i);
    if i >= 9
        [X_worse,FE_worse,y,delP] = channelmodel_full(Objective_low(9),Ly,Objective_low(10),vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
        [X_better,FE_better,y,delP] = channelmodel_full(Objective_high(9),Ly,Objective_high(10),vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
        Objective_low_sum(i)=fh(Objective_low, X_worse, FE_worse);
        Objective_high_sum(i)=fh(Objective_high, X_better, FE_better);
    else
        Objective_low_sum(i)=fh(Objective_low, X, FE);
        Objective_high_sum(i)=fh(Objective_high, X, FE);
    end
    low(i)=Objective_low_sum(i);
    high(i)=Objective_high_sum(i); 

    % The base value is where the y axis is centered
    Objective_base_value=fh(data,X,FE);
end
names_Objective=names;
% Create a figure and plot the low and high horizontally
figure
h = barh(Objective_high_sum);
hold on
xmin=min([min(Objective_low_sum),min(Objective_high_sum)]);
xmax=max([max(Objective_low_sum),max(Objective_high_sum)]);
xlim([1.025*xmin 0.975*xmax])
barh(Objective_low_sum,'r')
bh = get(h,'BaseLine');
set(bh,'BaseValue',Objective_base_value);
title('Sensitivities')

set(gca,'yticklabel',names)
set(gca,'Ytick',[1:length(names)],'YTickLabel',[1:length(names)])
set(gca,'yticklabel',names_Objective)

xlabel('NPV ($ millions)')

end