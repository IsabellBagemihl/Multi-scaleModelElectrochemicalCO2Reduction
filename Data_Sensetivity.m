%% Data 
% This file contains the relevant data which is shared for all simulations,
% hence physical constants, fixed dimensions, operating conditions etc.
%% Fixed constants 
 const.F    = 96485;                         %Faraday constant, [C mol^-1]
 const.R    = 8.314;                         %Universal gas constant
 const.T    = 298.15;                           %Temperature, [K]
 const.p    = 1e5;                           %Pressure in channel [Pa]
 
 L    = 1e-3;                          %channel thickness, [m] (-> Denoted as H in the thesis)
 Lw   = 1e-2;                          %channel width, [m] (-> Denoted as W in the thesis)
 L_c  = 3e-6;                          %Catalyst layer thickness, [m] (-> Denoted as H_c in the thesis)
 Ly = 0.1;
 Re = 1100;
 Lm = 115e-6;                           %Thickness of the membrane
 vis_H2O = 0.893e-6;                    %Dynamic viscosity (p=1bar, T = 298,15K) [m^2 s^-1]
 dh    = 2*Lw*L/(L+Lw);                 %Hydraulic diameter https://www.engineeringtoolbox.com/hydraulic-equivalent-diameter-d_458.html
 vL = Re/dh*(vis_H2O);                  %Liquid velocity [m/s]
 
 rnp  = 60e-9;                         %Radius nanoparticles, [m]
 por  = 0.7;                           %Porosity, [-]
 a    = 1e7;                           %Specific area for reaction, [m^-1]
 y0   = [const.p/(8.314*const.T) 0 0];             %Inlet concentrations
 
 j0 = 0.22;
 alpha_c = 0.25;
 
 sigma_el = 5.5;                       %Electrolyte conductivity [S m^-1]
 sigma_m = 9.3;
 
 Electrolyte = 1;                       % 1M KHCO3
 [CO2, HCO3, CO3, OH] = Concentration_Electrolyte(const, Electrolyte);
 pH = 14+log10(OH*10^(-3));
 
 E0_C2H4 = 0.08;
%% Data structs
% Reaction constants for carbon equilibrium reaction (Gupta)
k = struct('f1'  , 5.93,...      Basic Diss CO2, [m^3 (mol s)^-1]
           'r1'  , 1.34*1e-4,... Basic Diss CO2, [s^-1]
           'f2'  , 1e5,...       Diss HCO3- [m^3 (mol s)^-1]
           'r2'  , 2.15*1e4...   Diss HCO3- [s^-1]
           );     
% Diffusion coefficients in water      
D = struct('CO2'   , 1.91e-9, ... All in [m^2 s^-1]
           'CO'    , 2.03e-9, ...
           'H2'    , 4.5e-9,  ...
           'H'     , 9.3e-9,  ...
           'OH'    , 5.3e-9,  ...
           'HCOOH' , 1.5e-9,  ...
           'CO3'   , 0.92e-9, ...
           'HCO3'  , 1.19e-9, ...
           'CH4'   , 1.88e-9, ...
           'C2H4'  , 1.09e-9, ...
           'EtOH'  , 1.23e-9, ...
           'iBuOH' , 0.95e-9  ...
            );
        
c_int = struct('CO2' , CO2,           ...
               'OH'  , OH,         ... 
               'HCO3', HCO3,            ...
               'CO3' , CO3,            ...  
               'C2H4', 0);                        

 H = struct('CO2' ,   0.85,   ... All [-] (Henry constant in concentration form) Conversion from [mol/(kg bar)] by multplying with 25
            'C2H4',   0.01,   ...
            'H2'  ,   0.0195, ...
            'CO'  ,   0.02475);
 
 
 %% Calculations
 dh    = 2*Lw*L/(L+Lw);                %Hydraulic diameter https://www.engineeringtoolbox.com/hydraulic-equivalent-diameter-d_458.html

 %% Function for electrolyte concentration
 function [CO2, HCO3, CO3, OH] = Concentration_Electrolyte(const, Electrolyte)
 
hg = -0.0172-0.000338*(298-298.15);
Ks = (0.0922+hg)+(0.0967+hg);
CO2water = const.p*0.03406e3*10^-5;
CO2 = CO2water * 10^(-Ks*Electrolyte);

%homogenous reaction rate constants
Ionic_Strength = 0.5*(Electrolyte*1*1^2+Electrolyte*1*(-1)^2);
S = 1000*Ionic_Strength/(19.92+1.0049*Ionic_Strength);
pK01 = -126.34048 + 6320.813/const.T + 19.568224*log(const.T);   % http://www.sciencedirect.com/science/article/pii/S0304420305001921   Dissociation constants of carbonic acid in seawater as a function of salinity and temperature
A1 = 13.4191*S^0.5+0.0331*S-5.33E-5*S^2;
B1 = -530.123*S^0.5-6.103.*S;
C1 = -2.06950*S^0.5;

pK1 = pK01 + A1+B1/const.T+C1*log(const.T);

A2 = 21.0894*S^0.5+0.1248*S-3.687E-4*S^2;    B2 = -772.483*S^0.5-20.051*S;    C2 = -3.3336*S^0.5; 
pK02= -90.18333+5143.692/const.T+14.613358*log(const.T);
pK2 = pK02 + A2+B2/const.T+C2*log(const.T);

Kw = 1e-14;                 
K1 = 10^(-pK1);             
K2 = 10^(-pK2);

%Electrolyte
HCO3 = Electrolyte*1e3;
H = K1*CO2/HCO3;
OH = Kw/H*10^3;
CO3 = K1*K2*CO2/H^2;

end
