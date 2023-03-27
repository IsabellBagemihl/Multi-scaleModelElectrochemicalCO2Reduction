clc;
clear all;

FilePath = 'U:\E2CBWP21\00_Papers\02_MultiscaleModel\01_Written\ACS_SustainableChemistry&Engineering\Model - Copy\DataPoints';
files = get_files(FilePath);

Data;

%% Data from Tan
% Gas flow rate in mL
x_Tan = [5; 10; 20; 40];

Ly = 1.4e-2; %Channel length  [m]
L  = 1e-2;   %Channel height, [m] (-> Denoted as H in the manuscript)
Lw = 1.4e-2; %Channel width,  [m] (-> Denoted as W in the manuscript)
% Liquid flow rate of 12 ml/min
vL = 12*1e-6/60*1/L/Lw; %[m/s]
% Affect of flow rate
CD = 1000;
for i = 1:length(x_Tan)
    v = x_Tan(i)*1e-6/60*1/L/Lw; %[m/s]
    % Ethylene
    [X,FE,y,delP] = channelmodel_full(CD,Ly,v,vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
    
    FE_array_Tan_100(i) = FE;
    Matlab.conversion_Tan_100(i) = X.het;
    Matlab.consumption_Tan_100(i) = X.hom;
end

CD = 2000;
for i = 1:length(x_Tan)
    v = x_Tan(i)*1e-6/60*1/L/Lw; %[m/s]
    % Ethylene
    [X,FE,y,delP] = channelmodel_full(CD,Ly,v,vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
    
    FE_array_Tan_200(i) = FE;
    Matlab.conversion_Tan_200(i) = X.het;
    Matlab.consumption_Tan_200(i) = X.hom;
end

%% Current density plots Tan
figure(1);
errorbar(files.Tan_100(1:4,1),files.Tan_100(1:4,2),files.Tan_100(5:8,2)-files.Tan_100(1:4,2),files.Tan_100(9:12,2)-files.Tan_100(1:4,2),'-or');
hold on;
plot(x_Tan, 100*FE_array_Tan_100,'--r');
hold on;
errorbar(files.Tan_200(1:4,1),files.Tan_200(1:4,2),files.Tan_200(5:8,2)-files.Tan_200(1:4,2),files.Tan_200(9:12,2)-files.Tan_200(1:4,2),'-ob');
hold on;
plot(x_Tan, 200*FE_array_Tan_200,'--b');
xlim([0 45]);
ylim([50 200]);
xlabel('Flow rate in sccm');
ylabel('Partial current density in mA/cm2');
hold off;
%% Data Choi
%Gas flow rate in mL
x_Choi = [20; 30; 40; 50; 100];
Ly = 1.4e-2; %channel length [m]
L = 1e-2;
Lw = 1.4e-2;
%Liquid flow rate of 12 ml/min
vL = 1.2*1e-6/60*1/L/Lw; %[m/s]
CD = 1000;
for i = 1:length(x_Choi)
    v = x_Choi(i)*1e-6/60*1/L/Lw; %[m/s]
    % Ethylene
    [X,FE,y,delP] = channelmodel_full_CO(CD,Ly,v,vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
    
    FE_array_Choi(i) = FE;
    Matlab.conversion_Choi(i) = X.het;
    Matlab.consumption_Choi(i) = X.hom;
end

%% Conversion plots
figure(2);
%Data Tan
plot(files.Tan_Conversion(:,1),files.Tan_Conversion(:,2),'-or');
hold on;
plot(x_Tan, Matlab.conversion_Tan_200*100,'--r');
hold on;
%Data Choi
errorbar(files.Choi_Conversion(1:5,1),files.Choi_Conversion(1:5,2),files.Choi_Conversion(1:5,2)-files.Choi_Conversion(6:10,2),files.Choi_Conversion(11:15,2)-files.Choi_Conversion(1:5,2),'-ob');
hold on;
plot(x_Choi, Matlab.conversion_Choi*100,'--b');
xlim([0 110]);
ylim([0 25]);
xlabel('Flow rate in sccm');
ylabel('Heterogenous conversion in %');
hold off;

%% Data Kas
%Cathode potential vs RHE
E_array = linspace(-0.9,-1.35,7);

%Gas flow rate 5 mL/min
v = 5*1e-6/60*1/L/Lw; %[m/s]

Ly = 1e-2; %channel length [m]
L = 1e-3;
Lw = 1e-2;
%Liquid flow rate of 1 ml/min
vL = 1*1e-6/60*1/L/Lw; %[m/s]

for i = 1:length(E_array)
    E_appl =E_array(i);
    [X,FE,~,~,CD] = channelmodel_full_Ag(E_appl,Ly,v,vL,c_int,k,H,BV,const,L,Lw,y0,por,D,L_c,a);
    
    CD_array_5(i) = CD;
    FE_array_5(i) = FE;
    Matlab.conversion_5(i) = X.het;
    Matlab.consumption_5(i) = X.hom;
    Matlab.conversion_Kas_5(i) = X.het_Kas;
    Matlab.consumption_Kas_5(i) = X.hom_Kas;
end

x = 10; %[mL/min]
v = x*1e-6/60*1/H/W; %[m/s]
for i = 1:length(E_array)
    E_appl = E_array(i);
    [X,FE,~,~,CD] = channelmodel_full_Ag(E_appl,Ly,v,vL,c_int,k,H,BV,const,L,Lw,y0,por,D,L_c,a);

    CD_array_10(i) = CD;
    FE_array_10(i) = FE;
    Matlab.conversion_10(i) = X.het;
    Matlab.consumption_10(i) = X.hom;
    Matlab.conversion_Kas_10(i) = X.het_Kas;
    Matlab.consumption_Kas_10(i) = X.hom_Kas;
end

x = 15; %[mL/min]
v = x*1e-6/60*1/H/W; %[m/s]
for i = 1:length(E_array)
    E_appl = E_array(i);
    [X,FE,~,~,CD] = channelmodel_full_Ag(E_appl,Ly,v,vL,c_int,k,H,BV,const,L,Lw,y0,por,D,L_c,a);

    CD_array_15(i) = CD;
    FE_array_15(i) = FE;
    Matlab.conversion_15(i) = X.het;
    Matlab.consumption_15(i) = X.hom;
    Matlab.conversion_Kas_15(i) = X.het_Kas;
    Matlab.consumption_Kas_15(i) = X.hom_Kas;
end

%% Plot Data from Kas

figure(4);
plot(CD_array_5*0.1,Matlab.conversion_5*100,'--k');
hold on; 
plot(files.conversion_5(1:21,1)*(-1),files.conversion_5(1:21,2),'k');
hold on;
plot(CD_array_10*0.1,Matlab.conversion_10(1:length(Matlab.consumption_5))*100,'--r');
hold on; 
plot(files.conversion_10(1:13,1)*(-1),files.conversion_10(1:13,2),'r');
hold on;
plot(CD_array_15*0.1,Matlab.conversion_15(1:length(Matlab.consumption_5))*100,'--b');
hold on; 
plot(files.conversion_15(1:13,1)*(-1),files.conversion_15(1:13,2),'b');
hold off;

figure(5);
plot(CD_array_5(1:length(Matlab.consumption_5))*0.1,Matlab.consumption_5*100,'--k');
hold on; 
plot(files.consumption_5(1:5,1)*(-1),files.consumption_5(1:5,2),'k');
hold on;
plot(CD_array_10(1:length(Matlab.consumption_5))*0.1,Matlab.consumption_10(1:length(Matlab.consumption_5))*100,'--r');
hold on; 
plot(files.consumption_10(1:5,1)*(-1),files.consumption_10(1:5,2),'r');
hold on;
plot(CD_array_15(1:length(Matlab.consumption_5))*0.1,Matlab.consumption_15(1:length(Matlab.consumption_5))*100,'--b');
hold on; 
plot(files.consumption_15(1:5,1)*(-1),files.consumption_15(1:5,2),'b');
hold off;

function files = get_files(FilePath)
    
%Data from Kas et al.
    mL_data = [FilePath '\5mL_conversion.txt'];
    [files.conversion_5,delimiterOut]=importdata(mL_data);
    mL_data = [FilePath '\5mL_consumption.txt'];
    [files.consumption_5,delimiterOut]=importdata(mL_data);

    mL_data = [FilePath '\10mL_conversion.txt'];
    [files.conversion_10,delimiterOut]=importdata(mL_data);
    mL_data = [FilePath '\10mL_consumption.txt'];
    [files.consumption_10,delimiterOut]=importdata(mL_data);

    mL_data = [FilePath '\15mL_conversion.txt'];
    [files.conversion_15,delimiterOut]=importdata(mL_data);
    mL_data = [FilePath '\15mL_consumption.txt'];
    [files.consumption_15,delimiterOut]=importdata(mL_data);

%Data from Tan et al.
    Current_data = [FilePath '\Tan_Flowrate_100mA.txt'];
    [files.Tan_100,delimiterOut]=importdata(Current_data);

    Current_data = [FilePath '\Tan_Flowrate_200mA.txt'];
    [files.Tan_200,delimiterOut]=importdata(Current_data);
    
    Conversion_data = [FilePath '\Tan_Conversion.txt'];
    [files.Tan_Conversion,delimiterOut]=importdata(Conversion_data);
    
%Data from Choi et al. for formic acid
    Conversion_data = [FilePath '\Choi_Conversion.txt'];
    [files.Choi_Conversion,delimiterOut]=importdata(Conversion_data);
    
end