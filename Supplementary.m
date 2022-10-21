%function Supplementary(v)

v =  10/60*10^-6/(10^-2*10^-3); %8 mL/min

Data;
% Values Jouney
X.hom = 0;
X.het = 0.5;
FE = 0.9;

%% Figure a.) variable current density
CD_array = linspace(50, 2500, 100);
%CD_array = 2000;
V_cell_Jouny = 2.0;
V_cell_fixed = 3.69;

Ec = -const.T*const.R/(alpha_c*const.F)*log(CD_array./j0)+E0_C2H4;
eta.actA = const.R*const.T/(0.5*const.F)*asinh(CD_array./(2*1e-7));
eta.ohm = CD_array.*(L/sigma_el+Lm/sigma_m);
V_cell_array = abs(Ec) + eta.actA + eta.ohm + 1.23;

for i=1:length(CD_array)
    
    CD = CD_array(i);
    V_cell = V_cell_array(i);
    
    [NPV.Jouny(i)]=Finances1(X,FE,CD,Ly,v,V_cell_Jouny,const);
    [NPV.fixed(i)]=Finances1(X,FE,CD,Ly,v,V_cell_fixed,const);
    [NPV.variable(i)]=Finances1(X,FE,CD,Ly,v,V_cell,const);
    
end

conversion_array = linspace(0.1,0.9,100);
CD = 2000;
for i=1:length(conversion_array)
    
    X.het = conversion_array(i);
    
    [NPV.Jouny1(i)]=Finances1(X,FE,CD,Ly,v,V_cell_Jouny,const);
    [NPV.fixed1(i)]=Finances1(X,FE,CD,Ly,v,V_cell_fixed,const);
    
end

for i=1:length(CD_array)
     CD = CD_array(i);
    FE = 0.7;
    [NPV.variable1(i)]=Finances1(X,FE,CD,Ly,v,V_cell,const);
    X.het = 0.3;
    [NPV.variable2(i)]=Finances1(X,FE,CD,Ly,v,V_cell,const);
    
end

v_min = 2.6/60*10^-6/(1e-3*1e-2);
v_max =54/60*10^-6/(1e-3*1e-2);

vg = linspace(v_min, v_max,40);
%fixed current density of 100 mA/cm2
CD = 1000;
Ec = -const.T*const.R/(alpha_c*const.F)*log(CD/j0)+E0_C2H4;

for i = 1:length(vg)
    v = vg(i);
    
    %Model 2
    CDhom = 500;
    [X,FE] = channelmodel_simp(Ly,v, const.F,y0,CD,CDhom, Lw,L);
    Hom(i) = X.hom;
end

figure(1);
plot(CD_array*0.1,NPV.Jouny*-1,'b','LineWidth',.5);
hold on;
plot(CD_array*0.1,NPV.variable*-1,'r','LineWidth',.5);
hold on;
plot(CD_array*0.1,NPV.fixed*-1,'g','LineWidth',.5);
hold off;
xlim([50 250]);
ylim([-3e+7 1e+7]);
xlabel('Current density i_{tot} [mA cm-2]')
ylabel('NPV [$]')
pbaspect([1 1 1])
legend('V_{cell}=2V','V_{cell}=f(i_{tot})','V_{cell}=3.53V');

figure(2);
plot(conversion_array,NPV.Jouny1*-1,'b','LineWidth',.5);
hold on;
plot(conversion_array,NPV.fixed1*-1,'g','LineWidth',.5);
hold off;
xlim([0 1]);
ylim([-2.5e+7 1e+7]);
xlabel('Heterogenous conversion')
ylabel('NPV [$]')
pbaspect([1 1 1])
legend('V_{cell}=2V','V_{cell}=3.53V');

figure(3);
plot(CD_array*0.1,NPV.variable1*-1,'b','LineWidth',.5);
hold on;
plot(CD_array*0.1,NPV.variable2*-1,'r','LineWidth',.5);
hold off;
%xlim([50 250]);
%ylim([-2.5e+7 1e+7]);
xlabel('Current density i_{tot} [mA cm-2]')
ylabel('NPV [$]')
pbaspect([1 1 1])
legend('0.5','0.3');
%end

figure(4);
plot(vg*60/10^-6*(10^-2*10^-3),Hom);
xlabel('Gas flow rate in sscm min^{-1}')
ylabel('Homogenous conversion')
%Model 1
function [NPV]=Finances1(X,FE,CD,Ly,v,V_cell,const)
    Data;
    SetupBest;
    
    %Base values
    ProdRate    = 10000;        % kg d^-1
    lifetime    = 20;           % a
    optime      = 350;          % d per a
    IR          = 0.1;          % Interest rate
    
    %Flow rates
    Fl.Out = ProdRate/0.028;         %Moles C2H4 d^-1
    A = Fl.Out*12*const.F/(CD*FE*3600*24); %Electrolyzer area 
    Fl.tot = Fl.Out*2/X.het;         %Total flowrate through electrolyzer
    Fl.elec = Fl.tot*X.hom;          %CO2 lost in hom. RX
    Fl.In = 2*Fl.Out + Fl.elec;      %Total CO2 feed to system
    Fl.rec = Fl.tot - Fl.In;         %Recycle flowrate
    n = A/(Lw*Ly);                   %Channel number
    Fl.tot2 = n*v*Lw*L*y0(1);        %As sanity check for Fl.tot
    CO2Rate = 44/28*2*ProdRate ...
            + ProdRate*X.hom/X.het;  %Annual CO2 consumption
    V = Fl.tot/(y0(1)*3600*24);      %Conversion from mole/d to m^3/s
    
    %Overpotential
    eta.tot = V_cell;
    
    Vol = eta.tot;
    Power = A*CD*Vol; 

    %Economics
    %Investment cost
    Price.Inv  = (A*Price.A)*(1+35/65) + (Price.PSA*(V*3600/1000)^scale);        
    %Revenue 
    Price.Rev  =  (ProdRate*Price.Prod - CO2Rate*Price.CO2)*optime; 
    %Operating cost 
    Price.Op   = ((Power)*8.4 + 0.25*V*3600*24*optime)*Price.Elec*1.025; 

    CF = zeros(1,lifetime-1);
    for i = 1:lifetime-1
        CF(i) = (Price.Rev - Price.Op)/(1 + IR)^i; 
    end 
    CCF = sum(CF(:));
    NPV = -(- Price.Inv + CCF); 
end