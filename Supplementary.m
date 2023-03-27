%function Supplementary(v)

v =  10/60*10^-6/(10^-2*10^-3); %8 mL/min

% Values this work
Data;
V_cell_fixed = 3.69;
% Values Jouney
X.hom = 0;
X.het = 0.5;
FE = 0.9;
V_cell_Jouny = 2.0;

%% Figure a.) variable current density
CD_array = linspace(50, 2500, 100);
ZeroLine = linspace(0, 0, 100);

Ec = -const.T*const.R/(alpha_c*const.F)*log(CD_array./j0)+E0_C2H4;
eta.actA = const.R*const.T/(0.5*const.F)*asinh(CD_array./(2*1e-7));
eta.ohm = CD_array.*(L/sigma_el+Lm/sigma_m);
V_cell_array = abs(Ec) + eta.actA + eta.ohm + 1.23;

for i=1:length(CD_array)
    
    CD = CD_array(i);
    V_cell = V_cell_array(i);
    
    [NPV.Jouny(i),~]=Finances(X,FE,CD,Ly,v,V_cell_Jouny,const,1);
    [NPV.fixed(i)]=Finances(X,FE,CD,Ly,v,V_cell_fixed,const,1);
    [NPV.variable(i)]=Finances(X,FE,CD,Ly,v,V_cell,const,1);
    
end

[x,y]       = polyxpoly(CD_array,NPV.fixed,CD_array,NPV.variable);
Ec          = -const.T*const.R/(alpha_c*const.F)*log(x/j0)+E0_C2H4;
eta.actA    = const.R*const.T/(0.5*const.F)*asinh(x/(2*1e-7));
eta.ohm     = x*(L/sigma_el+Lm/sigma_m);
V_cell_Intersect = abs(Ec) + eta.actA + eta.ohm + 1.23;

conversion_array = linspace(0.1,0.9,100);
CD = 2000;
for i=1:length(conversion_array)
    
    X.het = conversion_array(i);
    
    [NPV.Jouny1(i)]=Finances(X,FE,CD,Ly,v,V_cell_Jouny,const,1);
    [NPV.fixed1(i)]=Finances(X,FE,CD,Ly,v,V_cell_fixed,const,1);
    
end

for i=1:length(CD_array)
     CD = CD_array(i);
    FE = 0.7;
    [NPV.variable1(i)]=Finances(X,FE,CD,Ly,v,V_cell,const,1);
    X.het = 0.3;
    [NPV.variable2(i)]=Finances(X,FE,CD,Ly,v,V_cell,const,1);
    
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
    [X,FE] = channelmodel_simp(Ly,v, const.F,y0,CD,CDhom,L);
    Hom(i) = X.hom;
end

figure(1);
plot(CD_array*0.1,NPV.Jouny,'b','LineWidth',.5);
hold on;
plot(CD_array*0.1,NPV.variable,'r','LineWidth',.5);
hold on;
plot(CD_array*0.1,NPV.fixed,'g','LineWidth',.5);
plot(CD_array*0.1, ZeroLine,'--k','LineWidth',.5);
hold off;
xlim([50 250]);
%ylim([-3e+7 1e+7]);
xlabel('Current density i_{tot} [mA cm-2]')
ylabel('NPV [$]')
pbaspect([1 1 1])
legend('V_{cell}=2V','V_{cell}=f(i_{tot})','V_{cell}=3.53V');

figure(2);
plot(conversion_array,NPV.Jouny1,'b','LineWidth',.5);
hold on;
plot(conversion_array,NPV.fixed1,'g','LineWidth',.5);
hold off;
xlim([0 1]);
%ylim([-2.5e+7 1e+7]);
xlabel('Heterogenous conversion')
ylabel('NPV [$]')
pbaspect([1 1 1])
legend('V_{cell}=2V','V_{cell}=3.53V');

figure(3);
plot(CD_array*0.1,NPV.variable1,'b','LineWidth',.5);
hold on;
plot(CD_array*0.1,NPV.variable2,'r','LineWidth',.5);
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
