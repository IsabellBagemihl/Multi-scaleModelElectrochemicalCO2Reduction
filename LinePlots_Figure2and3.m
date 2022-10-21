function [NPV_lim]=LinePlots_Figure2and3
%% Get parameters from Data.m
Data;

% Fixed performance variables for model M1
X.hom = 0;
X.het = 0.5;
FE = 0.7;
V_cell = 3.69;

%Fixed gas flow rate of 10 mL/min
v =  10/60*10^-6/(10^-2*10^-3); %velocity in m/s

%% Limiting NPV for model M1
diff = 0;

CD_1 = 10;
n = 0;
[NPV_old]=Finances(X,FE,CD_1,Ly,v,V_cell,const,1);
while diff < 0.999
%Model 1
    CD_1 =CD_1 + 100;
    X.hom = 0;
    X.het = 0.5;
    FE = 0.7;
    V_cell = 3.69;

    [NPV_new]=Finances(X,FE,CD_1,Ly,v,V_cell,const,1);
    diff = NPV_new/NPV_old;
    NPV_old = NPV_new;
    n = n+1;
    Array_CD(n) = CD_1;
    Array_NPV(n) = NPV_new;
end
NPV_lim = NPV_old;
figure(1);
plot(Array_CD/10, 10^7*Array_NPV);
xlabel('Current density [mA/cm2]');
ylabel('NPV [M$]');

%% Plots for sensitivity towards CD
%Base conditions for the plots
Ec = linspace(-const.T*const.R/(alpha_c*const.F)*log(2200/j0)+E0_C2H4,-const.T*const.R/(alpha_c*const.F)*log(500/j0)+E0_C2H4,40);
j = j0*exp(-alpha_c.*(Ec-E0_C2H4)*const.F/(const.T*const.R));

for i = 1:length(j)
    CD = j(i);
    %Model 3 
%     [X,FE] = channelmodel_full(CD,Ly,v,vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
% 
%     NPV = Finances(X,FE,CD,Ly,v,Ec(i),const,3);
% 
%     Model3.Xhom_CD(i) = X.hom;
%     Model3.Xhet_CD(i) = X.het; 
%     Model3.NPV_CD(i)  = NPV;
%     Model3.FE_CD(i)   = FE;

    %Model 2 
    CDhom = 500;                %Fixed loss current denisty of 50 mA/cm^2 here in A/m^2
    [X,FE] = channelmodel_simp(Ly,v, const.F,y0,CD,CDhom,L);

    NPV = Finances(X,FE,CD,Ly,v,Ec(i),const,2);

    Model2.Xhom_CD(i) = X.hom;
    Model2.Xhet_CD(i) = X.het; 
    Model2.NPV_CD(i)  = NPV;
    Model2.FE_CD(i)   = FE;

    %Model 1
    X.hom = 0;
    X.het = 0.5;
    FE = 0.7;
    V_cell = 3.69;

    [NPV]=Finances(X,FE,CD,Ly,v,V_cell,const,1);

    Model1.Xhom_CD(i) = X.hom;
    Model1.Xhet_CD(i) = X.het; 
    Model1.NPV_CD(i)  = NPV;
    Model1.FE_CD(i)   = FE;
end
    
%% Plots for sensitivity towards vg
v_min = 2.6/60*10^-6/(1e-3*1e-2);
v_max =54/60*10^-6/(1e-3*1e-2);

vg = linspace(v_min, v_max,40);
%fixed current density of 100 mA/cm2
CD = 1000;
Ec = -const.T*const.R/(alpha_c*const.F)*log(CD/j0)+E0_C2H4;

for i = 1:length(vg)
    v = vg(i);
    
    %Full channel model (M3)
    
%     [X,FE] = channelmodel_full(CD,Ly,v,vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
%  
%     NPV = Finances(X,FE,CD,Ly,v,Ec,const,3);
%  
%     Model3.Xhom_v(i) = X.hom;
%     Model3.Xhet_v(i) = X.het; 
%     Model3.NPV_v(i)  = NPV;
%     Model3.FE_v(i)   = FE;
    
    %Simplistic channel model (M2)
    CDhom = 500; %Fixed loss current denisty of 50 mA/cm^2 here in A/m^2
    
    [X,FE] = channelmodel_simp(Ly,v, const.F,y0,CD,CDhom,L);

    NPV = Finances(X,FE,CD,Ly,v,Ec,const,2);

    Model2.Xhom_v(i) = X.hom;
    Model2.Xhet_v(i) = X.het; 
    Model2.NPV_v(i)  = NPV;  
    Model2.FE_v(i)   = FE;
    
    %No channel model (M1)
    X.hom = 0;
    X.het = 0.5;
    FE = 0.7;
    V_cell = 3.69;
    
    NPV = Finances(X,FE,CD,Ly,v,V_cell,const,1);

    Model1.Xhom_v(i) = X.hom;
    Model1.Xhet_v(i) = X.het; 
    Model1.NPV_v(i)  = NPV;
    Model1.FE_v(i)   = FE;
end 

plots(Model1, Model2, Model3,j,vg,NPV_lim);
save('Line_Model3.mat','Model3');
end

%% Plotting
function plots(Model1, Model2, Model3,j,vg,NPV_lim)
  
    %Figure 2
    f1 = figure(2);
    f1.Position = [100 100 360 650];
    %Current density
    subplot(2,1,1)
    plot(j*0.1,NPV_lim./Model1.NPV_CD,'color','#8C8C8C','LineWidth',.5);
    hold on;
    plot(j*0.1,NPV_lim./Model2.NPV_CD,'--','color','#4166D2','LineWidth',.5);
    hold on;
    plot(j*0.1,NPV_lim./Model3.NPV_CD,'color','#41D25A','LineWidth',.5);
    hold off;
    xlim([50 250]);
    xlabel('Current density [mA cm^{-2}]')
    ylabel('NPV_lim/NPV')
    pbaspect([1 1 1])
    legend({'Model 1','Model 2','Model 3'},'Location','southeast')
    %Gas velocity
    subplot(2,1,2)
    plot(vg*10^-5*60*10^6,NPV_lim./Model1.NPV_v,'color','#A0A0A0','LineWidth',.5);
    hold on;
    plot(vg*10^-5*60*10^6,NPV_lim./Model2.NPV_v,'--','color','#0A0AAF','LineWidth',.5);
    hold on;
    plot(vg*10^-5*60*10^6,NPV_lim./Model3.NPV_v,'color','#0A0AAF','LineWidth',.5);
    hold off;
    xlabel('Gas flow rate [sccm min^{-1}]')
    ylabel('NPV_lim/NPV')
    pbaspect([1 1 1])
    annotation('textbox', [0.075, 0.97, 0, 0], 'string', 'a)')
    annotation('textbox', [0.075, 0.5, 0, 0], 'string', 'b)')
    
    % Figure 3
    f2 = figure(3);
    f2.Position = [100 100 80 350];
    subplot(1,3,1)
    plot(j*0.1,Model1.FE_CD,'color','#A0A0A0','LineWidth',.5) 
    hold on;
    plot(j*0.1,Model2.FE_CD,'--','color','#0A0AAF','LineWidth',.5);
    hold on;
    plot(j*0.1,Model3.FE_CD,'color','#0A0AAF','LineWidth',.5);
    hold off;
    ylim([0.4 1]);
    xlim([50 250]);
    xlabel('Current density [mA cm^{-2}]')
    ylabel('FE [-]')
    pbaspect([1 1 1])
    legend({'Model 1','Model 2','Model 3'},'Location','northeast');
    
    
    subplot(1,3,2)
    plot(j*0.1,Model1.Xhet_CD,'color','#A0A0A0','LineWidth',.5) 
    hold on;
    plot(j*0.1,Model2.Xhet_CD,'--','color','#0A0AAF','LineWidth',.5);
    hold on;
    plot(j*0.1,Model3.Xhet_CD,'color','#0A0AAF','LineWidth',.5);
    hold off;
    ylim([0 0.6]);
    xlim([50 250]);
    xlabel('Current density [mA cm^{-2}]')
    ylabel('X_{het} [-]')
    pbaspect([1 1 1])
    
    subplot(1,3,3)
    hold on;
    plot(j*0.1,Model2.Xhom_CD,'--','color','#0A0AAF','LineWidth',.5);
    hold on;
    plot(j*0.1,Model3.Xhom_CD,'color','#0A0AAF','LineWidth',.5);
    hold off;
    xlabel('Current density [mA cm^{-2}]')
    ylabel('X_{hom} [-]')
    xlim([50 250]);
    pbaspect([1 1 1])
    
end

