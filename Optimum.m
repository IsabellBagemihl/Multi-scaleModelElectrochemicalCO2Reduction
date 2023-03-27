function [Model2, Model3] = Optimum
%% Get parameters from Data.m
    SetupBest;
    Data;
%% Optimization Solver
    %options = optimset('PlotFcns',@optimplotfval,'MaxIter', 30);
    options = optimset('PlotFcns',@optimplotfval);
%     guess_f  = [0.03 -0.8];
%     [Model3(1,:)] = fminsearch(@(x) obj_func_full(x),guess_f,options);
%     Model3(1,3)     = j0*exp(-alpha_c*(Model3(1,2)-E0_C2H4)*const.F/(const.T*const.R));
%     [X,FE,y,delP]   = channelmodel_full(Model3(1,3),Ly,Model3(1,1),vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
%     Model3(1,4)     = X.het;
%     Model3(1,5)     = X.hom;
%     Model3(1,6)     = FE;
%     [Model3(1,7),Model3(1,8),Model3(1,9),Model3(1,10),Model3(1,11),Model3(1,12)]=Finances(X,FE,Model3(1,3),Ly,Model3(1,1),Model3(1,2),const,3);
% 
%     save('Optimum_Model3.mat','Model3');

    % Simplistic channel model M2
    guess_s = [0.01 -0.9];             %Initial guess, [v and Ec]
    CDhom   = 500;                      %Fixed loss current denisty of 50 mA/cm^2 here in A/m^2
    
    [Model2(1,:)]   = fmincon(@(x) obj_func_simple(x), guess_s(1,:), [],[],[],[],[0.01, -0.9], [10, -0.55],[],options);
    Model2(1,3)     = j0*exp(-alpha_c*(Model2(1,2)-E0_C2H4)*const.F/(const.T*const.R));
    [X,FE,delP]     = channelmodel_simp(Ly,Model2(1,1), const.F,y0,Model2(1,3),CDhom,L);
    Model2(1,4)     = X.het;
    Model2(1,5)     = X.hom;
    Model2(1,6)     = FE;
    [Model2(1,7),Model2(1,8),Model2(1,9),Model2(1,10),Model2(1,11),Model2(1,12)]=Finances(X,FE,Model2(1,3),Ly,Model2(1,1),Model2(1,2),const,2);
    
%% Literature values
    Jouny = [250 0.5];
    Verma = [200 0.5];
    
%% Plot
    figure(10);
    scatter(0.5,250,'o', 'MarkerFaceColor', 'g');
    hold on;
    scatter(Model2(1,4),Model2(1,3)/10,'o', 'MarkerFaceColor', 'b');
    hold on;
    scatter(Model3(1,4),Model3(1,3)/10,'o', 'MarkerFaceColor', 'r');
    hold on;
    scatter(0.5,250,'s', 'MarkerFaceColor', 'k');
    hold on;
    scatter(0.5,200,'s', 'MarkerFaceColor', 'k');
    hold off;
    ylim([0 250]);
    xlim([0 0.5]);
    ylabel('Current density [mA cm^{-2}]')
    xlabel('Conversion')
    pbaspect([1 1 1])
end

