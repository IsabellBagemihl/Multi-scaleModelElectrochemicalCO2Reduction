function Optimum
%% Get parameters from Data.m
    SetupBest;
    Data;
%% Optimization Solver    
    options = optimoptions('fmincon');
    
    % Full channel model M3
    guess_f  = [0.02 -0.79];            %Initial guess, [v and Ec]
    
    [Model3(1,:)]   = fmincon(@(x) obj_func_full(x), guess_f(1,:), [],[],[],[],[0.01, -0.85], [10, -0.7],[],options);
    Model3(1,3)     = j0*exp(-alpha_c*(Model3(1,2)-E0_C2H4)*const.F/(const.T*const.R));
    [X,FE,y,delP]   = channelmodel_full(Model3(1,3),Ly,Model3(1,1),vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
    Model3(1,4)     = X.het;
    
    % Simplistic channel model M2
    guess_s = [0.01 -0.85];             %Initial guess, [v and Ec]
    CDhom   = 500;                      %Fixed loss current denisty of 50 mA/cm^2 here in A/m^2
    
    [Model2(1,:)]   = fmincon(@(x) obj_func_simple(x), guess_s(1,:), [],[],[],[],[0.01, -0.85], [10, -0.55],[],options);
    Model2(1,3)     = j0*exp(-alpha_c*(Model2(1,2)-E0_C2H4)*const.F/(const.T*const.R));
    [X,FE,delP]     = channelmodel_simp(Ly,Model2(1,1), const.F,y0,Model2(1,3),CDhom,L);
    Model2(1,4)     = X.het;
    
    save('Optimum_Model3.mat','Model3');
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

