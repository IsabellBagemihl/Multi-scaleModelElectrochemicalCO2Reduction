function contourplots(NPV_lim)
%% Generation of colormap
% Orange - Grey
map(1,:) = [230,146,0];
map(2,:) = [230, 162, 52];
map(3,:) = [230,178,105];
map(4,:) = [230, 194,158];
map(5,:) = [230,210,210];
map(6,:) = [230, 220, 220];
map(7,:) = [230, 230, 230];

[X1,Y1]     = meshgrid([1:3],[1:90]);
map         = interp2(X1([1,15,30,45,60,75,90],:),Y1([1,15,30,45,60,75,90],:),map,X1,Y1);
map         = uint8(map);
%% Get Data
SetupBest;
Data;

%% Boundaries for gas velcoity and current density
v_min = 2.5/60*10^-6/(L*Lw);
v_max =50/60*10^-6/(L*Lw);
Mesh(1,:) = linspace(v_min,v_max,20);              %Velocity [m/s]
Ec_min = -const.T*const.R/(alpha_c*const.F)*log(500/j0)+E0_C2H4;
Ec_max = -const.T*const.R/(alpha_c*const.F)*log(2500/j0)+E0_C2H4;
Mesh(2,:) = linspace(Ec_min,Ec_max,20);           %Cathode voltage [V]

%Full channel model (M3)
vL = Re/dh*(vis_H2O);
for i = 1:size(Mesh,2)
    for j = 1:size(Mesh,2)
        x(1) = Mesh(1,i);
         x(2) = Mesh(2,j);
         
         CD = j0*exp(-alpha_c*(x(2)-E0_C2H4)*const.F/(const.T*const.R));
         
         [X,FE,y,delP] = channelmodel_full(CD,Ly,x(1),vL,c_int,k,H,const.F,L,Lw,y0,por,D,L_c,a);
         
         [Model3(i,j)]= Finances(X,FE,CD,Ly,x(1),x(2),const,3);
    end 
end 
save('Contour_Model3.mat','Model3');

% Simplistic channel model (M2)
for i = 1:size(Mesh,2)
    for j = 1:size(Mesh,2)
        x(1) = Mesh(1,i);
        x(2) = Mesh(2,j);
        
        CD = j0*exp(-alpha_c*(x(2)-E0_C2H4)*const.F/(const.T*const.R));
        CDhom = 500;                     %Fixed loss current denisty of 50 mA/cm^2 here in A/m^2
        
        [X,FE] = channelmodel_simp(Ly,x(1), const.F,y0,CD,CDhom,L);
        
        [Model2(i,j)]= Finances(X,FE,CD,Ly,x(1),x(2),const,2);
    end 
end 

%No channel model (M1)
for i = 1:size(Mesh,2)
    for j = 1:size(Mesh,2)
        x(1) = Mesh(1,i);
        x(2) = Mesh(2,j);
        
        X.hom = 0;
        X.het = 0.5;
        FE = 0.7;
        V_cell = 3.69;
        CD = j0*exp(-alpha_c*(x(2)-E0_C2H4)*const.F/(const.T*const.R));
        
        [Model1(i,j)]= Finances(X,FE,CD,Ly,x(1),V_cell,const,1);
    end 
end 


%% Figure
    figure(7)
    %Model1
    CD = j0*exp(-alpha_c*(Mesh(2,:)-E0_C2H4)*const.F/(const.T*const.R));
    s=pcolor(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model1);
    s.FaceColor = 'interp';
    set(s, 'EdgeColor', 'none');
    hold on 
    %contour(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model1,'k--','Linewidth',0.1,'ShowText','on')
    contour(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model1,'k--','Linewidth',0.1)
    colormap(map)
    caxis([0 1.2])
    ylim([5 50]);
    xlim([50 250]);
    xlabel('Current density [mA cm^{-2}]')
    ylabel('Flow rate [sccm min^{-1}]')
    pbaspect([1 1 1])

    %Model2
    figure(8)
    CD = j0*exp(-alpha_c*(Mesh(2,:)-E0_C2H4)*const.F/(const.T*const.R));
    s=pcolor(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model2);
    s.FaceColor = 'interp';
    set(s, 'EdgeColor', 'none');
    hold on 
    contour(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model2,'k--','Linewidth',0.1)
    %contour(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model2,'k--','Linewidth',0.1,'ShowText','on')
    colormap(map)
    caxis([0 1.2])
    ylim([5 50]);
    xlim([50 250]);
    xlabel('Current density [mA cm^{-2}]')
    ylabel('Flow rate [sccm min^{-1}]')
    pbaspect([1 1 1])

    figure(9)
    CD = j0*exp(-alpha_c*(Mesh(2,:)-E0_C2H4)*const.F/(const.T*const.R));
    s=pcolor(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model3);
    s.FaceColor = 'interp';
    set(s, 'EdgeColor', 'none');
    hold on 
    contour(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model3,'k--','Linewidth',0.1)
    %contour(CD*0.1,Mesh(1,:).*(10^-5*60*10^6),NPV_lim./Model3,'k--','Linewidth',0.1,'ShowText','on')
    colormap(map)
    caxis([0 1.2])
    ylim([5 50]);
    xlim([50 250]);
    xlabel('Current density [mA cm^{-2}]')
    ylabel('Flow rate [sccm min^{-1}]')
    pbaspect([1 1 1])
end
