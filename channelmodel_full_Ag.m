function [X,FE,y,delP,CD] = channelmodel_full_Ag(E_appl,L,v,vL,c_int,k,Henry,BV,const,H,W,y0,por,D,H_c,a)
%% Flowchannel-CatalystLayer-BoundaryLayer
% This file shows a running example of the full channel model (M3).
% In this model, the flowchannel and the catalyst layer are explicitly
% modeled, whereas the boundary layer evolving in the liquid electrolyte is
% inferred from emperical correlations, depending on the flowregime. A
% method of lines (MOL) technique is employed to solve the concentration
% field in 2D, coupled with a continuation technique to allow for solving
% the BVP in the catalyst layer as the boundary layer evolves. 

%% Prepare script
%Turn off warnings
warning('off','all')
%Deval returns a warning on every iteration. If the results seem off,
%turn on warnings to check whether bvp4c returned a warning. 

%% Prepare simulation 
%Solver for IVP in flowchannel
solver = 'Heun';        %Chose between Euler and Heun (2nd order RK)
%Set number of steps
nstep.init = 50;         %Number of steps with initial size
nstep.fin  = 100;        %Number of steps with final size 

%% Run simuation 
[y,X,FE,delP,~] = FLC_CL_BL(L,W,H,H_c,vL,v,por,D,k,Henry,a,E_appl,const,c_int,y0,nstep,solver);

%% Individual functions
function [y,X,FE,delP,pos] = FLC_CL_BL(L,W,H,H_c,vL,v,por,D,k,Henry,a,E_appl,const,c_int,y0,nstep,solver)
    %This function uses the aforementioned MOL and continuation technique to
    %solve the concentration profiles in the flow channel, catalyst layer and
    %boundary layer. 

    %To solve the BVP in the catalyst layer, a continutation technique must be
    %used. Hereby, the solver uses as initial guess the
    %previous solution. This is especially handy for the evolving boundary
    %layer, as the length scales become increasingly disparate. For this to
    %work however, the initial step size must be chosen such that the boundary
    %layer evolves slowly. Once the boundary layer is sufficiently evolved, a
    %larger step size can be chosen. 

    %Allocate step sizes 
    hmin  = (H_c/1.022)^3*vL/(H*D.HCO3);   %Min. step size s.t. BL = CL @iter 1  
    hmax  = (L-nstep.init*hmin)/nstep.fin;
    hh    = [hmin*ones(1,nstep.init) hmax*ones(1,nstep.fin)];
    ny    = size(hh,2);

    %Pre-allocate gas channel solutions 
    y.CO2(:)  = zeros(1,ny);
    y.C2H4(:) = zeros(1,ny);
    y.H2(:)   = zeros(1,ny);
    y.CO2(1)  = y0(1);

    %Pre-allocate other outcomes
    nhom(:)     = zeros(1,ny);
    nhet(:)     = zeros(1,ny);
    pHsurf      = zeros(ny,50);
    pos         = zeros(1,ny);

    oldsol = 0;          %For first iteration, in which continutation not used
    x      = hmin;       %First position in channel

    switch solver 
       case 'Euler'
        for i = 2:ny
            h = hh(i-1);
            %Solve cat. layer 
            [n,oldsol,pH] = CL_BL(H,H_c,vL,h,por,D,k,a,i,E_appl,const,y.CO2(i-1)*Henry.CO2,...
                            oldsol,x,c_int);
            %Gas phase 
            f.co2     = -(n.trans.CO2)/(v*H);
            f.c2h4    = n.trans.C2H4/(v*H);
            f.h2      = n.trans.H2/(v*H);
            y.CO2(i)  = y.CO2(i-1) + h*f.co2;
            y.C2H4(i) = y.C2H4(i-1) + h*f.c2h4;
            y.H2(i)   = y.H2(i-1) + h*f.h2; 
            if y.CO2(i) < 0              %Break if all CO2 depleted
                y.CO2(end) = y.CO2(i-1);
                y.C2H4(end) = y.C2H4(i-1);
                break 
            end 
            x = x+h;
            pos(i) = x;
            %Track other results
            nhom(i) = n.hom;
            nhet(i) = n.het.CO2;
            pHsurf(i,:) = pH(:);
        end 

        case 'Heun'
        for i = 2:ny
            % 1. Euler step 
            h = hh(i);
            [n,oldsol,pH] = CL_BL(H,H_c,vL,h,por,D,k,a,i,E_appl,const,y.CO2(i-1)*Henry.CO2,...
                         oldsol,x,c_int);
            fp.co2     = -(n.trans.CO2)/(v*H);
            fp.c2h4    = n.trans.C2H4/(v*H);
            fp.h2      = n.trans.H2/(v*H);
            yp.CO2(i)  = y.CO2(i-1) + h*fp.co2;
            yp.C2H4(i) = y.C2H4(i-1) + h*fp.c2h4;
            yp.H2(i)   = y.H2(i-1) + h*fp.h2; 
            % 2. Trapezoidal step 
            x = x+h;
            pos(i) = x;
            [n,oldsol,pH] = CL_BL(H,H_c,vL,h,por,D,k,a,i,E_appl,const,yp.CO2(i)*Henry.CO2,...
                         oldsol,x,c_int);
            f.co2     = -(n.trans.CO2)/(v*H);
            f.c2h4    = n.trans.C2H4/(v*H);
            f.h2      = n.trans.H2/(v*H);
            y.CO2(i) = y.CO2(i-1) + h/2*(fp.co2 + f.co2);
            y.C2H4(i) = y.C2H4(i-1) + h/2*(fp.c2h4 + f.c2h4);
            y.H2(i)   = y.H2(i-1) + h/2*(fp.h2 + f.h2); 
            if y.CO2(i) < 0             %Break if all CO2 depleted
                y.CO2(end) = y.CO2(i-1);
                y.C2H4(end) = y.C2H4(i-1);
                break 
            end 
            %Track other results
            nhom(i) = n.hom;
            nhet(i) = n.het.CO2;
            pHsurf(i,:) = pH(:);
        end 
    end 

    %Calculate performance metrics 
    X.tot = (y.CO2(1) - y.CO2(end))/y.CO2(1);
    X.het = (y.C2H4(end)*BV.mCO2)/y.CO2(1);
    X.hom = X.tot - X.het;

    V_punkt = v*H*W;
    F_CO2 = V_punkt*1.91e3/44;

    X.het_Kas = n.het.CO2*W*L/F_CO2;
    X.hom_Kas = n.hom*W*L/F_CO2;

    CO2m = mean(y.CO2(1:end))*Henry.CO2;
    CD = BV.iCO*CO2m/BV.CO2ref*exp(-1*(BV.CO*const.F/(const.R*const.T)*(E_appl+BV.ECO)))...
        + BV.iH2*exp(-1*(BV.H2*const.F/(const.R*const.T)*E_appl));
    %FE = y.C2H4(end)*BV.nCO*BV.mCO2*v*const.F*L/(CD*Ly);   
    FE =BV.iCO*y.CO2(end)*Henry.CO2/BV.CO2ref*exp(-1*(BV.CO*const.F/(const.R*const.T)*(E_appl+BV.ECO)))/CD;
    %Calculate pressure drop and flow regime
    dh    = 2*W*H/(H+W); 
    Re    = dh*vL/(0.893e-6);
    fd    = 64/Re;
    delP  = fd*500*vL^2/dh*L;
end 

%% Calculations in catalyst layer and boundary layer in electrolyte
function [n,oldsol,pH] = CL_BL(H,H_c,vL,h,por,D,k,a,i,E_appl,const,yco2,oldsol,x,c_int)
    %In this function the catalyst layer and boundary layer are solved. The
    %function calls the solution oldsol to calculate the new concentrations
    %profiles at a further developed boundary layer.

    %Allocate extent of the boundary layer 
    BL = 1.022*(H*D.HCO3*x/vL)^(1/3);

    %%%%%%%%%  Solve BVP in catalyst layer at new BL coordinate %%%%%%%%%%%%%%
    %Deliver the jacobian and set stringent error tolerance
    options = bvpset('FJacobian', @jac, 'RelTol',1e-8);
    %For first iteration, chose small number of collocation points. All BVPs
    %are solved with the multipoint BVP setting delivered in Matlab. 
    if i == 2
        xmesh = [linspace(0,H_c,3),linspace(H_c,BL,3)];
        yinit = zeros(1,8);
        sol = bvpinit(xmesh,yinit);
        sol = bvp4c(@(x,c,r) cat(x,c,r,D,k,a,const,E_appl,H_c,BV), @(ca,cb) cat_bc2(ca,cb), sol,options);
        oldsol = sol;
    else 
        oldsol = bvpxtend(oldsol,BL,'constant');
        sol = bvp4c(@(x,c,r) cat(x,c,r,D,k,a,const,E_appl,H_c,BV), @(ca,cb) cat_bc2(ca,cb), oldsol,options);
        oldsol = sol;
    end 
    %Get pH vector at 50 equidistant point in catalyst layer
    solval = deval(sol,linspace(0,H_c,50));
    pH     = -log10(10^-11./solval(3,:));   

    %Get reaction rates and calculate transfer to/from gas phase
    n = get_reacrate(sol,k,H_c,por,E_appl);
    n.trans.CO2 = n.hom + n.het.CO2;
    n.trans.C2H4 = n.het.CO2*1/BV.mCO2;
    n.trans.H2 = n.het.H2;
%% Subfunctions 

%%%%%%%%%%%%%%%%%%%  BVP function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function dcdx = cat(x,c,region,D,k,a,const,E_appl,H_c,BV)
        %Switch regions between catalyst and boundary layer 
        cd_CO = BV.iCO*exp(-1*(BV.CO*const.F/(const.R*const.T)*(E_appl+BV.ECO)));
        cd_H2 = BV.iH2*exp(-1*(BV.H2*const.F/(const.R*const.T)*E_appl));
        switch region %Diff-RX in catalyst layer
            case 1
                dcdx = [c(2) 
                (cd_CO*c(1)/BV.CO2ref/const.F/BV.nCO/H_c + k.f1*c(1)*c(3) - k.r1*c(5))/(D.CO2*por^1.5)  
                c(4) 
                (k.f1*c(1)*c(3) - k.r1*c(5) + k.f2*c(5)*c(3) - k.r2*c(7)  - (cd_CO*c(1)/BV.CO2ref/const.F/H_c+cd_H2/const.F/H_c))/(D.OH*por^1.5) 
                c(6)
                (-k.f1*c(1)*c(3) + k.r1*c(5) + k.f2*c(5)*c(3) - k.r2*c(7))/(D.HCO3*por^1.5) 
                c(8) 
                (-k.f2*c(5)*c(3) + k.r2*c(7))/(D.CO3*por^1.5)];
            case 2 %Only diffusion in boundary layer
                dcdx = [c(2);0;c(4);0;c(6);0;c(8);0];
        end 
    end 

    %%%%%%%%%    Boundary conditions overall systems   %%%%%%%%%%%%%%%
    function res = cat_bc2(ca,cb)
        coh = c_int.OH;
        chco3 = c_int.HCO3;
        cco3 = c_int.CO3;
        res = [ca(1,1) - yco2;           %CO2
            cb(1,1) - ca(1,2)
            cb(2,1) - ca(2,2)
            cb(2,end)
            %%%
            ca(4,1)                 %OH
            cb(3,1) - ca(3,2)
            cb(4,1) - ca(4,2)
            cb(3,end) - coh
            %%%%
            ca(6,1)                 %HCO3
            cb(5,1) - ca(5,2)
            cb(6,1) - ca(6,2)
            cb(5,end) - chco3
            %%%
            ca(8,1)                 %CO3
            cb(7,1) - ca(7,2)
            cb(8,1) - ca(8,2)
            cb(7,end) - cco3];
    end 

%%%%%%%%%%%   Getting reaction rates from solution   %%%%%%%%%%%%%%%%%%%%
    function n = get_reacrate(sol,k,H_c,por,E_appl)
        %In this subfunction the concentration profiles are use to calculate the
        %reaction rates in the catalyst layer. For this, the concentrations in the
        %catalyst layer are averaged and the usual equations solved. 

        %Get solution vector
        ysub   = deval(sol,sol.x);
        loc = find(sol.x > H_c,1);    %Last point in solution vector that is in CL

        %Average species concentration in catalyst boundary
        CO2 = mean(ysub(1,1:loc));
        OH = mean(ysub(3,1:loc));
        HCO3 = mean(ysub(5,1:loc));

        %Get homogeneous reaction rate
        n.hom = (k.f1*CO2*OH - k.r1*HCO3)*H_c*por;
        if n.hom < 0
            n.hom = 0;           %CO2 from hom. RX does not enter gas phase
        end 

        %Get heterogeneours reaction rate 
        CO2m = mean(ysub(1,1:loc));
        cd_CO = BV.iCO*CO2m/BV.CO2ref*exp(-1*(BV.CO*const.F/(const.R*const.T)*(E_appl+BV.ECO)));
        cd_H2 = BV.iH2*exp(-1*(BV.H2*const.F/(const.R*const.T)*E_appl));
        n.het.CO2 = cd_CO/(BV.nCO*const.F);
        %n.het.H2  = CD*(BV.CO2ref-CO2m)/(BV.CO2ref*BV.nH2*const.F);
        n.het.H2  = cd_H2/(BV.nH2*const.F);
    end 


    function dfdy = jac(x,c,region)
        %The analytical jacobian is delivered for better accuracy in the prediction
        %of the partial derivatives. Also here, two regions need to be taken into
        %account. 
        cd_CO = BV.iCO*exp(-1*(BV.CO*const.F/(const.R*const.T)*(E_appl+BV.ECO)));

        switch region 
        case 1 
            dfdy = [0                                 1 0                               0 0                         0 0            0
                    (cd_CO/BV.CO2ref/const.F/BV.nCO/H_c- k.f1*c(3))/(D.CO2*por^1.5)  0 k.f1*c(1)/(D.CO2*por^1.5)                0 k.r1/(D.CO2*por^1.5)                0 0            0
                    0                                 0 0                               1 0                         0 0            0
                    k.f1*c(3)/(D.OH*por^1.5)                   0 (k.f1*c(1)+ k.f2*c(5))/(D.OH*por^1.5)     0 (k.f2*c(3)-k.r1)/(D.OH*por^1.5)     0 -k.r2/(D.OH*por^1.5)   0
                    0                                 0 0                               0 0                         1 0            0
                    -k.f1*c(3)/(D.HCO3*por^1.5)                 0 (k.f2*c(5)-k.f1*c(1))/(D.HCO3*por^1.5)    0 (k.r1 + k.f2*c(3))/(D.HCO3*por^1.5) 0 -k.r2/(D.HCO3*por^1.5) 1
                    0                                 0 0                               0 0                         0 0            1
                    0                                 0 (-k.f2*c(5))/(D.CO3*por^1.5)              0 -k.f2*c(3)/(D.CO3*por^1.5)          0 k.r2/(D.CO3*por^1.5)   0];
        case 2
            dfdy = [0 1 0 0 0 0 0 0 
                    0 0 0 0 0 0 0 0 
                    0 0 0 1 0 0 0 0 
                    0 0 0 0 0 0 0 0 
                    0 0 0 0 0 1 0 0 
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 1
                    0 0 0 0 0 0 0 0];
        end 
    end  
end 

end 