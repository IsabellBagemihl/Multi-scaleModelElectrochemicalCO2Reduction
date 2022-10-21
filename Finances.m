function [NPV]=Finances(X,FE,CD,Ly,v,Ec,const,model)
    Data;
    SetupBest;
    
    %Base values
    ProdRate    = 10000;        % kg d^-1
    lifetime    = 20;           % years
    optime      = 350;          % days per year
    IR          = 0.1;          % Interest rate
    
    %ConversionValues
    hours_day       = 24;             %hours per day
    minutes_hour    = 60;          %minutes per hour
    seconds_minute  = 60;        %seconds per minute
    
    %%%%%%%%%%%%%%Electrolyser Scale%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_dot_C2H4_target   = ProdRate*optime;                              %Yearly production target in kg C2H4 per year
    F_dot_C2H4_target   = ProdRate/(MC2H4/1000)...
                        /hours_day/minutes_hour/seconds_minute;         %Daily production target in moles C2H4 per day
    A_r                 = F_dot_C2H4_target*12*const.F/(CD*FE);         %Electrolyzer area in m^2
    V_dot_r             = F_dot_C2H4_target*2/X.het/y0(1);              %Total flowrate through electrolyzer in m^3 per second
    m_dot_CO2           = (MCO2/MC2H4*2 + X.hom/X.het)*ProdRate*optime; %Annual CO2 consumption in kg per year
    
    %Sanity check
    n       = A_r/(Lw*Ly);                   %Channel number
    Fl.tot2 = n*v*Lw*L*y0(1);          %this should equal F_dot_C2H4_target
    
    %Overpotential
    if model > 1
        eta.actA    = const.R*const.T/(0.5*const.F)*asinh(CD/(2*1e-7));
        eta.ohm     = CD*(L/sigma_el+Lm/sigma_m);
        eta.tot     = abs(Ec) + eta.actA + eta.ohm + 1.23;
    else
        eta.tot     = Ec;
    end

    P_r     = A_r*CD*eta.tot;               %Power consumption in Watt
    Power   = P_r*10^-3*optime*hours_day;   %Power consumption in kWh per year

    %%%%%%%%%%%%%%%%%%Process scale%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Total capital investment costs
    Price.TCI  = (A_r*Price.A)*(1+35/65) + (Price.PSA*(V_dot_r*3.6)^scale);        
    %Revenue 
    Price.Rev  = m_dot_C2H4_target*Price.Prod;
    %Operating costs
    Price.Op   = m_dot_CO2*Price.CO2 + Price.Elec*(Power + V_dot_r*7.56*10^6);
    %Maintenance costs
    Price.m    = 0.025*A_r*Price.A;
    
    CF = zeros(1,lifetime);
    for i = 1:lifetime
        CF(i) = (Price.Rev - Price.Op - Price.m)/(1 + IR)^i; 
    end 
    NPV = -Price.TCI + sum(CF(:)); 
end
