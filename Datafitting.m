clc;
clear;
close all;
old = true;
npoints = 20; %Number of datapoints for the plotted U-I curve 

Data;

mdl = @(b,x)(-b(1)*exp(-b(2)*(x-(E0_C2H4))*const.F/(const.T*const.R)));
y = [-1000,-2000,-3000];
beta0 = [1e-8,1];


switch old 
%Exact values found by PlotDigitizer (https://plotdigitizer.com/app) from
%Fig. 4A rightmost panel in Tan2020
%(https://doi.org/10.1016/j.joule.2020.03.013)
%X = [-0.7607 -0.818 -0.8775]
    case true 
        opts = statset('nlinfit');
        opts.RobustWgtFun = 'bisquare';
         X = [-0.775 -0.825 -0.875];          %Applied voltage
         [beta,R,J,CovB,MSE] = nlinfit(X,y,mdl,beta0,opts);
    case false 
         E0 = ones(1,3)*0.08;
         X  = [-0.7607 -0.818 -0.8775] + E0;  %Plus or minus??
         [beta,R,J,CovB,MSE] = nlinfit(X,y,mdl,beta0);      
end 


ci = nlparci(beta,R,'Jacobian',J);
Eapp = linspace(-0.4,-0.9,npoints);
[ypred,delta] = nlpredci(mdl, Eapp, beta, R,'Jacobian',J,'MSE',MSE,'SimOpt','on');
lower = ypred - delta;
upper = ypred + delta;
figure()
hold on
plot(Eapp,ypred,'b')
plot(Eapp,[upper;lower], 'r--')
for i = 1:3
    scatter(X(i),y(i),'filled','dr')
end 


for i = 1:npoints 
    cd(i) = -beta(1)*exp(-beta(2)*(Eapp(i)-(E0_C2H4))*const.F/(const.T*const.R));
end 
plot(Eapp,cd)
%ylim([-4000 0])
xlim([-0.9 -0.4])
hold on
for i = 1:3
    scatter(X(i),y(i),'filled','dr')
end 