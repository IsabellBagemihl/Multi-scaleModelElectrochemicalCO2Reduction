function [X,FE,delP] = channelmodel_simp(Ly,v, F,y0,CD,CDhom,L)
%This function solves the simplistic channel model (M2). This entails the 
%assumption of constant homogeneous consumption of 50 mA/cm^2, captured in the variable 
%CDhom.

c = solvesystem(v,F,L,Ly,y0,CD,CDhom);
X.tot = (y0(1) - c(end,1))/y0(1);
X.het = (c(end,2)*2)/y0(1);
X.hom = X.tot-X.het;
FE = c(end,2)*12*v*F*L/(CD*Ly);
delP = 4500;                              %Assuming operation at delP limit

function c = solvesystem(v,F,L,Ly,y0,CD,CDhom)
[~,c] = ode15s(@odefun,[0 Ly], [y0(1) 0]);
function dcdx = odefun(~,c)
dcdx(1) = -CD/(6*F*v*L)*c(1)/y0(1)- CDhom/(6*F*v*L);
%Ethylene
dcdx(2) = CD/(12*F*v*L)*c(1)/y0(1);
dcdx = dcdx';
end 
end
end