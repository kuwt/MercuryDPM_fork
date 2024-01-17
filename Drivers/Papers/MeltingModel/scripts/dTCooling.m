function dTdt=dTCooling(t,T,cSolid,cLiquid,dT,R,Tm,m,L) 
if T<Tm-dT/2
    c = cSolid
elseif T>Tm+dT/2
    c = cLiquid
else
    c = (cSolid+cLiquid)/2+L/dT
end
Ta = Tm-2*dT
A = pi*R^2;
Q = 150*A*(T-Ta)+0.9*56.7e-9*A*(T^4-Ta^4);
dTdt = -Q/m/c;
end