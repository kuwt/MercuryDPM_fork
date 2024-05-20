format long
format compact
Sf = 1e6;
R = 50e-6;
rho = 1050*Sf;
m = rho*(4/3*pi*R^3);
E = 2.36e9;
poisson = 0.4;
Gs = 0.5*E/(1+poisson);
Es = 0.5*E/(1-poisson);
zeta = 0.9;
sigma = 34.4e-3*sqrt(Sf);
c = 1200/Sf;
k = 0.12;
h = 150;
eps = 0.9;
sigB = 56.7e-9;
T = 451.15;
Ea = 32.5e3;
Rg = 8.31432;
w = 150e-6;
alpha = 0.33;
P = 2.3;
eta = 0.01*exp(Ea/Rg/T)*sqrt(Sf); %
te = sqrt(m/Es/R)
td = sqrt(m/zeta^2/Es/R) %crit damped = min
ts = sqrt(m/sigma)
tv = m/eta/R
tsv = eta*R/sigma
tcd = m*c/k/R
tcv = m*c/h/R/R
tr = m*c/eps/sigB/R/R/T/T/T
th = m*c*w*w*T/alpha/R/R/P
tRayleigh = pi*R*sqrt(rho/Gs)/(0.1631*poisson+0.8766);
dt = 0.1*tRayleigh