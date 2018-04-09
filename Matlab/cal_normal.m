syms xa ya za
syms xb yb zb

syms theta

A=[cos(theta)*xa-sin(theta)*za;ya;sin(theta)*xa+cos(theta)*za];

B=[cos(theta)*xb-sin(theta)*zb;yb;sin(theta)*xb+cos(theta)*zb];

D=B-A;

C=[cos(theta) 0 -sin(theta)];

E=cross(C,D)

dot(E,D)
dot(E,C)


keyboard;