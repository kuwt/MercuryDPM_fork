function [m,b]=findlinearfit(x_, y_)
% finds linear fit of type y_est = m*x+b
global x y
x = x_;
y = y_;
coeff = fminsearch(@myfun, [0 0]);
m = coeff(1);
b = coeff(2);
end

function r=myfun(coeff)
global x y
yest = coeff(1) * x + coeff(2);
r = norm(y-yest,2);
end