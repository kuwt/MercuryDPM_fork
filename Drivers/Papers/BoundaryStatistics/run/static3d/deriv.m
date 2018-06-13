function dy = deriv(y,x)
dy1 = diff(y)./diff(x);
dy=(dy1([1 1:end])+dy1([1:end end]))/2;
return
