function int = int2inf(y,x)
int1 = (cumsum(y)-sum(y))*diff(x(1:2));
return
