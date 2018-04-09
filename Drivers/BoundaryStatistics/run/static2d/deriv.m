function dydx = deriv(y,x)

if size(x,2)==1|(x(1,1)==x(1,2)) 
  if x(1,1)==x(2,1)
    dydx = zeros(size(x));
  else
    %vertical differencing
    dydx = diff(y,1,1)./diff(x,1,1);
    dydx= (dydx([1:end end],:)+dydx([1 1:end],:))/2;
  end
else
  %horizontal differencing
  dydx = diff(y,1,2)./diff(x,1,2);
  dydx= (dydx(:,[1:end end])+dydx(:,[1 1:end]))/2;
end
return
