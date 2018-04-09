function same_scale(plots,ynum,xnum)
%use same scale

subplot(ynum,xnum,plots(1));
axis tight
a = axis;
for i = plots(2:end)
  subplot(ynum,xnum,i);
  axis tight
  a_local = axis;
  if (a_local(3)<a(3)); a(3) = a_local(3); end
  if (a_local(4)>a(4)); a(4) = a_local(4); end
end
for i = plots
  subplot(ynum,xnum,i);
  axis(a);
end

return