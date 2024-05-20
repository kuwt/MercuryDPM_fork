Rd=1;
% (a_n+1 - a_n) = da = Rd/2(a+da)*dt
% -> da^2+a*da-Rd/2 = 0
% -> da = -a/2+sqrt(a^2/4+Rd/2*dt)
t = linspace(0,10,100);
dt = t(2)-t(1);
da = @(a) -a/2+sqrt((a/2)^2+Rd*dt/2);
a = zeros(size(t));
for i=1:length(t)-1
  a(i+1) = a(i)+da(a(i));
end
da2 = @(t,a) Rd/(2*a);
[t2,a2] = ode45(da2,[0,10],1e-13)
plot(t,a,'.',t2,a2,'x',t,sqrt(Rd*t))