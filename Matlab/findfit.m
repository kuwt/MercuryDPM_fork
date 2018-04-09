function x=findfit(theta_stop, h_stop)
global tantheta h
tantheta = tan(theta_stop*pi/180);
h = h_stop;
%h = 3 * (.5 - tantheta)./(tantheta - .2);
x = fminsearch(@myfun, [1 .1 .2]);
%plot
plot(theta_stop, x(1) * (x(3) - tantheta)./(tantheta - x(2)));
end

function r=myfun(x)
global tantheta h
hest = x(1) * (x(3) - tantheta)./(tantheta - x(2));
r = norm(h-hest,2);
end

% function x=findfit(theta_stop, h_stop)
% global theta h
% theta = theta_stop;
% h = h_stop;
% x = fminsearch(@myfun, [0 0 0]);
% %plot
% A = x(1);
% tan1 = tan(x(2)*pi/180);
% tan2 = tan(x(3)*pi/180);
% tan_ = tan(theta*pi/180);
% plot(theta_stop, A * (tan2 - tan_)./(tan_ - tan1));
% end
% 
% function r=myfun(x)
% global theta h
% A = x(1);
% tan1 = tan(x(2)*pi/180);
% tan2 = tan(x(3)*pi/180);
% tan_ = tan(theta*pi/180);
% hest = A * (tan2 - tan_)./(tan_ - tan1);
% r = norm(h-hest);
% end
