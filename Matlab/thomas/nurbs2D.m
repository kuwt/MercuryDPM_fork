% Implementation of nurbs curves in 2D
clear variables
close all
%% Example 1: a simple line
% n control points (x,y), weights w
x = [0 1];
y = [0 0];
w = [1 1];
% m+1 knot vectors U (0=U_0<=U_i<=U_i+1<=U_m=1)
U=[0 0 1 1];
% create nurbs surface
[xC,yC] = nurbs2(x,y,w,U);
% plot nurbs surface and control points
figure(1)
plot(xC,yC,'.-',x,y,'o')
axis equal
%% Example 2: a quarter sphere
% n control points (x,y), weights w
x = [1 1 0];
y = [0 1 1];
w = [1 1 2];
% m+1 knot vectors U (0=U_0<=U_i<=U_i+1<=U_m=1)
U=[0 0 0 1 1 1];
% create nurbs surface
[xC,yC] = nurbs2(x,y,w,U);
% plot nurbs surface and control points
figure(2)
plot(xC,yC,'.-',x,y,'o')
axis equal
%% Example 3: a sphere
% n control points (x,y), weights w
x = [1 1 -1 -1 -1  1 1];
y = [0 1  1  0 -1 -1 0];
w = [2 1  1  2  1  1 2]/2;
% m+1 knot vectors U (0=U_0<=U_i<=U_i+1<=U_m=1)
U=[0 0 0 1 2 2 3 4 4 4]/4;
% create nurbs surface
[xC,yC] = nurbs2(x,y,w,U);
% plot nurbs surface and control points
figure(3)
plot(xC,yC,'.-',x,y,'o')
axis equal
%% NURBS2 returns 100 points on a nurbs curve
function [xC,yC] = nurbs2(x,y,w,U)
n = length(x);
m = length(U);
p = m-n-1;
% create 100 points on a line
u=linspace(0,max(U));
N = zeros(n,1);
% compute the corresponding points on the nurbs surface
xC = zeros(size(U));
yC = zeros(size(U));
for k=1:length(u)
    u_ = mod(u(k),1);
    %compute basis functions N of order p
    N = 1*(U(1:n)<=u_ & u_<U(2:n+1));
    for q=1:p
        N1 = (u_-U(1:n))./(U(q+(1:n))-U(1:n)).*N;
        N2 = (U(q+1+(1:n))-u_)./(U(q+1+(1:n))-U(1+(1:n))).*N([2:end 1]);
        N1(~isfinite(N1))=0;
        N2(~isfinite(N2))=0;
        N = N1+N2;
    end
    % define the points on the nurbs surface
    d = sum(N.*w);
    xC(k) = sum(N.*w.*x)/d;
    yC(k) = sum(N.*w.*y)/d;
end
end