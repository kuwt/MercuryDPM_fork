% Implementation of nurbs surfaces in 3D
clear variables
close all
%% Example 1: a quarter cylinder
% n control points (x,y,z), weights w
x = [1 1 0;1 1 0]';
y = [0 1 1;0 1 1]';
z = [0 0 0;2 2 2]';
w = [1 1 2;1 1 2]';
% m+1 knot vectors
U=[0 0 0 1 1 1];
V=[0 0 1 1];
% create nurbs surface
[xC,yC,zC] = nurbs3(x,y,z,w,U,V);
% plot nurbs surface and control points
figure(1)
surf(xC,yC,zC)
hold on
plot3(x,y,z,'o')
hold off
axis equal
%% Example 2: a sphere
% n control points (x,y,z), weights w
x = [0 0 0 0 0 0 0
     2 2 -2 -2 -2 2 2
     2 2 -2 -2 -2 2 2
     0 0 0 0 0 0 0];
y = [0 0 0 0 0 0 0
     0 4 4 0 -4 -4 0
     0 4 4 0 -4 -4 0
     0 0 0 0 0 0 0];
z = [1 1 1 1 1 1 1
     1 1 1 1 1 1 1
     -1 -1 -1 -1 -1 -1 -1
     -1 -1 -1 -1 -1 -1 -1];
w = [9 3 3 9 3 3 9
     3 1 1 3 1 1 3 
     3 1 1 3 1 1 3
     9 3 3 9 3 3 9]/9;
% m+1 knot vectors
U=[0 0 0 0 1 1 1 1];
V=[0 0 0 0 1 1 1 2 2 2 2]/2;
% create nurbs surface
[xC,yC,zC] = nurbs3(x,y,z,w,U,V);
% plot nurbs surface and control points
figure(2)
surf(xC,yC,zC)
hold on
plot3(x,y,z,'o')
hold off
axis equal
%% NURBS3 returns 100x100 points on a nurbs surface
function [xC,yC,zC] = nurbs3(x,y,z,w,U,V)
% n control points (x,y), weights w
[nu,nv] = size(x);
% m+1 knot vectors U (0=U_0<=U_i<=U_i+1<=U_m=1)
mu=length(U);
mv=length(V);
% polynomial degree p (m=n+(p+1))
pu = mu-nu-1;
pv = mv-nv-1;
% create 100x100 grid of points (uu,vv) in a square
u=linspace(0,max(U));
v=linspace(0,max(V));
[uu,vv]=meshgrid(u,v);
% compute the corresponding points on the nurbs surface
xC = zeros(size(uu));
yC = zeros(size(uu));
zC = zeros(size(uu));
for k=1:length(uu(:))
    u_ = uu(k);
    v_ = vv(k);
    %compute basis functions Nu(u) of order pu
    Nu = 1*(U(1:nu)<=u_ & u_<U(2:nu+1));
    for q=1:pu
        N1 = (u_-U(1:nu))./(U(q+(1:nu))-U(1:nu)).*Nu;
        N2 = (U(q+1+(1:nu))-u_)./(U(q+1+(1:nu))-U(1+(1:nu))).*Nu([2:end 1]);
        N1(~isfinite(N1))=0;
        N2(~isfinite(N2))=0;
        Nu = N1+N2;
    end
    %compute basis functions Nv(v) of order pv
    Nv = 1*(V(1:nv)<=v_ & v_<V(2:nv+1));
    for q=1:pv
        N1 = (v_-V(1:nv))./(V(q+(1:nv))-vv(1:nv)).*Nv;
        N2 = (V(q+1+(1:nv))-v_)./(V(q+1+(1:nv))-V(1+(1:nv))).*Nv([2:end 1]);
        N1(~isfinite(N1))=0;
        N2(~isfinite(N2))=0;
        Nv = N1+N2;
    end
    %compute weighted tensorial basis functions N(u,v)  
    N = (Nu'*Nv).*w;
    % define the points on the nurbs surface
    d = sum(sum(N));
    xC(k) = sum(sum(N.*x))/d;
    yC(k) = sum(sum(N.*y))/d;
    zC(k) = sum(sum(N.*z))/d;
end
end