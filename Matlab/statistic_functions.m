clear all;
close all;

xmax = 3;
ymax = 2;
zmax = 2
C = [1.8 1 1];
w = .5;
wBox = [.5 .25 .25];
Q = [1 1 1];

x = 0:xmax/200:xmax;
y = 0:ymax/200:ymax;
z = 0:zmax/200:zmax;

%% 1D code

X = x;
PC = abs(X-C(1));
PQ = abs(X-Q(1));
QC = abs(Q(1)-C(1));
wLim = abs([w-QC w+QC]);

% Heaviside(P,Q;w) if P is less than distance w from Q
CG{1} = (w-PQ>0)/(2*w);
% HeavisideBox(P,Q;wBox) if each coordinate P_i is less than distance wBox_i from Q_i, where wBox_i is the mesh size
CG{2} = (wBox(1)-abs(X-Q(1))>0)/prod(2*wBox(1));
% HeavisideRadial(P;wLim,C) if the distance of P and Q to C=0.5*(PmaP-Pmin) differs bQ less than omega
CG{3} = (min(PC-wLim(1),wLim(2)-PC)>0)/(2*(wLim(2)-wLim(1)));
% Gaussian(P,Q;w)
CG{4} = exp(-PQ.^2/(2*w^2))/sqrt(2*pi*w^2) / (erf((xmax-Q(1))/(sqrt(2)*w))-erf((0-Q(1))/(sqrt(2)*w)))*2;

figure(1);
for i=1:4
  subplot(2,2,i);
  plot(X,CG{i});
  title(['Volume:' num2str(sum(CG{i}(:))*x(2))]);
end

%% 2D code

[X,Y] = meshgrid(x,y);
PC2 = (X-C(1)).^2 + (Y-C(2)).^2;
PQ2 = (X-Q(1)).^2 + (Y-Q(2)).^2;
QC2 = (Q(1)-C(1)).^2 + (Q(2)-C(2)).^2;
wLim = [w-sqrt(QC2) w+sqrt(QC2)].^2;

% Heaviside(P,Q;w) if P is less than distance w from Q
CG{1} = ((w^2-PQ2)>0)/(pi*w^2);
% HeavisideBox(P,Q;wBox) if each coordinate P_i is less than distance wBox_i from Q_i, where wBox_i is the mesh size
CG{2} = (min(wBox(1)-abs(X-Q(1)), wBox(2)-abs(Y-Q(2)))>0)/prod(2*wBox(1:2));
% HeavisideRadial(P;wLim,C) if the distance of P and Q to C=0.5*(PmaP-Pmin) differs bQ less than omega
CG{3} = (min(PC2-wLim(1),wLim(2)-PC2)>0)/(pi*(wLim(2)-wLim(1)));
% Gaussian(P,Q;w)
CG{4} = exp(-PQ2/(2*w^2))/(2*pi*w^2) ...
  / (erf((xmax-Q(1))/(sqrt(2)*w))-erf((0-Q(1))/(sqrt(2)*w)))*2 ...
  / (erf((ymax-Q(2))/(sqrt(2)*w))-erf((0-Q(2))/(sqrt(2)*w)))*2;

figure(2);
for i=1:4
  subplot(2,2,i);
  contourf(X,Y,CG{i});
  colormap cool
  colorbar
  axis equal
  title(['Volume:' num2str(sum(CG{i}(:))*x(2)*y(2))]);
end

%% 3D code

[X,Y,Z] = meshgrid(x,y,z);
PC2 = (X-C(1)).^2 + (Y-C(2)).^2 + (Z-C(3)).^2;
PQ2 = (X-Q(1)).^2 + (Y-Q(2)).^2 + (Z-Q(3)).^2;
QC2 = (Q(1)-C(1)).^2 + (Q(2)-C(2)).^2 + (Q(3)-C(3)).^2;
wLim = [w-sqrt(QC2) w+sqrt(QC2)].^2;

% Heaviside(P,Q;w) if P is less than distance w from Q
CG{1} = ((w^2-PQ2)>0)/(4/3*pi*w^3);
% HeavisideBox(P,Q;wBox) if each coordinate P_i is less than distance wBox_i from Q_i, where wBox_i is the mesh size
CG{2} = (min(min(wBox(1)-abs(X-Q(1)), wBox(2)-abs(Y-Q(2))), wBox(3)-abs(Z-Q(3)))>0)/prod(2*wBox);
% HeavisideRadial(P;wLim,C) if the distance of P and Q to C=0.5*(PmaP-Pmin) differs bQ less than omega
CG{3} = (min(PC2-wLim(1),wLim(2)-PC2)>0)/(4/3*pi*(wLim(2)^1.5-wLim(1)^1.5));
% Gaussian(P,Q;w)
CG{4} = exp(-PQ2/(2*w^2))/(2*pi*w^2)^1.5 ...
  / (erf((xmax-Q(1))/(sqrt(2)*w))-erf((0-Q(1))/(sqrt(2)*w)))*2 ...
  / (erf((ymax-Q(2))/(sqrt(2)*w))-erf((0-Q(2))/(sqrt(2)*w)))*2 ...
  / (erf((zmax-Q(3))/(sqrt(2)*w))-erf((0-Q(3))/(sqrt(2)*w)))*2;

figure(3);
for i=1:4
  subplot(2,2,i);
  ind = floor(length(z)/2);
  [C,h] = contourf(X(:,:,ind),Y(:,:,ind),CG{i}(:,:,ind));
  colormap cool
  colorbar
  axis equal
  title(['Volume:' num2str(sum(CG{i}(:))*x(2)*y(2)*z(2))]);
end

