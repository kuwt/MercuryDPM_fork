%% path of the loadstatistics.m file
addpath('~/Mercury/Trunk/Matlab/thomas')
%% load statistics for full flow and species 0 and 1 for x-and y-mixing
T=loadstatistics('TriolietDietFeeder.stat');
Y0=loadstatistics('TriolietDietFeeder.0.stat');
Y1=loadstatistics('TriolietDietFeeder.1.stat');
X0=loadstatistics('TriolietDietFeeder.2.stat');
X1=loadstatistics('TriolietDietFeeder.3.stat');
%% load high-res time averaged data
H=loadstatistics('TriolietDietFeeder.timeAveraged.stat');
%% load torque data 
raw=importdata('TriolietDietFeeder.torques');
M.time = raw.data(:,1);
M.torqueLeft = raw.data(:,2:4);
M.torqueRight = raw.data(:,5:7);

%% compute mixing in x/y direction at each time step
figure(1)
t=zeros(1,length(T));
phiY=zeros(1,length(T));
phiX=zeros(1,length(T));
for i=1:length(T)
    t(i)=T{i}.time(1);
    n=T{i}.Nu(:);

    x0=Y0{i}.Nu(:)./n;
    x1=Y1{i}.Nu(:)./n;
    s=x0.*log(x0)+x1.*log(x1);
    s(n==0|x0==0|x1==0)=0;
    S=mean(s.*n);
    S_mix = log(0.5)*mean(n);
    S_seg = 0;
    phiY(i)=(S-S_seg)/(S_mix-S_seg);

    x0=X0{i}.Nu(:)./n;
    x1=X1{i}.Nu(:)./n;
    s=x0.*log(x0)+x1.*log(x1);
    s(n==0|x0==0|x1==0)=0;
    S=mean(s.*n);
    S_mix = log(0.5)*mean(n);
    S_seg = 0;
    phiX(i)=(S-S_seg)/(S_mix-S_seg);
end
plot(t,phiY,'x',t,phiX,'o')
legend('mixing in y-direction','mixing in x-direction')
xlabel('time [s]')
ylabel('Mixing index')
ylim([0 1])

%% plot the final mixing state
figure(2)
subplot(2,1,1)
val=Y0{end};
contourf(val.x,val.y,val.Nu,20,'EdgeColor','none')
axis image
title('Mixing in y-direction')
subplot(2,1,2)
val=X0{end};
contourf(val.x,val.y,val.Nu,20,'EdgeColor','none')
axis image
title('Mixing in x-direction')

%%
figure(3)
plot(M.time,M.torqueLeft(:,3),M.time,M.torqueRight(:,3))
legend('left screw','right screw')
xlabel('time [s]')
ylabel('Torque [Nm]')

%%
figure(4)
contourf(H.x,H.y,H.Temperature,20,'EdgeColor','none')
axis image
colorbar
title('velocity fluctuations [(m/s)^2] (averaged over height)')

%%
figure(5)
contourf(H.x,H.y,H.Density,20,'EdgeColor','none')
axis image
colorbar
title('density [kg/m^3] (averaged over height)')

%%
figure(6)
H.Speed = sqrt(H.MomentumX.^2+H.MomentumY.^2+H.MomentumZ.^2)./H.Density;
H.Speed(~isfinite(H.Speed))=0;
contourf(H.x,H.y,H.Speed,20,'EdgeColor','none')
h=colormap('jet');
colormap(0.5+0.5*h)
hold on
ix=1:2:size(H.x,1);
iy=1:2:size(H.x,2);
quiver(H.x(ix,iy),H.y(ix,iy),H.VelocityX(ix,iy),H.VelocityY(ix,iy),1.5,'Color','k')
hold off
axis image
title('velocity [m/s] (averaged over height)')
