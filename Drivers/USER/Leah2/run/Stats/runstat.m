j=0;
mymap=colormap;
mymap(1,1)=0;
mymap(1,2)=0;
mymap(1,3)=0;
colormap(mymap);
keyboard;
for i=0:0.1:3;
    
   
filename=strcat('CinderDriver.5.t=',num2str(i),'.stat');

%%dataX=loadstatistics('CinderDriver.1.X.stat')
%dataZ=loadstatistics('CinderDriver.1.Z.stat')
dataXZ=loadstatistics(filename)


%% plot volume fraction
figure(1); clf
%set(gcf,'Position',get(0,'ScreenSize').*[0 0 1 .8])

subplot(2,1,1)
pcolor(dataXZ.x,dataXZ.z,dataXZ.Nu)
shading interp;
xlabel('x')
ylabel('z')
title('Volume fraction \rho/\rho_p')
%axis equal
colorbar('Location','East')

%subplot(2,4,5:7)
%plot(dataX.x,dataX.Nu)
%xlabel('x')
%ylabel('Volume fraction \rho/\rho_p')
%axis tight

%subplot(2,4,4)
%plot(dataZ.Nu,dataZ.z)
%ylabel('z')
%xlabel('Volume fraction \rho/\rho_p')
%axis tight

%% plot speed
subplot(2,1,2)
%set(gcf,'Position',get(0,'ScreenSize').*[0 0 1 .8])

%subplot(2,4,1:3)
Speed = sqrt(dataXZ.VelocityX.^2+dataXZ.VelocityY.^2+dataXZ.VelocityZ.^2);
pcolor(dataXZ.x,dataXZ.z,Speed)
shading interp;
xlabel('x')
ylabel('z')
title('Speed |u|')
%axis equal
colorbar('Location','East')

%subplot(2,4,5:7)
%Speed = sqrt(dataX.VelocityX.^2+dataX.VelocityY.^2+dataX.VelocityZ.^2);
%plot(dataX.x,Speed)
%xlabel('x attempt')
%ylabel('Speed |u|')
%axis tight

%subplot(2,4,4)
%Speed = sqrt(dataZ.VelocityX.^2+dataZ.VelocityY.^2+dataZ.VelocityZ.^2);
%plot(Speed,dataZ.z)
%ylabel('z')
%xlabel('Speed |u|')
%axis tight


%s=size(Speed);
%s=s(1);
%for i=1:s,
 %   if (Speed(i) > 0)
  %      fprintf('Max Position XZ: %d \n', dataX.x(i));
   % end
%end

% fprintf('Max Position X: %d \n', max(dataX.x));
j=j+1;
M(j)=getframe(gcf);

end

movieview(M);
movie2avi(M,'CinderMovie.avi');
