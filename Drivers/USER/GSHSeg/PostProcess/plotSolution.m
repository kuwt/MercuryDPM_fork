function plotSolution(num)

allName=strcat('GSHCube.',num2str(num),'.stat');
smallName=strcat('GSHCube.',num2str(num),'.small.stat');
largeName=strcat('GSHCube.',num2str(num),'.large.stat');

all=loadstatistics(allName);

figure('Name','Total Temperature');
pcolor(all.x,all.z,all.Temperature);
shading flat;
colorbar;
title('Temparture');

figure('Name','Total Pressure');
pcolor(all.x,all.z,all.Pressure);
shading flat;
colorbar;
title('Total Pressure');

figure('Name','Volume Fraction');
pcolor(all.x,all.z,all.VolumeFraction);
shading flat;
colorbar;
title('Volume Fraction');

figure('Name','Kinetic Pressure');
pcolor(all.x,all.z,all.KineticPressure);
shading flat;
colorbar;
title('Kinetic Pressure');

figure('Name','Contact Stress ZZ');
pcolor(all.x,all.z,all.ContactStressZZ);
shading flat;
colorbar;
title('Contact Stress ZZ');

figure('Name','Contact Stress XZ');
pcolor(all.x,all.z,all.ContactStressXZ);
shading flat;
colorbar;
title('Contact Stress XZ');

figure('Name','Contact Stress XX');
pcolor(all.x,all.z,all.ContactStressXX);
shading flat;
colorbar;
title('Contact Stress XX');


small=loadstatistics(smallName);
large=loadstatistics(largeName);

figure('Name','Phi');
phi = (small.Density./all.Density);
pcolor(all.x,all.z,phi);
shading flat;
colorbar;
title('Phi');

keyboard;

%figure('Name','Velocity Field')
%quiver(all.x,all.z,all.Vel


