% if ~exist('data','var')
  addpath('~/code/MD/matlab/thomas')
  data=loadstatistics('*.stat');
% end
figure(1); clf;
plot(data{1}.z,deriv(data{1}.StressZZ,data{1}.z))
hold on
%plot(data{2}.z,deriv(data{2}.StressZZ,data{2}.z),'r-')
plot(data{2}.z,deriv(data{2}.StressZZ,data{2}.z)+data{2}.TractionZ,'r--')
plot(data{2}.z,data{2}.Density*data{2}.Gravity(3),'k:')
title('StressZZ vs TractionZ')

figure(2); clf;
plot(data{1}.z,deriv(data{1}.StressXZ,data{1}.z))
hold on
%plot(data{2}.z,deriv(data{2}.StressXZ,data{2}.z),'r-')
plot(data{2}.z,deriv(data{2}.StressXZ,data{2}.z)+data{2}.TractionX,'r--')
plot(data{2}.z,data{2}.Density*data{2}.Gravity(1),'k:')
title('StressXZ vs TractionX')

figure(3); clf;
plot(data{1}.z,deriv(data{1}.StressXX,data{1}.z))
hold on
plot(data{2}.z,deriv(data{2}.StressXX,data{2}.z),'r-')
plot(data{2}.z,deriv(data{2}.StressXX,data{2}.z)+data{2}.TractionZ,'r--')
plot(data{2}.z,data{2}.TractionZ,'r--')
title('StressXX vs TractionZ')
