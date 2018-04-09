addpath('~/code/MD/matlab/thomas')

% ene=cell(1,4);
% ene{1}=load_ene('Test/ContinuousTiltMu0.1Mur0.02.ene');
% ene{2}=load_ene('Test/ContinuousTiltMu0.1Mur0.ene');
% ene{3}=load_ene('Test/SteppedTiltMu0.1Mur0.02.ene');
% ene{4}=load_ene('Test/SteppedTiltMu0.1Mur0.ene');
% 
c=lines(4)

figure(1); clf
for i=1:4
  semilogy(ene{i}.field9(10:end),smooth(ene{i}.ene_kin(10:end)),'.','Color',c(i,:),'Displayname',ene{i}.name)
  hold on
end
legend('show')

figure(2); clf
for i=1:4
  plot(ene{i}.t(10:end),ene{i}.field9(10:end),'.','Color',c(i,:),'Displayname',ene{i}.name)
  hold on
end
legend('show')

figure(3); clf
for i=1:4
  plot(ene{i}.field9(10:end),ene{i}.X_COM(10:end),'.','Color',c(i,:),'Displayname',ene{i}.name)
  hold on
end
xlim([20,30])
legend('show','Location','Best')
