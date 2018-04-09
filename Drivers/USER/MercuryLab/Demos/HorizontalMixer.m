%cd ~/Mercury/CleanTrunk/cmake-build-debug/Drivers/USER/MercuryLab/Demos/
cd ~/MercuryDPM/Trunk/Build/Drivers/USER/MercuryLab/Demos/
%% read in first file, to get the initial positions
f = fopen('HorizontalMixerParticle_0.vtu');
% header
line = textscan(f,'%s %s %s %s %s %s %s',1,'Delimiter','\n');
% number of particles
N = textscan(line{5}{1}(24:end),'%d',1); N=N{1};
% positions
P = textscan(f,'%f %f %f',N);
%scatter(P{1},P{2})
fclose(f);
%% define a new speciesIndex, based on position, to color particles
index = (P{1}>mean(P{1})) + 2*(P{2}>mean(P{2}));
%% read in second file, a write out again with modified index
%f = fopen('bkup/HorizontalMixerDemoParticle_8.vtu');
f = fopen('HorizontalMixerParticle_499.vtu');
g = fopen('mod.vtu','w');
% header
line = textscan(f,'%s',3*N+15,'Delimiter','\n');
for i=1:length(line{1}), fprintf(g,'%s\n',line{1}{i}); end
% i/o indSpecies
textscan(f,'%f',N,'Delimiter','\n');
fprintf(g,'%f\n',index);
% footer
line = textscan(f,'%s','Delimiter','\n');
for i=1:length(line{1}), fprintf(g,'%s\n',line{1}{i}); end
fclose(f);
fclose(g);
