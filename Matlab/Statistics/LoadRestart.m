function data = LoadRestart(filename)

filename = [filename '.restart'];

fid=fopen(filename);
if (fid==-1), disp([filename ' not found']); end
stringData = textscan(fid,'%s');
fclose(fid);

fid=fopen(filename);
doubleData = textscan(fid,'%f');
fclose(fid);
% 
% %this assumes monodispersed particles
% if (length(tdata{1})==1) 
% 	%old restart files
% 	data.ParticleDensity = str2double(tdata(29));
% 	data.d = 2*str2double(tdata(28));
% 	data.Gravity = str2double(tdata(2:4))';
% 	data.N = str2double(tdata(end-3));
% 	data.Domain=str2double(tdata(5:10))';
% else
% 	%new restart version
% 	i = find(strcmp(tdata,'rho'));
% 	data.ParticleDensity = str2num(tdata{i+1});
% 	i = find(strcmp(tdata,'Particles'));
% 	data.d = 2*str2num(tdata{i+8});
% 	i = find(strcmp(tdata,'gravity'));
% 	data.Gravity = str2double(tdata(i+(1:3)))';
% 	i = find(strcmp(tdata,'Particles'));
% 	data.N = str2num(tdata{i+1});
% % 	data.Radii = str2double(tdata(i+(8:15:15*data.N)))';
% % 	data.Speed = sqrt( ... 
% %       str2double(tdata(i+(5:15:15*data.N)))'.^2 ...
% %     + str2double(tdata(i+(6:15:15*data.N)))'.^2 ...
% %     + str2double(tdata(i+(7:15:15*data.N)))'.^2);
% % 	data.Depth = str2double(tdata(i+(4:15:15*data.N)))';
% 	i = find(strcmp(tdata,'xmin'));
% 	data.Domain = str2double(tdata(i+(1:2:11)))';
% 	i = find(strcmp(tdata,'FixedParticleRadius'));
% 	data.FixedParticleRadius = str2double(tdata{i+1});
% 	i = find(strcmp(tdata,'MinInflowParticleRadius'));
% 	data.InflowParticleRadius = str2double(tdata(i+[1 3]))';
% end
% data.ParticleVolume = pi/6*data.d^3;
% data.DomainVolume = prod(data.Domain([2 4 6])-data.Domain([1 3 5]));
% data.ChuteAngle = round(atand(-data.Gravity(1)/data.Gravity(3))*400)/400;
return