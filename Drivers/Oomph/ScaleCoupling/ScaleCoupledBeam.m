nFiles = length(dir("ScaleCoupledBeamParticle_*.vtu"))

figure(1)
timeStep = 29;
dpm = readDPM(num2str(timeStep,'ScaleCoupledBeamParticle_%d.vtu'));
[~,ix] = sort(dpm.x);
plot(dpm.x(ix),dpm.vx(ix),'.-');
hold on
fem= readFEM(num2str(timeStep,'ScaleCoupledBeamFEM_%d.vtu'));
[~,ix] = sort(fem.x);
plot(fem.x(ix),fem.vx(ix),'.-');
hold off
legend('DPM','FEM')
xlabel('x [m]')
ylabel('v_x [m]')
title(dpm.time)

figure(2)
timeStep = 65;
dpm = readDPM(num2str(timeStep,'ScaleCoupledBeamParticle_%d.vtu'));
[~,ix] = sort(dpm.x);
plot(dpm.x(ix),dpm.vx(ix),'.-');
hold on
fem= readFEM(num2str(timeStep,'ScaleCoupledBeamFEM_%d.vtu'));
[~,ix] = sort(fem.x);
plot(fem.x(ix),fem.vx(ix),'.-');
hold off
legend('DPM','FEM')
xlabel('x [m]')
ylabel('v_x [m]')
title(dpm.time)


function data = readFEM(file)
% open file
fid = fopen(file);
% read time
header = textscan(fid,'%s',8,'Delimiter','\n');
npos = strfind(header{1}{2},'--');
data.time = str2double(header{1}{2}(npos(1)+8:npos(2)-1));
% read nPoints
npos = strfind(header{1}{5},'"');
nPoints = str2double(header{1}{5}(npos(1)+1:npos(2)-1));
% read positions
pos = textscan(fid,'%f %f %f',nPoints);
data.x = pos{1};
data.y = pos{2};
data.z = pos{3};
% read velocities
textscan(fid,'%s',6,'Delimiter','\n');
vel = textscan(fid,'%f %f %f',nPoints);
data.vx = vel{1};
data.vy = vel{2};
data.vz = vel{3};
% close file
fclose(fid);
end

function data = readDPM(file)
% open file
fid = fopen(file);
% read time
header = textscan(fid,'%s',7,'Delimiter','\n');
npos = strfind(header{1}{2},'--');
data.time = str2double(header{1}{2}(npos(1)+8:npos(2)-1));
% read nPoints
npos = strfind(header{1}{5},'"');
nPoints = str2double(header{1}{5}(npos(1)+1:npos(2)-1));
% read positions
pos = textscan(fid,'%f %f %f',nPoints);
data.x = pos{1};
data.y = pos{2};
data.z = pos{3};
% read velocities
textscan(fid,'%s',5,'Delimiter','\n');
vel = textscan(fid,'%f %f %f',nPoints);
data.vx = vel{1};
data.vy = vel{2};
data.vz = vel{3};
% close file
fclose(fid);
end

% freitag 8, 16 mai Di 17:00, 17 mai 13:15