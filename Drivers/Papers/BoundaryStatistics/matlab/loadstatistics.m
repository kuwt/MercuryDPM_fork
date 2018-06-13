% loadstatistics (version 1)
%
% loads data from a .stat or .data file
%
% it accepts a single filename, a cell of filenames, or a string that will
% be inserted into ls
%
% Note for .stat files: should work for any stattype
%
% 
function data=loadstatistics(filenames)

if iscell(filenames) % if argument is a cell of filenames
	data=cell(size(filenames));
  for i=1:length(filenames)
		data{i} = load_file(filenames{i});
	end
else 
  if (isempty(strfind(filenames,'?'))&&isempty(strfind(filenames,'*'))) % if argument is a single filename
  	data = load_file(filenames);
  else % if argument contains * or ?, run it through ls to procuce a cell of filenames
    data = loadstatistics(strread(ls([filenames ' -1']),'%s ')); 
  end
end

return

%load a single .stat or .data file
function data = load_file(filename)
disp(['loading data from ' filename]);

% distinguishes between stat and data files
if (strcmp(filename(end-4:end),'.data'))
	data = load_data_file(filename);
else
	data = load_stat_file(filename);
end

% \todo{why is this needed?}
%data.label = data.name;

data = make_readable(data);

data = get_standard_variables(data);
if (data.nz>1), data = get_depth_variables(data); end

%[data.MomentumEquationsRemainder, data.MomentumEquationsMaximum]=get_momentum_equation(data);

return

%loads all variables stored in a stat file
function data = load_stat_file(filename)
data.name = filename(1:end-5);

% load raw data from file, with the first three lines as header files
% (such that matlab recognizes the array size)
rawdata = importdata(filename,' ',3);

if iscell(rawdata), return; end

if (size(rawdata.data,2)==26)
  % for old .stat files
  disp('WARNING:outdated stat file')
  rawdata.data = [ ...
    rawdata.data(:,1:14) ...
    zeros(size(rawdata.data,1),3) ...
    rawdata.data(:,15:end) ...
    zeros(size(rawdata.data,1),10)];
  data.variable_names = { ...
    'Nu'; 'Density'; ...
    'MomentumX'; 'MomentumY'; 'MomentumZ'; ...
    'MomentumFluxXX'; 'MomentumFluxXY'; 'MomentumFluxXZ'; 'MomentumFluxYY'; 'MomentumFluxYZ'; 'MomentumFluxZZ'; ...
    'EnergyFluxX'; 'EnergyFluxY'; 'EnergyFluxZ'; ...
    'NormalStressXX'; 'NormalStressXY'; 'NormalStressXZ'; 'NormalStressYY'; 'NormalStressYZ'; 'NormalStressZZ'; ...
    'TangentialStressXX'; 'TangentialStressXY'; 'TangentialStressXZ'; 'TangentialStressYY'; 'TangentialStressYZ'; 'TangentialStressZZ'; ...
    'FabricXX'; 'FabricXY'; 'FabricXZ'; 'FabricYY'; 'FabricYZ'; 'FabricZZ'; ...
    'CollisionalHeatFluxX'; 'CollisionalHeatFluxY'; 'CollisionalHeatFluxZ'; ...
    'Dissipation'; ...
    };
else
  % reads variable names from line 1 and text output from line 2
  variable_names = textscan(rawdata.textdata{1},'%s ');
  data.variable_names = variable_names{1};
end

%also allows input from the restart and ene files if they exist
data = read_restart(filename,data);
data = read_ene(filename,data);

%reads header file of stat file
text = textscan(rawdata.textdata{2},'%s ');
data.text = text{1};

% writes the raw data into bits of data for each timestep
index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
data.time = cell2mat(textscan(rawdata.textdata{3},'%f '));
data.w = cell2mat(textscan(rawdata.textdata{2}(3:end),'%f '));

data.coordinates = rawdata.data(1:index_time(1)-1,1:3);
data.variables = rawdata.data(1:index_time(1)-1,4:end);
if (length(index_time)>1)
  data.variance = rawdata.data(index_time(1)+1:index_time(2)-1,4:end);
end

% \todo{why is this needed?}
%data.description = ['Goldhirsch w/d=' num2str(str2double(data.text{2})/1e-3)];

return

%loads from data file the same variables that are present in the
%stat file (for comparison to Stefan's code)
function data = load_data_file(filename)
data.name = filename(1:end-5);

% load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',17);

% reads variable names from line 1 and text output
% from line 2
data.variable_names = { ...
	'Nu'; 'Density'; ...
	'MomentumX'; 'MomentumY'; 'MomentumZ'; ...
	'MomentumFluxXX'; 'MomentumFluxXY'; 'MomentumFluxXZ'; 'MomentumFluxYY'; 'MomentumFluxYZ'; 'MomentumFluxZZ'; ...
	'EnergyFluxX'; 'EnergyFluxY'; 'EnergyFluxZ'; ...
	'NormalStressXX'; 'NormalStressXY'; 'NormalStressXZ'; 'NormalStressYY'; 'NormalStressYZ'; 'NormalStressZZ'; ...
	'TangentialStressXX'; 'TangentialStressXY'; 'TangentialStressXZ'; 'TangentialStressYY'; 'TangentialStressYZ'; 'TangentialStressZZ'; ...
	'FabricXX'; 'FabricXY'; 'FabricXZ'; 'FabricYY'; 'FabricYZ'; 'FabricZZ'; ...
	'CollisionalHeatFluxX'; 'CollisionalHeatFluxY'; 'CollisionalHeatFluxZ'; ...
	'Dissipation'; ...
	};
text = textscan([rawdata.textdata{8} ' ' rawdata.textdata{9}],'%s ');
data.text = [text{1}{5} num2str(str2double(text{1}{9})*100)];
rho_text = textscan(rawdata.textdata{11},'%s ');
rho = str2double(rho_text{1}{end});

% writes the raw data into bits of data for each timestep
%index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
data.time = 0;
data.coordinates = rawdata.data(:,4:6);
nu = max(rawdata.data(:,10),1e-200*ones(size(rawdata.data,1),1));
data.variables = [nu ...
	nu*rho ...
	rawdata.data(:,51:53).*nu(:,[1 1 1])*rho ...
	rawdata.data(:,[41:43 45 46 49]) ...
	nan(size(rawdata.data,1),3) ...
	rawdata.data(:,[21:23 25 26 29]) ...
	rawdata.data(:,[31:33 35 36 39]) ...
	rawdata.data(:,[71:73 75 76 79]) ...
	nan(size(rawdata.data,1),4) ...
	];

%data.description = ['Luding ' data.text];

return

% rewrites coordinates and variables into n dimensional shape where n is
% the number of dimensions in stattype (easy for plotting) 
function [data] = make_readable(data)

data.nx=length(data.coordinates(:,1))/sum(data.coordinates(:,1)==data.coordinates(1,1));
data.ny=length(data.coordinates(:,2))/sum(data.coordinates(:,2)==data.coordinates(1,2));
data.nz=length(data.coordinates(:,3))/sum(data.coordinates(:,3)==data.coordinates(1,3));
shape=size(squeeze(zeros(data.nz,data.ny,data.nx)));
data.x = reshape(data.coordinates(:,1),shape);
data.y = reshape(data.coordinates(:,2),shape);
data.z = reshape(data.coordinates(:,3),shape);
data = rmfield(data,{'coordinates'});


if isfield(data,'variance')
	for i=1:length(data.variable_names)
  	data.(['variance_' sscanf(data.variable_names{i},'%s ',16)]) ...
      =reshape(data.variance(:,i),shape);
	end
	data = rmfield(data,{'variance'});
end

for i=1:length(data.variable_names)
	data.(sscanf(data.variable_names{i},'%s ',16)) ...
    = reshape(data.variables(:,i),shape);
end
data = rmfield(data,{'variables','variable_names'});

return

% extracts a load of extra variables (such as stress) from the basic
% microscopic fields 
function data = get_standard_variables(data)

data.MomentumX(isnan(data.MomentumX)) = 0;
data.MomentumY(isnan(data.MomentumY)) = 0;
data.MomentumZ(isnan(data.MomentumZ)) = 0;

data.VelocityX = data.MomentumX./data.Density;
data.VelocityY = data.MomentumY./data.Density;
data.VelocityZ = data.MomentumZ./data.Density;

data.VelocityX(isnan(data.VelocityX)) = 0;
data.VelocityY(isnan(data.VelocityY)) = 0;
data.VelocityZ(isnan(data.VelocityZ)) = 0;

data.TractionX = data.NormalTractionX + data.TangentialTractionX;
data.TractionY = data.NormalTractionY + data.TangentialTractionY;
data.TractionZ = data.NormalTractionZ + data.TangentialTractionZ;

data.KineticStressXX = data.MomentumFluxXX - data.Density.*data.VelocityX.*data.VelocityX;
data.KineticStressXY = data.MomentumFluxXY - data.Density.*data.VelocityX.*data.VelocityY;
data.KineticStressXZ = data.MomentumFluxXZ - data.Density.*data.VelocityX.*data.VelocityZ;
data.KineticStressYY = data.MomentumFluxYY - data.Density.*data.VelocityY.*data.VelocityY;
data.KineticStressYZ = data.MomentumFluxYZ - data.Density.*data.VelocityY.*data.VelocityZ;
data.KineticStressZZ = data.MomentumFluxZZ - data.Density.*data.VelocityZ.*data.VelocityZ;
data.StressXX = data.NormalStressXX + data.TangentialStressXX + data.KineticStressXX;
data.StressXY = data.NormalStressXY + data.TangentialStressXY + data.KineticStressXY; 
data.StressXZ = data.NormalStressXZ + data.TangentialStressXZ + data.KineticStressXZ;
data.StressYX = data.NormalStressYX + data.TangentialStressYX + data.KineticStressXY;
data.StressYY = data.NormalStressYY + data.TangentialStressYY + data.KineticStressYY;
data.StressYZ = data.NormalStressYZ + data.TangentialStressYZ + data.KineticStressYZ;
data.StressZX = data.NormalStressZX + data.TangentialStressZX + data.KineticStressXZ;
data.StressZY = data.NormalStressZY + data.TangentialStressZY + data.KineticStressYZ;
data.StressZZ = data.NormalStressZZ + data.TangentialStressZZ + data.KineticStressZZ;
data.Temperature = (...
	(data.MomentumFluxXX + data.MomentumFluxYY + data.MomentumFluxZZ)./data.Density ...
	- (data.VelocityX.^2 + data.VelocityY.^2 + data.VelocityZ.^2) )/3;
data.Pressure =  (data.StressXX + data.StressYY + data.StressZZ)/3;

data.CoordinationNumber=(data.FabricXX+data.FabricYY+data.FabricZZ)./data.Nu;
%data.TotalCoordinationNumber=sum(data.FabricXX+data.FabricYY+data.FabricZZ)/sum(data.Nu(:));

% data.ShearRate = diff(data.VelocityX)./diff(data.z);
% data.ShearRate(end+1) = data.ShearRate(end);
% data.InertialParameter=data.ShearRate*data.d./sqrt(data.StressZZ/data.ParticleDensity);

%macroscopic friction coefficient, or deviator-stress ratio sD1:=(S1-S3)/(2*p)  or sD2 (to be discussed)
data.MacroFrictionCoefficient = zeros(size(data.x));
data.Sd = zeros(size(data.x));
data.SdStar = zeros(size(data.x));
for j=1:length(data.x(:))
	Stress = [ ...
		data.StressXX(j) data.StressXY(j) data.StressXZ(j); ...
		data.StressXY(j) data.StressYY(j) data.StressYZ(j); ...
		data.StressXZ(j) data.StressYZ(j) data.StressZZ(j) ];
	if (isequal(Stress,zeros(3))|| sum(sum(isnan(Stress))) )
		data.MacroFrictionCoefficient(j) = NaN;
		data.Sd(j) = nan;
		data.SdStar(j) = nan;
	else
		[v1,d1]=eig(Stress);
		[v,d]=dsort_ev(v1,d1);
    d = diag(d);
    p=sum(d)/3;
		data.MacroFrictionCoefficient(j) = -data.StressXZ(j)/data.StressZZ(j);
		data.Sd(j) = (d(1)-d(3))/(2*p);
		data.SdStar(j) = sqrt(((d(1)-d(3))^2+(d(2)-d(3))^2+(d(1)-d(2))^2) / (6*p^2));
	end
end

return

% extracts extra variables that only exist for stattype Z
function data = get_depth_variables(data)

dz=diff(data.z(1:2,1));
data.FlowNu = sum(data.Nu,1)*dz/diff(data.Domain([5 6]));
data.FlowDensity = sum(data.Density,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumX = sum(data.MomentumX,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumY = sum(data.MomentumY,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentumZ = sum(data.MomentumZ,1)*dz/diff(data.Domain([5 6]));
data.FlowVelocityX = data.FlowMomentumX./data.FlowDensity;
data.FlowVelocityY = data.FlowMomentumY./data.FlowDensity;
data.FlowVelocityZ = data.FlowMomentumZ./data.FlowDensity;

data.FlowMomentum(1,:) = sum(data.MomentumX,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentum(2,:) = sum(data.MomentumY,1)*dz/diff(data.Domain([5 6]));
data.FlowMomentum(3,:) = sum(data.MomentumZ,1)*dz/diff(data.Domain([5 6]));
data.FlowVelocity = data.FlowMomentum./data.FlowDensity([1,1,1],:);
FlowMomentumFlux(1,:) = sum(data.MomentumFluxXX,1)*dz/diff(data.Domain([5 6]));
FlowMomentumFlux(2,:) = sum(data.MomentumFluxYY,1)*dz/diff(data.Domain([5 6]));
FlowMomentumFlux(3,:) = sum(data.MomentumFluxZZ,1)*dz/diff(data.Domain([5 6]));
data.FlowVelocity2 = FlowMomentumFlux./data.FlowDensity([1,1,1],:);


data.BottomFrictionCoefficient = -data.StressXZ(1,:)/data.StressZZ(1,:);

% height and average density of the flow
if ~(data.nx==1&&data.ny==1), 
  ContactStressZZ=data.NormalStressZZ+data.TangentialStressZZ;
  data.FlowHeight=max(data.z.*(ContactStressZZ~=0),[],1);
else
  kappa=0.02;
  IntDensity = cumsum(data.Density);
  Base = min(data.z(IntDensity>kappa*max(IntDensity)));
  Surface = max(data.z(IntDensity<(1-kappa)*max(IntDensity)));
  data.FlowHeight=(Surface-Base)/(1-2*kappa);
  data.Base=Base-kappa*data.FlowHeight;
  data.Surface=Surface+kappa*data.FlowHeight;

  FlowHeightIndex = data.z>data.Base & data.z<data.Base+data.FlowHeight;
  if ~isempty(FlowHeightIndex)
    data.Froude=norm([data.FlowVelocityX data.FlowVelocityY data.FlowVelocityZ])/sqrt(data.FlowHeight*(-data.Gravity(3)));
  end
end

% % height and average density of the flow
% FlowHeightIndex = min([length(data.Nu),round(length(data.Nu)/2)-1+find(data.Nu(round(end/2):end)<0.3,1)]);
% if ~isempty(FlowHeightIndex)
% 	AveragingInterval = ceil(FlowHeightIndex/3):ceil(FlowHeightIndex*2/3);
% 	data.FlowDensity = sum(data.Density(AveragingInterval))/length(AveragingInterval);
% 	data.FlowNu = sum(data.Nu(AveragingInterval))/length(AveragingInterval);
% 	data.FlowHeight = data.z(1) + sum(data.Density)/length(data.z)/data.FlowDensity * (data.z(end)-data.z(1));
% 
%   index = data.z>0&data.z<data.FlowHeight;
% 	data.FlowVelocityX = average(data.VelocityX(index));
% 	data.FlowVelocityY = average(data.VelocityY(index));
% 	data.FlowVelocityZ = average(data.VelocityZ(index));
%   data.Froude=norm([data.FlowVelocityX data.FlowVelocityY data.FlowVelocityZ])/sqrt(data.FlowHeight*(-data.Gravity(3)));
% end


return

function [Vs, Ds] = dsort_ev(V, D)
    Ds=D;
    Vs=V;
%     tmpvec=[0;0;0];
    if((Ds(1,1))<(Ds(2,2)))
        tmp=Ds(1,1);
        Ds(1,1)=Ds(2,2);
        Ds(2,2)=tmp;
        tmpvec=Vs(:,1);
        Vs(:,1)=Vs(:,2);
        Vs(:,2)=tmpvec;
    end
    if((Ds(2,2))<(Ds(3,3)))
        tmp=Ds(2,2);
        Ds(2,2)=Ds(3,3);
        Ds(3,3)=tmp;
        tmpvec=Vs(:,2);
        Vs(:,2)=Vs(:,3);
        Vs(:,3)=tmpvec;
    end
    if((Ds(1,1))<(Ds(2,2)))
        tmp=Ds(1,1);
        Ds(1,1)=Ds(2,2);
        Ds(2,2)=tmp;
        tmpvec=Vs(:,1);
        Vs(:,1)=Vs(:,2);
        Vs(:,2)=tmpvec;
    end
    if((Ds(2,2))<(Ds(3,3)))
        tmp=Ds(2,2);
        Ds(2,2)=Ds(3,3);
        Ds(3,3)=tmp;
        tmpvec=Vs(:,2);
        Vs(:,2)=Vs(:,3);
        Vs(:,3)=tmpvec;
    end

return

function [remainder,maximum] = get_momentum_equation(data)
    
NablaMomentumX = nabla(data.ParticleDensity*data.Nu.*data.VelocityX,data.x,data.y,data.z);
NablaMomentumY = nabla(data.ParticleDensity*data.Nu.*data.VelocityY,data.x,data.y,data.z);
NablaMomentumZ = nabla(data.ParticleDensity*data.Nu.*data.VelocityZ,data.x,data.y,data.z);

VelocityDotNablaMomentum = [...
  sum([data.VelocityX;data.VelocityX;data.VelocityZ].*NablaMomentumX,1); ...
  sum([data.VelocityX;data.VelocityX;data.VelocityZ].*NablaMomentumY,1); ...
  sum([data.VelocityX;data.VelocityX;data.VelocityZ].*NablaMomentumZ,1) ];

NablaDotStress = nabla([data.StressXX;data.StressXY;data.StressXZ;data.StressYY;data.StressYZ;data.StressZZ],data.x,data.y,data.z);

DensityGravity = data.ParticleDensity * data.Gravity' * data.Nu;

remainder = VelocityDotNablaMomentum -DensityGravity +NablaDotStress;

maximum = max([max(VelocityDotNablaMomentum) max(DensityGravity) max(NablaDotStress)]);

function data = read_ene(statname,data)

%load ene data
filename = [statname(1:end-5) '.ene'];
if ~exist(filename,'file'), 
  %if unsuccessful, load ene data with last part of filename removed
  dots=strfind(statname,'.');
  filename = [statname(1:min(dots(end-1:end))) 'ene'];
  if ~exist(filename,'file')
    disp([filename ' not found']); 
  else
    %if unsuccessful, give up
    disp([statname(1:end-5) '.ene not found, using ' filename ' instead']); 
  end
end

if exist(filename,'file'),
  % load raw data from file, with the first three lines as header files
  % (such that matlab recognizes the array size)
  rawdata = importdata(filename,' ',1);
  % if tabs are used as delimiters
  if (size(rawdata.data,2)~=8)
    rawdata = importdata(filename,'\t',1);
  end

  data.Ene.Time = rawdata.data(:,1);
  data.Ene.Gra = rawdata.data(:,2);
  data.Ene.Kin = rawdata.data(:,3);
  data.Ene.Rot = rawdata.data(:,4);
  data.Ene.Ela = rawdata.data(:,5);
  data.Ene.ComX = rawdata.data(:,6);
  data.Ene.ComY = rawdata.data(:,7);
  data.Ene.ComZ = rawdata.data(:,8);
else
  disp([filename ' not found'])
end

return

function data = read_restart(statname,data)

%load restart data
filename = [statname(1:end-5) '.restart'];
fid=fopen(filename);
if (fid==-1), 
  %if unsuccessful, load restart data with last part of filename removed
  dots=strfind(statname,'.');
  filename = [statname(1:min(dots(end-1:end))) 'restart'];
  fid=fopen(filename);
  if (fid==-1),
    disp([filename ' not found']); 
  else
    %if unsuccessful, give up
    disp([statname(1:end-5) '.restart not found, using ' filename ' instead']); 
  end
end
tdata = textscan(fid,'%s');
tdata = tdata{1};
fclose(fid);

%this assumes monodispersed particles
if (length(tdata{1})==1) 
	%old restart files
	data.ParticleDensity = str2double(tdata(29));
	data.d = 2*str2double(tdata(28));
	data.Gravity = str2double(tdata(2:4))';
	data.N = str2double(tdata(end-3));
	data.Domain=str2double(tdata(5:10))';
else
	%new restart version
	i = find(strcmp(tdata,'rho'));
  if length(i)>1, disp('warning: two species'); end
	data.ParticleDensity = str2double(tdata(i+1));
	i = find(strcmp(tdata,'mu'));
	data.Mu = str2double(tdata(i+1));
	i = find(strcmp(tdata,'gravity'));
	data.Gravity = str2double(tdata(i+(1:3)))';
	i = find(strcmp(tdata,'Particles'));
	data.N = str2double(tdata{i+1});
% 	data.Radii = str2double(tdata(i+(8:15:15*data.N)))';
% 	data.Speed = sqrt( ... 
%       str2double(tdata(i+(5:15:15*data.N)))'.^2 ...
%     + str2double(tdata(i+(6:15:15*data.N)))'.^2 ...
%     + str2double(tdata(i+(7:15:15*data.N)))'.^2);
% 	data.Depth = str2double(tdata(i+(4:15:15*data.N)))';
	i = find(strcmp(tdata,'xmin'));
	data.Domain = str2double(tdata(i+(1:2:11)))';
	i = find(strcmp(tdata,'FixedParticleRadius'));
	if ~isempty(i), data.FixedParticleRadius = str2double(tdata{i+1}); end
	i = find(strcmp(tdata,'MinInflowParticleRadius'));
	if ~isempty(i), 
		data.InflowParticleRadius = str2double(tdata(i+[1 3]))'; 
		data.d = sum(data.InflowParticleRadius);
	else
		i = find(strcmp(tdata,'Particles'));
		data.d = 2*str2double(tdata{i+8});
	end
end
data.ParticleVolume = pi/6*data.d^3;
data.DomainVolume = prod(data.Domain([2 4 6])-data.Domain([1 3 5]));
data.ChuteAngle = round(atand(-data.Gravity(1)/data.Gravity(3))*400)/400;
return

function NablaVariable = nabla(variable,x,y,z)

% stores the grid dimensions
n(1) = length(x)/sum(x(1)==x);
n(2) = length(y)/sum(y(1)==y);
n(3) = length(z)/sum(z(1)==z);

% stores coordinates in grid
X = zeros(n(end:-1:1));
X(:) = x;
Y = zeros(n(end:-1:1));
Y(:) = y;
Z = zeros(n(end:-1:1));
Z(:) = z;
V = zeros(n(end:-1:1));

% stores differentials (small number if X/Y/Z is constant)
dX = max(1e-60,(X(:,:,[2:end end])-X(:,:,[1 1:end-1])));
dY = max(1e-60,(Y(:,[2:end end],:)-Y(:,[1 1:end-1],:)));
dZ = max(1e-60,(Z([2:end end],:,:)-Z([1 1:end-1],:,:)));

if size(variable,1)==1 %scalar variable
  V(:) = variable;
  dXV = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  dYV = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  dZV = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;
  
  NablaVariable = [dXV(:) dYV(:) dZV(:)]';
elseif size(variable,1)==3 %vector
  V(:) = variable(1,:);
  dXV = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  V(:) = variable(2,:);
  dYV = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  V(:) = variable(3,:);
  dZV = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;

  NablaVariable = dXV(:) + dYV(:) + dZV(:);
elseif size(variable,1)==6 %symmetric matrix
  V(:) = variable(1,:);
  dXVXX = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  V(:) = variable(2,:);
  dXVYX = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  V(:) = variable(3,:);
  dXVZX = (V(:,:,[2:end end])-V(:,:,[1 1:end-1])) ./dX;
  
  V(:) = variable(2,:);
  dYVXY = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  V(:) = variable(4,:);
  dYVYY = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  V(:) = variable(5,:);
  dYVZY = (V(:,[2:end end],:)-V(:,[1 1:end-1],:)) ./dY;
  
  V(:) = variable(3,:);
  dZVXZ = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;
  V(:) = variable(5,:);
  dZVYZ = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;
  V(:) = variable(6,:);
  dZVZZ = (V([2:end end],:,:)-V([1 1:end-1],:,:)) ./dZ;

  NablaVariable = [dXVXX(:)+dYVXY(:)+dZVXZ(:) dXVYX(:)+dYVYY(:)+dZVYZ(:) dXVZX(:)+dYVZY(:)+dZVZZ(:)]';
else 
  disp('Error');
end
return
