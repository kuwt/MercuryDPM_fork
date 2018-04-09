function [data, info] = load_data_file()
%load stat file data
filename='tav8.3400.data';
disp(['loading data from ' filename]);

% load raw data from file, with the first three lines as header files (such
% that matlab recognizes the array size)
rawdata = importdata(filename,' ',17);
  
% reads variable names from line 1 and text output
% from line 2
info.variable_names = {'PVF'; 'Velocity_x'; 'Velocity_y'; 'Velocity_z'; 'StressXX'; 'StressXY'; 'StressXZ'; 'StressYY'; 'StressYZ'; 'StressZZ'; 'Temperature'; 'Pressure'};
info.text = ' ';

% writes the raw data into bits of data for each timestep
index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
datafull.time = 0;
datafull.coordinates = rawdata.data(:,4:6);
datafull.epsilon = rawdata.data(:,(61:69));
datafull.sigma = (rawdata.data(:,[21:29]) + rawdata.data(:,[31:39]));

% stores the grid dimensions
for i=1:3
  n(i) = size(datafull.coordinates,1) / sum(datafull.coordinates(1,i)==datafull.coordinates(:,i));
end

% stores coordinates in grid
for i=1:3
  if (length(n)==1)
    X{i} = datafull.coordinates(:,i);
  else
    X{i} = zeros(n(end:-1:1));
    X{i}(:) = datafull.coordinates(:,i);
  end
end
Z = zeros(size(X{1}));

dim = boolean([0 0 1]);
dimavg = find(~dim);

data.z = squeeze(sum(sum(X{dim},4-dimavg(2)),4-dimavg(1))/prod(n(~dim)));
for i=1:9
  Z(:) = datafull.epsilon(:,i);
  data.epsilon(:,i) = squeeze(sum(sum(Z,4-dimavg(2)),4-dimavg(1))/prod(n(~dim)));
  Z(:) = datafull.sigma(:,i);
  data.sigma(:,i) = squeeze(sum(sum(Z,4-dimavg(2)),4-dimavg(1))/prod(n(~dim)));
end

plot_epsilon(data);
plot_sigma(data);
plot_dev(data);
return

function plot_epsilon(data)

figure(1);
domain = [min(data.z) max(data.z) min(data.epsilon(:))/9 max(data.epsilon(:))];
for i=1:9
  subplot(3,3,i); 
  plot(data.z,data.epsilon(:,i));
  xlabel('z');
  title(['\epsilon_' num2str(i)]);
  axis(domain);
end

return

function plot_sigma(data)

figure(2);
domain = [min(data.z) max(data.z) min(data.sigma(:)) max(data.sigma(:))];
for i=1:9
  subplot(3,3,i); 
  plot(data.z,data.sigma(:,i));
  xlabel('z');
  title(['\sigma_' num2str(i)]);
  axis(domain);
end

return

function plot_dev(data)

for i=1:length(data.z)
  epsilon = zeros(3);
  epsilon(:) = data.epsilon(i,:);
  epsilon = .5*(epsilon + epsilon');

  sigma = zeros(3);
  sigma(:) = data.sigma(i,:);
  sigma = .5*(sigma + sigma');
  sigma = sigma - eye(3) * trace(sigma)/3;
  
  %eigendecomposition of sigma and epsilon (ve*de*ve'=epsilon)
  [vs,ds] = eig(sigma);
  [ve,de] = eig(epsilon);

  %rotate sigma and epsilon into the epsilon reference frame
  epsilon_eps = ve'*epsilon*ve;
  sigma_eps = ve'*sigma*ve;

  %calculate the deviation
  diage = diag(epsilon_eps);
  diags = diag(sigma_eps);
  E = dot(diage,diags) / dot(diage,diage);
  
  sigma_d_eps = sigma_eps - E * epsilon_eps;
  sigma_d = ve*sigma_d_eps*ve';

  %calculate the frobenius norm
  normStokes = norm(E*epsilon,'fro');%sqrt(sum(sum((E*epsilon).^2)));
  normSigma = norm(sigma_eps,'fro');%sqrt(sum(sum(sigma_eps.^2)));
  normNonStokes = norm(sigma_d,'fro');%sqrt(sum(sum(sigma_d.^2)));
  data.nonStokeness(i) = normNonStokes/normSigma;
  data.E(i) = E;
  data.normEpsilon(i) = norm(epsilon,'fro');%sqrt(sum(sum(epsilon.^2)));
  data.normSigma(i) = normSigma;
end

figure(3);
subplot(2,2,1);
plot(data.z,data.nonStokeness);
xlabel('z');
title(['non-Stokeness min_E |\sigma - E \epsilon|_{Frob}/|\sigma|_{Frob}']);
subplot(2,2,2);
plot(data.z,data.E);
xlabel('z');
title(['E']);
subplot(2,2,3);
plot(data.z,data.normEpsilon);
xlabel('z');
title(['|\epsilon|_{Frob}']);
v = axis();
axis([v(1:3) v(4)*.1]);
subplot(2,2,4);
plot(data.z,data.normSigma);
xlabel('z');
title(['|\sigma|_{Frob}']);

saveas(3 ,['figure' num2str(3 ) '.pdf']);


return

% #
% # columns 1-3
%              integer cell positions (iix,iiy,iiz) - starting from 0
% #
% # columns 4-6
%              real cell positions (xpos,ypos,zpos) - lower,left corner
% #
% # column 7
%              averaging volume of cell (Vxyz)-V_{cell}
% #
% # columns 8,9
%              total number of particles and contacts
% #
% # columns 10-12
%              volume fraction xnu=<nu> and xnustdev=sqrt(<nu^2>-<nu>^2), 
%                                       both based on mass-center average
%              and volume fraction xCnu_av=<nu>_{av} based on smoothing function
% #
% # columns 13-18
%              cc0 ... number of particles with      C == 0 contacts
%              cc3 ... number of particles with 1 <= C <= 3 contacts
%              cc4 ... number of particles with      C == 4 contacts
%              cc6 ... number of particles with 5 <= C <= 6 contacts
%              cc9 ... number of particles with 7 <= C <= 9 contacts
%              ccx ... number of particles with 10 <= C     contacts
% #
% # columns 19,20
%              static pressure p1=<p> and standard deviation pstdev=sqrt(<p^2>-<p>^2)
% # columns 21-29
%              normal static stress tensor sigma_{ab}=(1/V) \sum_c l_a f^n_b
%              indices are (1-9): 11, 12, 13, 21, 22, 23, 31, 32, 33
% #
% # column 30
%              tangential static pressure pt=<p^t> 
% # columns 31-39
%              tangential static stress tensor sigma^t_{ab}=(1/V) \sum_c l_a f^t_b
%              indices are (1-9): 11, 12, 13, 21, 22, 23, 31, 32, 33
% #
% # column 40
%              dynamic pressure pt=<p^d> 
% # columns 41-49
%              dynamic stress tensor sigma^t_{ab}=(1/V) \sum_p m_p v^p_a v^p_b
%              indices are (1-9): 11, 12, 13, 21, 22, 23, 31, 32, 33
% #
% # column 50
%              time t
% #
% # columns 51-53
%              average velocities <v_1>, <v_2>, <v_3>
%              for the smoothed velocity field, see col. 114-116
% #
% # columns 54-56
%              average squared-velocities <v_1^2>, <v_2^2>, <v_3^2>
% #
% # columns 57-59
%              average displacements per time-step <d_1>/dt, <d_2>/dt, <d_3>/dt
% #
% # column 60
%              isotropic velocity gradient 
% # columns 61-69
%              velocity gradient dva/db=(va(b+db)-va(b-db)/2db
%              (based on the displacement field - NOT on velocities)
%              indices are (1-9): 11, 12, 13, 21, 22, 23, 31, 32, 33
% #
% # column 70
%              isotropic fabric (contact number density)
% # columns 71-79
%              fabric F_ab = (1/V) \sum_pc V^pc n^c_a n^c_b
%              indices are (1-9): 11, 12, 13, 21, 22, 23, 31, 32, 33
% #
% # column 80
%              normal isotropic rel. strain
% # columns 81-89
%              normal rel. strain eps_n = \sum_pc (delta/l) n^c_a n^c_b
% #
% # column 90
%              tangential isotropic rel. strain
% # columns 91-99
%              tangential rel. strain eps_n = \sum_pc (delta/l) n^c_a n^c_b
% #
% # columns 100,101
%              smoothed/weighted average of particle- and contact-number
% #
% # columns 102-106
%              First five moments (smoothed) of the size distribution
%              <a>, <a^2>, <a^3>, <a^4>, <a^5>
% #
% # columns 107-109
%              systems size in 1,2,3-directions (see col. 117-119 for ref0)
% #
% # column 110
%              not-used
% #
% # columns 111-113
%              average spin <w_a>
% #
% # columns 114-116
%              smoothed velocity field <v_a>_{av} <-> col. 51-53
% #
% # columns 117-119
%              reference systems size in 1,2,3-directions 
%              (can be used to compute dimensionless change/strain
%               see col. 107-109 for actual size)

% # columns 201-215 (same as the CN-stat-line)
% # 201        <nu>
% # 202        N_total
% # 203        N_0
% # 204-223    N_{index=column-203}
% # 224        C_total
% # 225        C_0 (always=0)
% # 226-245    C_{index=column-225}
% #