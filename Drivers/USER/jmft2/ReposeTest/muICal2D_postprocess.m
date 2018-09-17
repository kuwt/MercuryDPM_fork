simdir = '~/MercuryDPM/MercuryBuild/Drivers/USER';
drivername = 'ReposeTest';

% seriesname = 'blasius-grit';
% speciesname = 'smoothbit';
% angles = [31:40];

seriesname = 'mixedspecies';
speciesname = 'ballotini-on-sand';
angles = [15:40];

seriesname = 'blasius-cali';
speciesname = 'smoothbit';
angles = [10:25];
% speciesname = 'roughbit';
% angles = [15:40];

domainlength = 1;
diam = 0.008;

%%
xmoms = zeros(size(angles));
depths = zeros(size(angles));
masses = zeros(size(angles));
phis = zeros(size(angles));
nframes = 1;
for i = 1:length(angles)
   probname = sprintf('%s-%.1f', ...
       speciesname, angles(i));
   muIfn = sprintf('%s/%s/%s/%s.muI', ...
       simdir, drivername, seriesname, probname);
   muIfile = importdata(muIfn);
   xmoms(i)  = mean(muIfile.data(end-nframes+1:end, 5));
   depths(i) = mean(muIfile.data(end-nframes+1:end, 3));
   masses(i) = mean(muIfile.data(end-nframes+1:end, 4));
   fprintf('For %s, t=%f\n', probname, muIfile.data(end,1));
end
phis = masses ./ (domainlength * depths);
clear probname enefn enefile

%% Try to fit the function mu(I). 
% We know the depth of the current and the flow rate,
% which are meant to be related by the formula
%   Q = (2 I sqrt(g cos(theta))) / (5 d) * h^(5/2)
% and so we can obtain the inertial number
%   I = (5dQ/2) / (h^(5/2) cos(theta)^(1/2))

% inertials = (5*diam/2) * ( xmoms ./ (depths.^(5/2) .* cos(angles * pi/180).^(1/2)) );
inertials = (5*diam)/(2*domainlength) * (xmoms) ./ (depths.^5 .* cos(angles * pi/180) .* phis.^3).^(1/2);
muIs = tan(angles * pi/180);

muI_fittype = fittype('a + (b-a)/(c/x + 1)');
muI_fit = fit(inertials', muIs', muI_fittype, ...
    'StartPoint', [min(muIs), max(muIs), 1], ...
    'Lower', [0 0 0] ...
)
mu1 = muI_fit.a;
mu2 = muI_fit.b;
I0 = muI_fit.c;
muI_fit = @(I) (mu1 + (mu2 - mu1)./(I0 ./ I + 1));

inertials_interpolated = linspace(min(0, min(inertials)), max(inertials));
muIs_interpolated = arrayfun(muI_fit, inertials_interpolated);

muICal2D_results = struct(...
    'seriesname', seriesname, 'speciesname', speciesname, ...
    'domainlength', domainlength, 'diam', diam, ...
    'xmoms', xmoms, 'depths', depths, 'masses', masses, ...
    'inertials', inertials, 'muIs', muIs, ...
    'mu1', mu1, 'mu2', mu2, 'I0', I0, 'muI_fit', muI_fit, ...
    'inertials_interpolated', inertials_interpolated, ...
    'muIs_interpolated', muIs_interpolated ...
);

%%
figure;

subplot(2,3,1);
plot(angles, xmoms, 'kx');
title(sprintf('%s : %s : %s', drivername, seriesname, speciesname));
xlabel('angle');
ylabel('x-momentum (flow rate)');

subplot(2,3,2);
plot(angles, depths, 'kx');
% title(sprintf('%s : %s : %s', drivername, seriesname, speciesname));
title('depth against angle');
xlabel('angle');
ylabel('depth (2 y-com)');

subplot(2,3,3);
% plot(angles, inertials, 'kx');
% xlabel('angle');
% ylabel('inertial number I');
plot(angles, phis, 'kx');
title('volume fraction against angle');
xlabel('angle');
ylabel('\phi');

%%
subplot(2,3,4);
plot(inertials, muIs, 'kx', ...
     inertials_interpolated, ...
     muIs_interpolated, 'k-'  ...
);
title(sprintf('\\theta_1 = %.2f, \\theta_2 = %.2f, I_0 = %f', ...
    atan(mu1) * 180/pi, atan(mu2) * 180/pi, I0));
xlabel('inertial number I');
xlim([0,inf]);
ylabel('$\mu(I) = \tan\theta$', 'Interpreter', 'LaTeX');
legend('data', 'Pouliquen fit', 'Location', 'southeast');
%%
subplot(2,3,5);
mid_inertials = (inertials(1:end-1) + inertials(2:end))/2;
d_muIs        = (muIs(2:end) - muIs(1:end-1)) ./ (inertials(2:end) - inertials(1:end-1));

fitting = polyfit(...
    log(mid_inertials(d_muIs > 0)), ...
    log(d_muIs(d_muIs > 0)), ...
1);

loglog(mid_inertials, d_muIs, 'kx', ...
       mid_inertials, exp( fitting(1) * log(mid_inertials) + fitting(2) ), 'k-' ...
);
title(sprintf('$\\mathrm{d}\\mu/\\mathrm{d}I \\sim I^{-\\alpha},\\, \\alpha = %f$', -fitting(1)), ...
     'Interpreter', 'LaTeX');
xlabel('inertial number I');
xlim([I0/100, I0*10]);
ylabel('$\mathrm{d}\mu/\mathrm{d}I$', 'Interpreter', 'LaTeX');
ylim([-inf, 100]);
legend('data', 'power law fit');

%%
subplot(2,3,6);
plot(inertials, phis, 'kx');
title('$\phi$ against $I$', 'Interpreter', 'LaTeX');
xlabel('inertial number I');
xlim([0,inf]);
ylabel('\phi');