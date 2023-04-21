%% CompareInsertedPSD
% This script can check if the inserted by an InsertionBoundary of
% MercuryDPM fits to the AnalyticalPSD one read in through a .csv file. The
% .data file generates all particle information needed for the
% simulationPSD.
% 
% Input:  SimulationPSD.csv  A .csv file with particle information copied
%                            from the .data file of MercuryDPM
%         AnalyticalPSD.csv  A .csv file created from the input 2D array
%                            used to feed the InsertionBoundary.
%                            (see PSD class)
%         
% Output: cumulative probability density function as histogram figure of
%         particles inserted into the simulation against a line plot
%         representing the analyticalPSD
%%
clc
close all
clear all

% Read distribution from particle simulations
T = readtable('simulationPSD.csv');
simulationPSD = table2array(T);
% set histogram variables
radius = simulationPSD(:,7);
numberOfBins = 3000;
% plot histogram of simulation radii
figure(3)
histogram(radius,numberOfBins,"Normalization", "cdf")
% xlim([0, 3e-3]);
hold on

% Read distribution from input to particle simulations
AnalyticalPSD = readmatrix('AnalyticalPSD.csv', 'delimiter', ',');
% validate cumulative probability density function
AnalyticalPSD = validateCDF(AnalyticalPSD);
% convert if necessary
% T2 = convertCDFtoPDF(T2);
% T2 = convertPVDFtoPNDF(T2);
% T2 = convertPDFtoCDF(T2);

% Conversion from CDF to PDF
% for j = 1:length(T2)-1
%     T2(j,2) = T2(j+1,2)-T2(j,2);
% end


% plot analyitical solution (scale the radii if necessary)
AnalyticalRadius = AnalyticalPSD(:,1);
AnalyticalProbability = AnalyticalPSD(:,2);
plot(AnalyticalRadius,AnalyticalProbability, "linewidth", 2);
legend("PSD MDPM", "Analytical PSD", "location", "northwest")
xlabel("Size [m]")
ylabel("Probability [%]")