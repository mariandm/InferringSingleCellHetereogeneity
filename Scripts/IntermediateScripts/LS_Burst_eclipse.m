% Code by Marian Dominguez-Mirazo, 2026
% Caution: This script may take several ours to run, the resulting file has
% been precomputed for convenience
% Based on RMSE find the parameter set that nest describes the cumulative
% burst size of the single cell protocol assuming variability in the latent
% period is present and accounts of a proportion of total latent period
% variability
% It requires the pre-computed files: 
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/IntermediateFiles/burstMean_perReplicate.csv
% Output: ~/IntermediateFiles/LS_Burst_eclipse.mat
%%
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/')

vvals = [0.0001,0.05:0.05:0.95];
%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));
[thisshape,thisscale] = gamma_fromstats(CI(1,1),CI(2,1));

%% Get effective burst size per replicate
% Load file
file = 'burstMean_perReplicate.csv';
tab = table2array(readtable(file, 'ReadVariableNames', false));
sample_time = tab(1,:);
efburst_per_replicate_og = tab(2:end,:);

%% Find best fit for the models
nrand = 100; % Number of initial search points
best_linear = zeros(numel(vvals),4);
best_mm = zeros(numel(vvals),5);
best_kannoly = zeros(numel(vvals),6);

for i = 1:numel(vvals)
    best_linear(i,1:3) = burst_bestfit_eclipse(efburst_per_replicate_og,sample_time,CI,vvals(i),'linear',nrand);
    best_linear(i,4) = vvals(i);
    best_mm(i,1:4) = burst_bestfit_eclipse(efburst_per_replicate_og,sample_time,CI,vvals(i),'mm',nrand);
    best_mm(i,5) = vvals(i);
    best_kannoly(i,1:5) = burst_bestfit_eclipse(efburst_per_replicate_og,sample_time,CI,vvals(i),'kannoly',nrand);
    best_kannoly(i,6) = vvals(i);
end

%% Select best parameter search for each model and store
[A,I] = min(best_linear(:,3));
bestlinear = best_linear(I,:);
[A,I] = min(best_mm(:,4));
bestmm = best_mm(I,:);
[A,I] = min(best_kannoly(:,5));
bestkannoly = best_kannoly(I,:);

bestlinear_allv = best_linear;

save('../../IntermediateFiles/LS_Burst_eclipse.mat','bestlinear','bestmm','bestkannoly','bestlinear_allv');