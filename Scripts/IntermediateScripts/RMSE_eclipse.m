% Code by Marian Dominguez-Mirazo, 2026
% This script calculates the RMSE for a range of parameters assuming a
% linear model of production period to burst size relationship, and
% assuming variability in the latent period. It stores values needed to
% produce SupplFig11.
% It requires the pre-computed files: 
% - ../../IntermediateFiles/MLE_CDF.mat
% - ../../IntermediateFiles/burstMean_perReplicate.csv
% Output: '../../IntermediateFiles/RMSE_eclipse.mat'
%% This scripts makes a figure to check sloppines of p parameter
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/');
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
%%
% Loop multiple values of varp (percentage of LP var explained by eclipse)
range_varp = [0.0001,0.05:0.05:0.95];
% And r, post-eclipse period to burst size slope
range_r = 10:2.5:70;
% And average eclipse period, 3 examples
range_d = [4.5,5.5,6.5];
% Initialize storage
RMSE_matrix = zeros(numel(range_d),numel(range_varp),numel(range_r));
model="linear";
% Average eclipse time
for k = 1:numel(range_d)
    d = range_d(k);
    for i = 1:numel(range_varp)
        varp = range_varp(i);
        for j = 1:numel(range_r)
            r = range_r(j);
            pars = [d,r];
            RMSE_matrix(k,i,j) = burst_RSME_eclipse(efburst_per_replicate_og,sample_time,CI,varp,pars,model);
        end
    end
end

%% Save the thing
save('../../IntermediateFiles/RMSE_eclipse.mat','RMSE_matrix','range_d','range_r','range_varp');