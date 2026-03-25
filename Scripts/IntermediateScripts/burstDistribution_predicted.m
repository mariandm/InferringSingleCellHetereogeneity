% Code by Marian Dominguez-Mirazo, 2025
% Caution: this script may take >5 min to run
% Calculate the predicted burst size distribution based on the predicted
% latent period distribution, and the predicted latent period to burst size
% relationship
% It requires the pre-computed files: 
% - ~/IntermediateFiles/LS_burst.mat
% - ~/IntermediateFiles/MLE_CDF.mat
% Output: ~/IntermediateFiles/burstDistribution_predicted.mat

%%
clear all; close all; clc;
addpath('../../Data');
addpath('../Functions/');
addpath('../../IntermediateFiles/');
%% Load best LP to burst size linear model
load('LS_burst.mat');

%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));
[thisshape,thisscale] = gamma_fromstats(CI(1,1),CI(2,1));

%% Get burst data
file = 'SingleCellData.csv';
tab = readtable(file, 'ReadVariableNames', true);
tab_tmp = table2array(tab(:,3:end));
tab_tmp(tab_tmp < 2) = NaN;

% Sampling time points are in the header
tab = table2array(readtable(file, 'ReadVariableNames', false));
sample_time = tab(1,3:end);

%% Calculate expected burst size distribution

% Starting at 9hr after initial phage-bacteria incubation
idx = sample_time >= 9;
times = sample_time(idx);

% Initialize storage
bursts = [];
bursts_raw = [];

for i = 1:numel(times)
    % Randomly sample the LP dist
    random_LP = gamrnd(thisshape,thisscale,1e8,1);
    % Get only the bursts we would be able to observe by that sampling
    % timepoint
    random_LP_corrected = random_LP(random_LP <= times(i));
    % Convert LP to bursts according to the predicted relationship
    prebursts = burst_linear(real_best_linear_expdecay,random_LP_corrected);
    prebursts_raw = burst_linear(real_best_linear_expdecay,random_LP);
    % Account for decay by sticking to well plastic
    bursts = [bursts;exp(-0.04 .* (times(i) - random_LP_corrected)) .* prebursts];
    bursts_raw = [bursts_raw,prebursts_raw];
end

%% Calculate distributions
% Get histogram for the experimental data
thisdata = tab_tmp(:,idx); 
dh_y = histcounts(thisdata,'BinWidth',10);

% Get expected histogram with both corrections (sampling and adhesion)
mh = histogram(bursts,'BinWidth',1,BinLimits=[0,250]);
mh_y = mh.Values;
mh_y_new = mh_y / sum(mh_y) * sum(dh_y)*10;
mh_minpoints = (mh.BinEdges(1:end-1)+mh.BinEdges(2:end))/2;

% Get expected histogram with no corrections
mh_raw = histogram(bursts_raw,'BinWidth',1,BinLimits=[0,250]);
mh_y_raw = mh_raw.Values;
mh_y_new_raw = mh_y_raw / sum(mh_y_raw) * sum(dh_y)*10;
mh_minpoints_raw = (mh_raw.BinEdges(1:end-1)+mh_raw.BinEdges(2:end))/2;

%%
save('../../IntermediateFiles/burstDistribution_predicted.mat',...
    'mh_minpoints_raw','mh_y_new_raw',...
    'mh_minpoints','mh_y_new',...
    'thisdata','dh_y');

