% Code by Marian Dominguez-Mirazo, 2026
% Calculate the predicted burst size distribution based on the predicted
% latent period distribution, and the predicted latent period to burst size
% relationship, assuming eclipse period variability
% It requires the pre-computed files: 
% - ~/IntermediateFiles/LS_burst.mat
% - ~/IntermediateFiles/MLE_CDF.mat
% Output: ~/IntermediateFiles/burstDistribution_predicted_eclipse.mat

%%
clear all; close all; clc;
addpath('../../Data');
addpath('../Functions/');
addpath('../../IntermediateFiles/');

%%
v_vals = [0.05, 0.25, 0.5, 0.75];

%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));

%% Load best LP to burst size linear model
load('LS_Burst_eclipse.mat');

%% Get burst data
file = 'SingleCellData.csv';
tab = readtable(file, 'ReadVariableNames', true);
tab_tmp = table2array(tab(:,3:end));
tab_tmp(tab_tmp < 2) = NaN;

% Sampling time points are in the header
tab = table2array(readtable(file, 'ReadVariableNames', false));
sample_time = tab(1,3:end);

%%
% Get histogram for the experimental data
idx = sample_time >= 9;
thisdata = tab_tmp(:,idx); 
dh_y = histcounts(thisdata,'BinWidth',10);

% initialize storage
mh_y_new_raw_all = zeros(numel(v_vals),250);
mh_minpoints_raw_all = zeros(numel(v_vals),250);

for i = 1:numel(v_vals)
    
    var_percentage = v_vals(i);
    idx_bl = find(bestlinear_allv(:,4)==var_percentage);
    
    % Calculate eclipse and post-eclipse distributions
    emean = bestlinear_allv(idx_bl,1);
    [CIdist_e,CIdist_pe] = splitDistr_eclipse(CI,var_percentage,emean);
    %
    [thisshape_e,thisscale_e] = gamma_fromstats(CIdist_e(1,1),CIdist_e(2,1));
    [thisshape_pe,thisscale_pe] = gamma_fromstats(CIdist_pe(1,1),CIdist_pe(2,1));

    % Get random post-eclipse values
    random_pe = gamrnd(thisshape_pe,thisscale_pe,1e8,1);
    % Calculate burst
    bursts_raw = burst_linear([0,bestlinear_allv(idx_bl,2)],random_pe);

    % Get expected histogram with no corrections
    mh_raw = histogram(bursts_raw,'BinWidth',1,BinLimits=[0,250]);
    mh_y_raw = mh_raw.Values;
    mh_y_new_raw = mh_y_raw / sum(mh_y_raw) * sum(dh_y)*10;
    mh_minpoints_raw = (mh_raw.BinEdges(1:end-1)+mh_raw.BinEdges(2:end))/2;

    % store 
    mh_y_new_raw_all(i,:) = mh_y_new_raw;
    mh_minpoints_raw_all(i,:) = mh_minpoints_raw;
end
%%
save('../../IntermediateFiles/burstDistribution_predicted_eclipse.mat',...
    'mh_y_new_raw_all',...
    'mh_minpoints_raw_all',...
    'thisdata','dh_y');