% Code by Marian Dominguez-Mirazo, 2026
% This script takes pre-computed latent period distribution with confidence
% intervals, pre-computed latent period to burst size linear model with 
% confidence intervals and computes the mean and CV of the burst size 
% distribution as predicted from the latent period to burst size
% relationship and latent period distribution
% Required pre-computed files:
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/IntermediateFiles/LS_burst.mat
% - ~/IntermediateFiles/bootstrap_burstmodel.mat
% Output: ~/IntermediateFiles/CI_burstdist.mat
%%
clear all; close all; clc;
addpath('../../Data');
addpath('../Functions/');
addpath('../../IntermediateFiles/');
%% Number of Lp distributions and model parameter combinations
nsamples = 1e3;
%% Load best LP to burst size linear model
load('LS_burst.mat');
%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));
%% Load linear confidence intervals
load('bootstrap_burstmodel.mat');
% Number of bootstrap we ran
n_bs = size(bs_linear,1);
% 95% parameter CI
idxs = round([n_bs*0.025 , n_bs*0.975]); 
% Get linear parameter range
for i = 1:size(bs_linear,2)
    this_sort = sort(bs_linear(:,i));
    param_CI_linear(i,1) = this_sort(idxs(1));
    param_CI_linear(i,2) = this_sort(idxs(2));
end

%%
means = zeros(nsamples,1);
cvs = zeros(nsamples,1);

for i = 1:nsamples
    % random sample latent period distribution
    thisCI = lhsu([CI(1,2),CI(2,2)],[CI(1,3),CI(2,3)],1)';
    [thisshape,thisscale] = gamma_fromstats(thisCI(1,1),thisCI(2,1));
    
    % random sample parameters of linear model
    thispars = lhsu([param_CI_linear(1,1),param_CI_linear(2,1)],[param_CI_linear(1,2),param_CI_linear(2,2)],1);

    % Get 1e5 random latent periods
    random_LP = gamrnd(thisshape,thisscale,1e5,1);
    % Compute burst sizes according to LP to burst size model
    prebursts_raw = burst_linear(thispars,random_LP);
    % Calculate mean and cv
    means(i) = mean(prebursts_raw);
    cvs(i) = std(prebursts_raw)/means(i);
end

%% Calculate mean and CV for best LP distribution and best parameter set
[thisshape,thisscale] = gamma_fromstats(CI(1,1),CI(2,1));
random_LP = gamrnd(thisshape,thisscale,1e8,1);
prebursts_raw = burst_linear(real_best_linear_expdecay,random_LP);
bestmean = mean(prebursts_raw);
bestcv = std(prebursts_raw)/bestmean;
%% Get 95% CI for mean and cv
% 95% parameter CI
idxs = round([nsamples*0.025 , nsamples*0.975]); 
% Get linear parameter range
sort_mean = sort(means);
sort_cv = sort(cvs);
%% Store data
CI_burstdist = [bestmean,sort_mean(idxs(1)),sort_mean(idxs(2));...
    bestcv, sort_cv(idxs(1)), sort_cv(idxs(2))];
%%
save('../../IntermediateFiles/CI_burstdist.mat','CI_burstdist');