% Code by: Marian Dominguez-Mirazo, 2025
% This script prepares visual CI for burst size fits by calculating the
% expected cumulative effective burst size using a random sample of
% parameters obtained from bootstrapping the original data
% It requires the pre-computed files: 
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/IntermediateFiles/bootstrap_burstmodel.mat
% Output: ~/IntermediateFiles/visualCI_burstfit.mat

clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/')

xs = 0.1:0.1:14;

%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));

%% Load bootstrap model predictions
load('bootstrap_burstmodel.mat');

%% Obtain curves for random parameter combinations within CI

% Number of randomly sample parameter combinations iterations
nrand = 1e3;
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

% Get mm parameter range
n_bs = size(bs_mm,1);
for i = 1:size(bs_mm,2)
    this_sort = sort(bs_mm(:,i));
    param_CI_mm(i,1) = this_sort(idxs(1));
    param_CI_mm(i,2) = this_sort(idxs(2));
end

% Get kannoly parameter range
n_bs = size(bs_kannoly,1);
for i = 1:size(bs_kannoly,2)
    this_sort = sort(bs_kannoly(:,i));
    param_CI_kannoly(i,1) = this_sort(idxs(1));
    param_CI_kannoly(i,2) = this_sort(idxs(2));
end

%%
% Randomly sample parameters from parameter range
betaeff_CI = zeros(nrand,numel(xs),3);
beta_CI = zeros(nrand,numel(xs),3);
for i = 1:nrand
    % Randomly sample LP distribution
    thisCI = lhsu([CI(1,2),CI(2,2)],[CI(1,3),CI(2,3)],1)';

    % Linear
    % Effective burst size
    % Randomly sample parameters
    this_params = lhsu(param_CI_linear(:,1),param_CI_linear(:,2),1);
    % eff burst size
    betaeff_CI(i,:,1) = calculate_betaeff(xs,thisCI,this_params,'linear');
    % Remove lines gone wrong
    if sum(betaeff_CI(i,:,1)<0)
        betaeff_CI(i,:,1) = NaN;
    end
    % burst size function 
    beta_CI(i,:,1) = burst_linear(this_params,xs);


    % mm
    % effective burst size
    this_params = lhsu(param_CI_mm(:,1),param_CI_mm(:,2),1);
    betaeff_CI(i,:,2) = calculate_betaeff(xs,CI,this_params,'mm');
    if sum(betaeff_CI(i,:,2)<0)
        betaeff_CI(i,:,2) = NaN;
    end
    % burst size function 
    beta_CI(i,:,2) = burst_mm(this_params,xs);
    if sum(beta_CI(i,:,2)<0)
        beta_CI(i,:,2) = NaN;
    end

    % kannoly
    % effective burst size
    this_params = lhsu(param_CI_kannoly(:,1),param_CI_kannoly(:,2),1);    
    betaeff_CI(i,:,3) = calculate_betaeff(xs,thisCI,this_params,'kannoly');
    % burst size function 
    beta_CI(i,:,3) = burst_kannoly(this_params,xs);

    if sum(beta_CI(i,:,3)<0) || sum(beta_CI(i,:,3)>120) || sum(betaeff_CI(i,:,3)<0)
        betaeff_CI(i,:,3) = NaN;
        beta_CI(i,:,3) = NaN;
    end
end

%%
% Get contour values to plot
minfill_betaeff = squeeze(min(betaeff_CI,[],1));
maxfill_betaeff = squeeze(max(betaeff_CI,[],1));

minfill_beta = squeeze(min(beta_CI));
maxfill_beta = squeeze(max(beta_CI));

%%
save('../../IntermediateFiles/visualCI_burstfit.mat','minfill_betaeff','maxfill_betaeff',...
    'minfill_beta','maxfill_beta');