% Code by Marian Dominguez-Mirazo, 2025
% Caution: This script takes several hours to run
% Bootstrap cumulative effective burst size data and calculate the Least
% Square parameter set that best describes the data given 3 models: linear,
% kannoly, and hill function
% Requires pre-computed files:
% - GMM_burst.mat
% - MLE_CDF.mat
% - burstMean_perReplicate.csv
% Output: ~/IntermediateFiles/bootstrap_burstmodel.mat
%%
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/')

%% Load GMM burst sizes
load('GMM_burst.mat');

%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));
[thisshape,thisscale] = gamma_fromstats(CI(1,1),CI(2,1));

%% Get effective burst size per replicate
% Load file
file = 'SingleCellData.csv';
tab = readtable(file, 'ReadVariableNames', true);
% Get sampling times (they are the file header)
read_times = table2array(readtable(file));
sample_time = read_times(1,3:end);

%% Bootstrap the parameter search using the uncertainty for the distribution

n_bs = 1e3; % number of bootstraps
nrand = 50; % number of initial parameter conditions for fminsearch

% store data
bs_linear = zeros(n_bs,3);
bs_kannoly = zeros(n_bs,5);
bs_mm = zeros(n_bs,4);

CI_og = CI; % original CI for latent period distribution

%% Linear bootstrap
for bs = 1:n_bs

    disp(bs);

    % Random sample data points
    rs = [randsample(30,30,true);...
    randsample(30,30,true) + 30;...
    randsample(30,30,true) + 60;...
    randsample(30,30,true) + 90];
    
    % Calculate per replicate mean for the bootstraps
    thistab = tab(rs,:);
    efburst_per_replicate = burstmean_perrep(thistab);

    % Get a LP distribution within the CI
    CI = lhsu([CI_og(1,2),CI_og(2,2)],[CI_og(1,3),CI_og(2,3)],1)';

    best_linear = burst_bestfit(efburst_per_replicate,sample_time,CI,'linear',nrand);
    bs_linear(bs,:) = best_linear;

end

%% Kannoly bootstrap
for bs = 1:n_bs

    disp(bs);

    % Random sample data points
    rs = [randsample(30,30,true);...
    randsample(30,30,true) + 30;...
    randsample(30,30,true) + 60;...
    randsample(30,30,true) + 90];
    
    % Calculate per replicate mean for the bootstraps
    thistab = tab(rs,:);
    efburst_per_replicate = burstmean_perrep(thistab);

    % Get a LP distribution within the CI
    CI = lhsu([CI_og(1,2),CI_og(2,2)],[CI_og(1,3),CI_og(2,3)],1)';

    best_kannoly = burst_bestfit(efburst_per_replicate,sample_time,CI_og,'kannoly',nrand);
    bs_kannoly(bs,:) = best_kannoly;
end

%% mm bootstrap
for bs = 1:n_bs

    disp(bs);

    % Random sample data points
    rs = [randsample(30,30,true);...
    randsample(30,30,true) + 30;...
    randsample(30,30,true) + 60;...
    randsample(30,30,true) + 90];
    
    % Calculate per replicate mean for the bootstraps
    thistab = tab(rs,:);
    efburst_per_replicate = burstmean_perrep(thistab);

    % Get a LP distribution within the CI
    CI = lhsu([CI_og(1,2),CI_og(2,2)],[CI_og(1,3),CI_og(2,3)],1)';

    best_mm = burst_bestfit(efburst_per_replicate,sample_time,CI_og,'mm',nrand);
    bs_mm(bs,:) = best_mm;

end
%%
save('../../IntermediateFiles/bootstrap_burstmodel.mat','bs_linear','bs_kannoly','bs_mm');



