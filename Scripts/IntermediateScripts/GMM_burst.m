% Code by Marian Dominguez-Mirazo, 2025
% This script performs GMM clustering to predict the individual bursts 
% that ocurred at half hour intervals througout the single-cell experiment
% It requires the pre-computed files: 
% Output: ~/IntermediateFiles/GMM_burst.mat
%%
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/');

%% Load burst size data
% Load file
file = 'SingleCellData.csv';
tab = readtable(file, 'ReadVariableNames', true);
tab = table2array(tab(:,3:end));
% Take only burst sizes larger than 1
tab(tab<2) = NaN;
data = tab;

% Get sampling times (they are the file header)
read_times = table2array(readtable(file));
sample_time = read_times(1,3:end);

% Get moments of effective burst size at sample timepoints
burstmean = nanmean(tab,1);
burststd = nanstd(tab,1);

%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));
[thisshape,thisscale] = gamma_fromstats(CI(1,1),CI(2,1));

%%
GMM_burst = [];
% Loop through sample points
% We start at 5.5 hr (i=4) because there's not enough burst sizes to cluster before that
for i = 4:16
    % Get burst sizes at sample timepoint
    sizes = data(~isnan(data(:,i)),i);
    % Because we are assuming 2 components, we need at least 3 data
    % points
    if numel(sizes)>2
        %% Set known parameters
        % Expected weigths of each burst size distribution drawn from LP distribution
        % CDF of main time point
        p2 = gamcdf(sample_time(i),thisshape,thisscale);
        % CDF of previous time point normalized to main timepoint CDF
        x1 = gamcdf(sample_time(i-1),thisshape,thisscale) / p2;
        % Probability of the half hour interval
        x2 = 1 - x1;

        % Observed mean and standard deviation for the first component
        mu_known = burstmean(i-1);
        sigma_known = burststd(i-1);

        %% Set up initial guess for other parameters
        % Initialize the second component using the whole distribution
        mu_init = burstmean(i);  % Initialize mean
        sigma_init = burststd(i); % Initialize the covariance to the std of the data
        
        % Initialize weigths
        pi_init = [x1; x2];  % Equal mixing weights for two components

        %% Initialize GMM parameters
        % Initial structure
        initParams = struct();
        % Initial means (known mean for the first component, approximation for the second)
        initParams.mu = [mu_known; mu_init];
        % Initial covariances (known variance for the first component, approximation for the second)
        initParams.Sigma = cat(3, sigma_known^2, sigma_init^2);
        % Initial weights per component
        initParams.PComponents = [x1; x2];

        %% Fit the GMM 
        options = statset('MaxIter', 1000);
        gmm = fitgmdist(sizes, 2, 'Start', initParams, 'Options', options,'RegularizationValue',1e-6);

        %% Force first cluster and WEIGHTS to be what we observe
        gm = gmdistribution([mu_known,gmm.mu(2)]',gmm.Sigma,[x1,x2]);

        %% Clustering
        % Cluster burst sizes
        labels = cluster(gm, sizes);
        this_sizes = sizes(labels == 2);
        % Get the average time for lysis events to occur
        this_time = mean([sample_time(i-1),sample_time(i)]);

        %% Store clustered burst sizes
        GMM_burst = [GMM_burst;[repmat(this_time,size(this_sizes)),this_sizes]];
    end
end

%% Save data
save('../../IntermediateFiles/GMM_burst.mat','GMM_burst');