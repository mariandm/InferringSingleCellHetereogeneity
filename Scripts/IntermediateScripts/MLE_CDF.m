% Code by Marian Dominguez-Mirazo, 2025
% This script finds the Maximum Likelihood Estimate for the
% latent period distribution based on data obtained via the
% 'single-cell lysis detection protocol'
% It requires the precomputed file: ~/IntermediateFiles/SingleCellInfections.csv
% Output: ~/IntermediateFiles/MLE_CDF.mat
%%
clc; clear all; close all;
addpath('../../IntermediateFiles/');
%% 
file = 'SingleCellInfections.csv';
% Column order is: (1) sample time, (2) replicate, (3) number of viable plates, (4) number
% of infected cells (plaque count > 0), (5) number of lysed cells (plaque count > 1)
exp = csvread(file);
% Remove NaN columns
exp = exp(~isnan(exp(:,3)),:);

%% From multiple initial conditions, find the minimum 
% Search across multiple models

nrand = 1e3; % Number of initial conditions
% Store data
s_gamma = zeros(nrand,3);
s_gamma_pl = zeros(nrand,4);
s_logn = zeros(nrand,3);
s_logn_pl = zeros(nrand,4);
% x0 range for distribution mean
Trange = [4,10];
% x0 range for CV
CVrange = [0.05,0.4];
% x0 range for plateau 
PLrange = [0.4,1];
% Randomly sample parameter values
randsamp = lhsu([Trange(1),CVrange(1),PLrange(1)],...
        [Trange(2),CVrange(2),PLrange(2)],...
        nrand);

for i = 1:nrand
    % Get an x0 sample
    x0 = randsamp(i,:);

    % Minimum search for all models

    % Gamma no plateau
    % Search functions
    [this_s_gamma, this_fval_gamma] = fminsearch(@(x)Likelihood_CDF(exp,x,0,"gamma"),x0(1:2));
    % Save output
    s_gamma(i,1:2) = this_s_gamma;
    s_gamma(i,3) = this_fval_gamma;
    
    % Gamma with plateau 
    % Search functions
    [this_s_gammapl, this_fval_gammapl] = fminsearch(@(x)Likelihood_CDF(exp,x,1,"gamma"),x0);
    % Save output
    s_gamma_pl(i,1:3) = this_s_gammapl;
    s_gamma_pl(i,4) = this_fval_gammapl;

    % Logn no plateau
    % Search functions
    [this_s_logn, this_fval_logn] = fminsearch(@(x)Likelihood_CDF(exp,x,0,"lognormal"),x0(1:2));
    % Save output
    s_logn(i,1:2) = this_s_logn;
    s_logn(i,3) = this_fval_logn;
    
    % Logn with plateau 
    % Search functions
    [this_s_lognpl, this_fval_lognpl] = fminsearch(@(x)Likelihood_CDF(exp,x,1,"lognormal"),x0);
    % Save output
    s_logn_pl(i,1:3) = this_s_lognpl;
    s_logn_pl(i,4) = this_fval_lognpl;
end

%% Manually check differences between minimums 
% Q: is there a global minimum? 
% A: Yes
[min(s_gamma,[],1); max(s_gamma,[],1)]
[min(s_gamma_pl,[],1); max(s_gamma_pl,[],1)]
[min(s_logn,[],1); max(s_logn,[],1)]
[min(s_logn_pl,[],1); max(s_logn_pl,[],1)]

%% Get Maximum Likelihood Estimates (MLE)
% We are using min() because we multiplied the likelihood by -1 
% to use the minimization function fminsearch

% Gamma no plateau
[A_gamma,I_gamma] = min(s_gamma(:,end));
MLE_gamma = s_gamma(I_gamma,1:2);
% Gamma with plateau 
[A_gamma_pl,I_gamma_pl] = min(s_gamma_pl(:,end));
MLE_gamma_pl = s_gamma_pl(I_gamma_pl,1:3);
% Logn no plateau
[A_logn,I_logn] = min(s_logn(:,end));
MLE_logn = s_logn(I_logn,1:2);
% Logn with plateau
[A_logn_pl,I_logn_pl] = min(s_logn_pl(:,end));
MLE_logn_pl = s_logn_pl(I_logn_pl,1:3);

%% Get likelihood for grid of values to calculate CI

% Range of parameter grid
Ts = Trange(1):0.1:Trange(2);
cvs = CVrange(1):0.01:CVrange(2);
pls = PLrange(1):0.01:PLrange(2);

% Combination of parameters
x_combns = table2array(combinations(Ts,cvs));
x_combns_pl = table2array(combinations(Ts,cvs,pls));

% Calculate maximum likelihoods
ll_gamma = Likelihood_CDF(exp,x_combns,0,"gamma");
ll_gamma_pl = Likelihood_CDF(exp,x_combns_pl,1,"gamma");
ll_logn = Likelihood_CDF(exp,x_combns,0,"lognormal");
ll_logn_pl = Likelihood_CDF(exp,x_combns_pl,1,"lognormal");

% Add information to the parameters tested
grid_ll_gamma = [x_combns,ll_gamma'];
grid_ll_gamma_pl = [x_combns_pl,ll_gamma_pl'];
grid_ll_logn = [x_combns,ll_logn'];
grid_ll_logn_pl = [x_combns_pl,ll_logn_pl'];

%% Calculate CIs
% Using ML_CI function
CI_gamma = [[MLE_gamma,1]', ML_CI(grid_ll_gamma,0,A_gamma)];
CI_gamma_pl = [MLE_gamma_pl', ML_CI(grid_ll_gamma_pl,1,A_gamma_pl)];
CI_logn = [[MLE_logn,1]', ML_CI(grid_ll_logn,0,A_logn)];
CI_logn_pl = [MLE_logn_pl', ML_CI(grid_ll_logn_pl,1,A_logn_pl)];

CI(1,:,:) = CI_gamma;
CI(2,:,:) = CI_gamma_pl;
CI(3,:,:) = CI_logn;
CI(4,:,:) = CI_logn_pl;

%% Save CI
save('../../IntermediateFiles/MLE_CDF.mat','CI');