% Code by Marian Dominguez-Mirazo, 2025
% This script takes an MCMC for viral trait prediction and calculates 
% 95% confidence intervals for the viral dynamics
% Output: ~/IntermediateFiles/visualCI_popLevel.mat

%% Plot population-level fit
clc; clear all; close all;
addpath('../../IntermediateFiles/');
addpath('../Functions/');
addpath('../../Data/');
%% Load data
% Viral timeseries
file = 'TimeseriesVirus.csv';
virus = readtable(file);
virus = virus{:,:};

%% Set host-only parameters
% Obtained from host parameter MCMC
pars.mu     = 0.04;    % Host growth, hr^-1
pars.K      = 3e8;    % Carrying capacity, CFU/ml

%% Read viral chain

% Viral chain
file = 'MCMCChain_020325.csv';
tab = readtable(file, 'ReadVariableNames', false);
tab = table2array(tab);
H0_chain = tab(:,1);
V0_chain = tab(:,2);
phi_chain = tab(:,3);
beta_chain = tab(:,4);
eta_chain = tab(:,5);
cv_chain = tab(:,6);

% Get 95% CI chains by getting quantiles
H0_chain_95 = H0_chain(find((H0_chain > quantile(H0_chain,0.025)) .* (H0_chain < quantile(H0_chain,0.975))));
V0_chain_95 = V0_chain(find((V0_chain > quantile(V0_chain,0.025)) .* (V0_chain < quantile(V0_chain,0.975))));
phi_chain_95 = phi_chain(find((phi_chain > quantile(phi_chain,0.025)) .* (phi_chain < quantile(phi_chain,0.975))));
beta_chain_95 = beta_chain(find((beta_chain > quantile(beta_chain,0.025)) .* (beta_chain < quantile(beta_chain,0.975))));
eta_chain_95 = eta_chain(find((eta_chain > quantile(eta_chain,0.025)) .* (eta_chain < quantile(eta_chain,0.975))));
cv_chain_95 = cv_chain(find((cv_chain > quantile(cv_chain,0.025)) .* (cv_chain < quantile(cv_chain,0.975))));

%% Get CI
nsteps = 1e4; % Number of parameter combinations to produce
t = 0:0.1:24;
options = odeset('AbsTol',1e-6,'RelTol',1e-6);

% Create parameter combinations
phi_chain_sample = randsample(phi_chain_95,nsteps,true);
beta_chain_sample = randsample(beta_chain_95,nsteps,true);
eta_chain_sample = randsample(eta_chain_95,nsteps,true);
cv_chain_sample = randsample(cv_chain_95,nsteps,true);
H0_chain_sample = randsample(H0_chain_95,nsteps,true);
V0_chain_sample = randsample(V0_chain_95,nsteps,true);

% Store data
Vsave_data = zeros(numel(t),nsteps);

% loop through chain steps
for i = 1:nsteps
    % Set parameters
    pars.phi    = 10^phi_chain_sample(i);   % Adsorption rate, ml/(CFUxhr)
    pars.beta   = beta_chain_sample(i);    % Burst size
    pars.eta = eta_chain_sample(i); % Lysis rate
    pars.cv = cv_chain_sample(i);
    pars.n = round(1/pars.cv^2 -1);
    pars.initS  = H0_chain_sample(i); % Initial host density, CFU/ml
    pars.initV  = V0_chain_sample(i); % Initial viral density, PFU/ml

    % initialize
    x0 = zeros(pars.n+3,1);
    x0(1) = pars.initS; x0(end) = pars.initV;
    % Run ODE
    [tsol,ysol] = ode45(@ODE_SEnIV_sink,t,x0,options,pars);
    % Save run
    Vsave_data(:,i) = ysol(:,end);
end
V_interval = [min(Vsave_data,[],2),max(Vsave_data,[],2)];

%%
save('../../IntermediateFiles/visualCI_popLevel.mat','V_interval');