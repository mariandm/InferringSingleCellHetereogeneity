% Code by Marian Dominguez-Mirazo, 2025
% Caution: This script takes ~20 min to run
% This script takes an incubation time (time to allow phage to infect 
% bacteria) for the single-cell lysis detection protocol and draws 
% random adsorption times and random latent periods, then calculates the
% observed average latent period and CV when including the adsorption times
% Required pre-computed file: ~/IntermediateFiles/MLE_CDF.mat
% Output: ~/IntermediateFiles/incubationError.mat

%%
clear all; close all; clc;
addpath('../Functions');
addpath('../../IntermediateFiles/');

%% Load predicted latent period distribution
load('MLE_CDF.mat');
gamma_CI = squeeze(CI(2,:,:));
%% 
S0 = 1.5e8; % Initial cell density in experiment
V0 = 1.5e8 * 1.5; % Initial phage density in experiment
params.phi = 1.6e-8; % adsorption rate, drawn from pop-level estimate

ncells = 1e4; % number of cells to simulate
nexp = 1e3; % number of experiments
%%
adsTimes = 5:1:30; % Incubation times to simulate (min)
% store resulting mean and CV
means = zeros(numel(adsTimes),1);
cvs = zeros(numel(adsTimes),1);

% for each simulated incubation time
for i = 1:numel(adsTimes)
    % Convert to hours
    expInfo.adsTime = adsTimes(i)/60;
    % for each experiment
    for j = 1:nexp

        %% Generate random adsorption times
        pars.phi = params.phi;
        pars.tf = expInfo.adsTime;
        random_adsTime = random_ads(rand(ncells,1),[S0,V0],pars);
    
        %% Generate random times of lysis for as many samples
        [thisshape,thisscale] = gamma_fromstats(gamma_CI(1,1),gamma_CI(2,1));
        random_LP = gamrnd(thisshape,thisscale,ncells,1);
    
        %% Get full distribution description
        full_time = random_adsTime + random_LP;
    
        thismean = mean(full_time);
        thiscv = std(full_time) / thismean;
    
        % Save
        means(i) = means(i) + thismean/nexp;
        cvs(i) = cvs(i) + thiscv/nexp;
    end
end

%%
save('../../IntermediateFiles/incubationError.mat','adsTimes','means','cvs');