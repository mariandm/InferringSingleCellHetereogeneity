% Code by Marian Dominguez-Mirazo, 2025
% Caution: This script takes a bit to run (~30 min)
% It finds the parameters for a latent period to burst size
% function (linear, kannoly or hill model) that best describe teh effective
% burst size data from the sigle-cell experiments, entangled with the
% predicted latent period distribution
% It requires the pre-computed files: 
% - ../../IntermediateFiles/MLE_CDF.mat
% - ../../IntermediateFiles/burstMean_perReplicate.csv
% Output: '../../IntermediateFiles/LS_Burst.mat'

%%
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/')
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

%% Find best fits for each model 
nrand = 100; % Number of initial search points
%% Linear model
real_best_linear_expdecay = burst_bestfit(efburst_per_replicate_og,sample_time,CI,'linear',nrand);
%% Kannoly model 
real_best_kannoly_expdecay = burst_bestfit(efburst_per_replicate_og,sample_time,CI,'kannoly',nrand);
%% Hill function
real_best_mm_expdecay = burst_bestfit(efburst_per_replicate_og,sample_time,CI,'mm',nrand);
%% Save
save('../../IntermediateFiles/LS_Burst.mat','real_best_linear_expdecay','real_best_kannoly_expdecay','real_best_mm_expdecay');