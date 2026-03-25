% Code by: Marian Dominguez-Mirazo, 2026
% Get the visual fit for the cumulative burst size based on the best set of
% parameters for three models accounting for eclipse period, and the visual
% representation of the latent period to burst size distribution. 
% Required pre-computed files:
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/IntermediateFiles/LS_Burst_eclipse.mat
%%
clear all; close all; clc;
%%
% time after infection
xs = 0:0.25:14;
models = ["linear","mm","kannoly"];
% number of parameters per model
nparams = [2,3,4];
% v values to compute for linear model
v_vals = [.05,.25,.50,.75];
%%
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/');

%% Load Latent period distribution preciction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));

%% Load the best parameter fits for the three models
load('LS_Burst_eclipse.mat')

%% Get the best fit for all models
thebest_betaeff = zeros(numel(xs),3);
thebest_beta = zeros(numel(xs),3);

% linear parameters
pars_linear = [0,bestlinear(2:nparams(1))];
emean_linear = bestlinear(1);
v_linear = bestlinear(end);
[CIdist_e_linear,CIdist_pe_linear] = splitDistr_eclipse(CI,v_linear,emean_linear);
% mm parameters 
pars_mm = [0,bestmm(2:nparams(2))];
emean_mm = bestmm(1);
v_mm = bestmm(end);
[CIdist_e_mm,CIdist_pe_mm] = splitDistr_eclipse(CI,v_mm,emean_mm);
% kannoly parameters
pars_kannoly = [0,bestkannoly(2:nparams(3))];
emean_kannoly = bestkannoly(1);
v_kannoly = bestkannoly(end);
[CIdist_e_kannoly,CIdist_pe_kannoly] = splitDistr_eclipse(CI,v_kannoly,emean_kannoly);

% Calculate for linear 
thebest_betaeff(:,1) = calculate_betaeff_eclipse(xs,CIdist_pe_linear,CIdist_e_linear,pars_linear,'linear');
thebest_beta(:,1) = LP_burstdistr_eclipse(xs,CIdist_e_linear,pars_linear,'linear');
% Calculate for mm
thebest_betaeff(:,2) = calculate_betaeff_eclipse(xs,CIdist_pe_mm,CIdist_e_mm,pars_mm,'mm');
thebest_beta(:,2) = LP_burstdistr_eclipse(xs,CIdist_e_mm,pars_mm,'mm');
% Calculate for kannoly
thebest_betaeff(:,3) = calculate_betaeff_eclipse(xs,CIdist_pe_kannoly,CIdist_e_kannoly,pars_kannoly,'kannoly');
thebest_beta(:,3) = LP_burstdistr_eclipse(xs,CIdist_e_kannoly,pars_kannoly,'kannoly');

%% Get the best fit of linear model for given values of v
for i = 1:numel(v_vals)
    var_percentage = v_vals(i);
    bestlinear_row = find(bestlinear_allv(:,nparams(1)+2) == var_percentage);
    % load parameters
    pars_linear = [0,bestlinear_allv(bestlinear_row,2:nparams(1))];
    emean_linear = bestlinear_allv(bestlinear_row,1);
    [CIdist_e_linear,CIdist_pe_linear] = splitDistr_eclipse(CI,var_percentage,emean_linear);
    % Calculate
    thebest_betaeff_linear(:,i) = calculate_betaeff_eclipse(xs,CIdist_pe_linear,CIdist_e_linear,pars_linear,'linear');
    thebest_beta_linear(:,i) = LP_burstdistr_eclipse(xs,CIdist_e_linear,pars_linear,'linear');
end

%% Save data
save('../../IntermediateFiles/visual_burstfit_eclipse.mat','xs','v_vals','models','thebest_betaeff','thebest_beta','thebest_betaeff_linear','thebest_beta_linear');