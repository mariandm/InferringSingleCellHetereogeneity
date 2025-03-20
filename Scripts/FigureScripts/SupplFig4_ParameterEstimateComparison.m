% Code by Marian Dominguez-Mirazo, 2025
% Compare parameter estimations (LP mean, CV and burst size) between the
% population-level estimation and the single-cell derived estimates
% Requires pre-computed files:
% - ~/IntermediateFiles/SingleCellData.csv
% - ~/MLE_CDF.mat
% - ~/MCMCChain_020325.csv
% Output: ~/Figures/SupplFig4_ParameterEstimateComparison.svg
%%
clear all; close all; clc;
addpath('../../IntermediateFiles/');
addpath('../../Data');
addpath('../Functions')
%% Load experiment data
file = 'SingleCellData.csv';
singlecelldata = readtable(file, 'ReadVariableNames', true);
singlecelldata = table2array(singlecelldata(:,3:end));
% Get burst sizes larger from more than 1 plaque 
singlecelldata(singlecelldata<2) = NaN;
% Times are in header
read_times = table2array(readtable(file));
sample_time = read_times(1,3:end);

%% Load  MLE and CI from single cell
load('MLE_CDF.mat');
% The best model (gamma with plateau) is stored in index 2
i=2;
bestCI = squeeze(CI(2,:,:));

%% Read viral chain
% Viral chain
file = 'MCMCChain_020325.csv';
tab = readtable(file, 'ReadVariableNames', false);
tab = table2array(tab);
phi_chain = tab(:,3);
beta_chain = tab(:,4);
eta_chain = tab(:,5);
cv_chain = tab(:,6);

% Get 95% CI chains
phi_chain_95 = phi_chain(find((phi_chain > quantile(phi_chain,0.025)) .* (phi_chain < quantile(phi_chain,0.975))));
beta_chain_95 = beta_chain(find((beta_chain > quantile(beta_chain,0.025)) .* (beta_chain < quantile(beta_chain,0.975))));
eta_chain_95 = eta_chain(find((eta_chain > quantile(eta_chain,0.025)) .* (eta_chain < quantile(eta_chain,0.975))));
cv_chain_95 = cv_chain(find((cv_chain > quantile(cv_chain,0.025)) .* (cv_chain < quantile(cv_chain,0.975))));

%% Load burst size model fit
load('LS_burst.mat');
betaeff = calculate_betaeff(12,bestCI,real_best_linear_expdecay,'linear');
load('visualCI_burstfit.mat');
betaeff_CI = [minfill_betaeff(end,1),maxfill_betaeff(end,1)];

%% Choosing a time to stop counting burst size
% By observing if there are statistical significance of burst size
% distribution at different times when compared to the last two sampling
% times (11.5 and 12 hr)

% Concatenate the last two times
lasttime = [singlecelldata(~isnan(singlecelldata(:,end)),end);singlecelldata(~isnan(singlecelldata(:,end-1)),end-1)];
% index list to substract from end
is = 2:15;
% Saving null hypothesis acceptance or rejection
hs = nan(numel(is,1));
% Iterate at previous times and save whether null hypothesis was rejected
for i = is
    % Get burst sizes of previous time
    this_time = singlecelldata(~isnan(singlecelldata(:,end-i)),end-i);
    % Wilcoxon rank-sum test
    [p,h] = ranksum(lasttime,this_time);
    % Save if null hypothesis accepted (0), or rejected (1)
    hs(i-1) = h;
end
% Get the times for which the null hypothsesis was accepted
sample_time(numel(sample_time) - is(~logical(hs)));
% They go way back to 7 hrs, so assuming 9 hrs should be fine

%% Calculate single-cell derived burst size mean and CI from data bootstraps
% Starting from sampling time point = 9hr
idx = find(sample_time >= 9);
% Reshape data to be a vector and remove NaNs
newtab = reshape(singlecelldata(:,idx),1,numel(singlecelldata(:,idx)));
newtab = newtab(~isnan(newtab));
burst_mean = mean(newtab);
nboots = 1e4; % Number of bootstraps 
burst_boot = zeros(nboots,1);
% Get a burst size mean
for i = 1:nboots
    burst_boot(i) = mean(randsample(newtab,numel(newtab),true));
end
% Find 95% CI
burst_boot_CI = [quantile(burst_boot,0.025), quantile(burst_boot,0.975)];

%% Start Figure
figure('Position',[10,10,1000,400]);
tiledlayout(1,3);

%% A) Plot LP mean comparison
nexttile(1);
% Single-cell level mean latent period
errorbar(1,bestCI(1,1),bestCI(1,1)-bestCI(1,2),bestCI(1,1)-bestCI(1,3),...
    'LineWidth',1.5,'Color','b'); hold on;
scatter(1,bestCI(1,1),'filled','b');
% Pop-level  mean latent period CI
errorbar(2,mean(1 ./ eta_chain_95),...
    mean(1 ./ eta_chain)-min(1 ./ eta_chain_95), ...
    mean(1 ./ eta_chain)-max(1 ./ eta_chain_95),...
    'LineWidth',1.5,'Color','r'); hold on;
scatter(2,mean(1 ./ eta_chain),'filled','r');
xlim([0.5,2.5]);
ylim([6,9]);
ylabel('Latent period mean, hr');

%% B) Plot LP CV comparison
nexttile(2);
% Single-cell level CV
errorbar(1,bestCI(2,1),bestCI(2,1)-bestCI(2,2),bestCI(2,1)-bestCI(2,3),...
    'LineWidth',1.5,'Color','b'); hold on;
scatter(1,bestCI(2,1),'filled','b');
% Pop-level CV CI
errorbar(2,mean(cv_chain_95),...
    mean(cv_chain)-min(cv_chain_95), ...
    mean(cv_chain)-max(cv_chain_95),...
    'LineWidth',1.5,'Color','r'); hold on;
scatter(2,mean(cv_chain),'filled','r');
xlim([0.5,2.5]);
ylim([0.05,0.35]);
ylabel('Coefficient of Variation (CV)');

%% C) Plot burst size comparison
nexttile(3);
% Single-cell level burst size
errorbar(1,burst_mean,burst_mean-burst_boot_CI(1),burst_mean-burst_boot_CI(2),...
    'LineWidth',1.5,'Color','b'); hold on;
scatter(1,burst_mean,'filled','b');
% Pop-level burst size
errorbar(2,mean(beta_chain_95),...
    mean(beta_chain)-min(beta_chain_95), ...
    mean(beta_chain)-max(beta_chain_95),...
    'LineWidth',1.5,'Color','r'); hold on;
scatter(2,mean(beta_chain),'filled','r');
ylabel('Effective burst size');

legend('Single-cell 95% CI','Single-cell best estimate', ...
    'Population-level 95% CI', 'Population-level best estimate');

xlim([0.5,2.5]);
ylim([40,80]);
legend('Location','eastoutside');
legend('boxoff');


%% Aesthetics 
for i = 1:3
    nexttile(i);
    xticklabels([]);
    xticks([]);
    ax=gca;
    ax.FontSize=17;
    set(gca,'FontName','Latin Modern Roman')
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
end

%%
saveas(gcf,'../../Figures/SupplFig4_ParameterEstimateComparison.svg');