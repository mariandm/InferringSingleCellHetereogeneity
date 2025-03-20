% Code by: Marian Dominguez-Mirazo, 2025
% This script creates all panels in Figure 3
% It requires the pre-computed files: 
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/IntermediateFiles/LS_burst.mat
% - ~/IntermediateFiles/GMM_burst.mat
% - ~/IntermediateFiles/burstMean_perReplicate.csv
% - ~/IntermediateFiles/visualCI_burstfit.mat
% - ~/IntermediateFiles/burstDistribution_predicted.mat
% Output: ~/Figures/Figure3_BurstSizeFunction.svg
%%
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/')
%%
xs = 0.1:0.1:14;

%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));
[thisshape,thisscale] = gamma_fromstats(CI(1,1),CI(2,1));

%% Load best burst size model predictions
load('LS_burst.mat');
thebest_betaeff = calculate_betaeff(xs,CI,real_best_linear_expdecay,'linear');
thebest_beta = burst_linear(real_best_linear_expdecay,xs);

%% Load GMM burst sizes
load('GMM_burst.mat');

%% Get effective burst size per replicate
% Load file
file = 'burstMean_perReplicate.csv';
tab = table2array(readtable(file, 'ReadVariableNames', false));
sample_time = tab(1,:);
efburst_per_replicate_og = tab(2:end,:);

% Get mean per time
efburst_og = nanmean(efburst_per_replicate_og);

%% Load visual CI
load('visualCI_burstfit.mat')

%% load burst distributions
load('burstDistribution_predicted.mat');

%% Initialize figure
figure('Position',[10,10,800,800]);
tiledlayout(2,2);
markershapes = ["o","square","^","v"];

%% A) Plot the cumulative effective burst size fit

nexttile(1);
% Confidence Interval
fill([xs';flipud(xs')], [minfill_betaeff(:,1);flipud(maxfill_betaeff(:,1))],'b', ...
    'FaceColor','#D95319',...
    'EdgeColor','None','FaceAlpha',0.3); hold on;
% Best estimate
plot(xs,thebest_betaeff,'Color','#D95319','LineWidth',1.5);

% Plot mean across all replicates
% In this order for legend purposes
scatter(sample_time,efburst_og,'filled','MarkerFaceColor','k');

% Plot single cell experiment data
for j = 1:size(efburst_per_replicate_og,1)
    scatter(sample_time,efburst_per_replicate_og(j,:),'Filled',...
        'MarkerFaceColor','#999999',...
        'SizeData',50,...
        'Marker',markershapes(j));
end

% Plot mean across all replicates
scatter(sample_time,efburst_og,'filled','MarkerFaceColor','k');

xlabel('Time after infection, hr')
ylabel({'Effective Burst Size of';'Cumulative Lysis Time'});
xlim([0,12.5]);
ylim([0,120]);

legend('CI','Best fit','Replicates mean','Replicates','Direction','reverse');
legend('boxoff');
legend('Location','northwest');

%% B) Plot the LP to burst size relationship

% Plot burst size function
nexttile(2);
% Confidence Interval
fill([xs';flipud(xs')], [minfill_beta(:,1);flipud(maxfill_beta(:,1))],'b', ...
    'FaceColor','#D95319',...
    'EdgeColor','None','FaceAlpha',0.3); hold on;
% Best estimate
plot(xs,thebest_beta,'Color','#D95319','LineWidth',1.5);

% Plot detangled burst sizes
scatter(GMM_burst(:,1),GMM_burst(:,2),'filled','MarkerFaceColor','#999999',...
    'Marker','d','SizeData',50);

xlabel('Latent period, hr');
ylabel('Expected Burst Size');

xlim([0,max(GMM_burst(:,1))]);
legend('CI','Best fit','Inferred lysis events','Direction','reverse');
legend('boxoff');
legend('Location','northwest');

%% C) Plot the predicted burst size distribution
nexttile(3,[1 2]);

histogram(thisdata,'BinWidth',10,'FaceColor','#999999','EdgeColor','#989898'); hold on;
plot(mh_minpoints_raw(2:end),mh_y_new_raw(2:end),'Color','#D95319','LineWidth',1.5);
plot(mh_minpoints(2:end),mh_y_new(2:end),'Color','k','LineWidth',1.5,'LineStyle','--');

xlim([0,250])
xlabel('Burst size');
ylabel('Counts');
legend('Single-cell experiment',...
    'Expected burst size distribution from latent period variation',...
    'Sampling and adhesion correction');
legend('boxoff');
%% Aesthetics
for i =1:3
    nexttile(i);
    box off;
    ax=gca;
    ax.FontSize=18;
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
end

%%
saveas(gcf,'../../Figures/Figure3_BurstSizeFunction.svg');