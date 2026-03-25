% Code by marian Dominguez-Mirazo, 2026
% This script plots three different production period to burst size models,
% assuming the latent period can be divided into a eclipse period and a
% production period. It also plots the linear model of production period to
% burst size using different values of parameter nu, the proportion of
% eclipse period variability that explains latent period variability. 
% It requires the precomputed files:
% - ~/IntermediateFiles/visual_burstfit_eclipse.mat
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/IntermediateFiles/GMM_burst.mat
% - ~/IntermediateFiles/burstMean_perReplicate.csv

%%
clear all; close all; clc;
%%
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/');

%% Load prerun best fits for burst models
load('visual_burstfit_eclipse.mat');
%% Load predited LP dsitribution
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));
%% Load predictd burst size distributions for linear model different values of v
load('burstDistribution_predicted_eclipse.mat');
%% Load GMM burst sizes
load('GMM_burst.mat');
% Average across time
times = unique(GMM_burst(:,1));
average_GMM = zeros(size(times));
for t = 1:numel(times)
    average_GMM(t) = mean(GMM_burst(GMM_burst(:,1)==times(t),2));
end
average_GMM = [times,average_GMM];
%% Get effective burst size per replicate
% Load file
file = 'burstMean_perReplicate.csv';
tab = table2array(readtable(file, 'ReadVariableNames', false));
sample_time = tab(1,:);
efburst_per_replicate_og = tab(2:end,:);

% Get mean per time
efburst_og = nanmean(efburst_per_replicate_og);

%% Initialize figure
figure('Position',[10,10,800,1200]);
tiledlayout(3,2);
markershapes = ["o","square","^","v"];
pcolors = ["#000000","#7879FF","#BFBFFF","#E5E5FF"];
plinestyles = ["--","-","-","-","-"];
colors = ["#000000","#D95319","#349a04"];
linestyles = ["--","-","-"];

%% Second column specifics

for i = [2,4]
    nexttile(i);
    % Plot detangled burst sizes
    scatter(GMM_burst(:,1),GMM_burst(:,2),'filled','MarkerFaceColor','#999999',...
        'Marker','d','SizeData',50); hold on;
    
    xlabel('Latent period, hr');
    ylabel('Expected Burst Size');
    xlim([0,max(GMM_burst(:,1))]);
    ylim([0,250]);
end

%% Second row linear model
% Plot fit to burst temporal data
nexttile(3);
for i = 4:-1:1 % just p_vals = [0.0500, 0.2500, 0.5000, 0.7500]
    plot(xs,thebest_betaeff_linear(:,i),'Color',pcolors(i),...
        'LineWidth',1.5,'LineStyle',plinestyles(i)); hold on;
end

% Plot corresponding fit on GMM predictions
nexttile(4);
for i = 1:4 % just p_vals = [0.0500, 0.2500, 0.5000, 0.7500]
    plot(xs,thebest_beta_linear(:,i),'Color',pcolors(i),...
        'LineWidth',1.5,'LineStyle',plinestyles(i)); hold on;
end

%% First row best linear vs best mm and best kannoly

nexttile(1);
for i = 3:-1:1
    plot(xs',thebest_betaeff(:,i),'Color',colors(i),...
        'LineWidth',1.5,'LineStyle',linestyles(i)); hold on;
end

nexttile(2);
for i = 1:3
    plot(xs,thebest_beta(:,i),'Color',colors(i),...
        'LineWidth',1.5,'LineStyle',linestyles(i)); hold on;
end

%% First column specifics

for i = [1,3]
    nexttile(i);
    % Plot mean across all replicates
    % In this order for legend purposes
    scatter(sample_time,efburst_og,'filled','MarkerFaceColor','k'); hold on;
    
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
end

%% Burst size distribution

nexttile(5,[1 2]);
histogram(thisdata,'BinWidth',10,'FaceColor','#999999','EdgeColor','#989898'); hold on;

for i = 1:4 
    plot(mh_minpoints_raw_all(i,2:end),mh_y_new_raw_all(i,2:end),'Color',pcolors(i),'LineWidth',1.5,...
        'LineStyle',plinestyles(i));
end

xlim([0,250]);
ylim([0,30]);

xlabel('Burst size');
ylabel('Counts');
legend('Single-cell experiment',...
    'v = 0.05',...
    'v = 0.25',...
    'v = 0.5',...
    'v = 0.75');
legend('boxoff');


%% Aesthetics
nexttile(1);
legend('Logistic growth','Hill function','Linear model',...
    'Replicates mean','Replicates',...
    'Direction','reverse');
legend('boxoff');
legend('Location','northwest');

nexttile(2);
legend('Inferred lysis events');
legend('boxoff');
legend('Location','northwest');

nexttile(3);
legend('v = 0.75', 'v = 0.5', 'v = 0.25', 'v = 0.05',...
    'Replicates mean','Replicates',...
    'Direction','reverse');
legend('boxoff');
legend('Location','northwest');

nexttile(4);
legend('Inferred lysis events');
legend('boxoff');
legend('Location','northwest');

for i =1:5
    nexttile(i);
    box off;
    ax=gca;
    ax.FontSize=13;
    set(gca,'FontName','Latin Modern Sans');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
end

%%
saveas(gcf,'../../Figures/SupplFig10_BurstSizeFunction_eclipse.svg');