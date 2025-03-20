% Code by Marian Dominguez-Mirazo, 2025
% This script plots the fits for different burst size models in the
% cumulative effective burst size and predicted latent period to burst size
% relationship
% Requires pre-computed files:
% - burstMean_perReplicate.csv
% - MLE_CDF.mat
% - LS_burst.mat
% - visualCI_burstfit.mat
% - GMM_burst.mat
% Output: ~/Figures/Suppl_BurstSizeFitModels_expdecay.svg
%%
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');
addpath('../../IntermediateFiles/');

%%
xs = 0.1:0.1:14;

%% Load effective burst size per replicate, and burst size mean per time
% Load file
file = 'burstMean_perReplicate.csv';
tab = table2array(readtable(file, 'ReadVariableNames', false));
% Get sampling times
sample_time = tab(1,:);
% get means
efburst_per_replicate_og = tab(2:end,:);

% Get mean per time
efburst_og = nanmean(efburst_per_replicate_og);

%% Load LP distribution prediction
load('MLE_CDF.mat')
CI = squeeze(CI(2,:,:));

%% Load and calculate best burst size model predictions
load('LS_burst.mat');

thebests_betaeff(1,:) = calculate_betaeff(xs,CI,real_best_linear_expdecay,'linear');
thebests_betaeff(2,:) = calculate_betaeff(xs,CI,real_best_mm_expdecay,'mm');
thebests_betaeff(3,:) = calculate_betaeff(xs,CI,real_best_kannoly_expdecay,'kannoly');

thebests_beta(1,:) = burst_linear(real_best_linear_expdecay,xs);
thebests_beta(2,:) = burst_mm(real_best_mm_expdecay,xs);
thebests_beta(3,:) = burst_kannoly(real_best_kannoly_expdecay,xs);

%% Load visual CI for plotting
load('visualCI_burstfit.mat')

%% Load GMM burst sizes
load('GMM_burst.mat');

%% Plot the three models with confidence intervals
figure('Position',[10,10,1200,800]);
tiledlayout(2,3);

titles = ["Linear model","Hill function","Kannoly model"]; % model order
markershapes = ["o","square","^","v"]; % for replicates

for i = 1:3
    % Plot beta effective
    nexttile(i);
    fill([xs';flipud(xs')], [minfill_betaeff(:,i);flipud(maxfill_betaeff(:,i))],'b', ...
        'FaceColor','#D95319',...
        'EdgeColor','None','FaceAlpha',0.3); hold on;
    plot(xs,thebests_betaeff(i,:),'Color','#D95319','LineWidth',1.5);

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
    title(titles(i));
    xlim([0,12.5]);

    % Plot burst size function
    nexttile(i+3);
    fill([xs';flipud(xs')], [minfill_beta(:,i);flipud(maxfill_beta(:,i))],'b', ...
        'FaceColor','#D95319',...
        'EdgeColor','None','FaceAlpha',0.3); hold on;
    plot(xs,thebests_beta(i,:),'Color','#D95319','LineWidth',1.5);
    % GMM detangled bursts at specific times
    scatter(GMM_burst(:,1),GMM_burst(:,2),'filled','MarkerFaceColor','#999999','Marker','diamond');

    xlim([0,11.5]);
end

nexttile(1);
ylabel({'Effective Burst Size of';'Cumulative Lysis Time'});
nexttile(4);
ylabel('Expected burst size');
nexttile(2);
xlabel('Time after infection, hr')
nexttile(5);
xlabel('Latent period, hr');

nexttile(3);
legend('CI','Best fit','Replicates mean','Replicates','Direction','reverse')
legend('boxoff');
legend('Location','eastoutside')

nexttile(6);
legend('CI','Best fit','Inferred lysis events','Direction','reverse')
legend('boxoff');
legend('Location','eastoutside')

% Aesthetics
for i=1:6
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
saveas(gcf,'../../Figures/SupplFig9_BurstSizeFitModels.svg');