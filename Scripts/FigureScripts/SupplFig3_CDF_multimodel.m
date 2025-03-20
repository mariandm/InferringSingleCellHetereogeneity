% Code by Marian Dominguez-Mirazo, 2025
% Plot fits for latent period CDF for four different models
% Requires pre-computed files:
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/SingleCellInfections.csv
% Output: ~/Figures/SupplFig3_CDFmodels.svg
%% 
clc; clear all; close all;
addpath('../../IntermediateFiles/');
addpath('../Functions/');

%% Load data
% MLE and CI
load('MLE_CDF.mat');
% Experiment data
file = 'SingleCellInfections.csv';
% Column order is: (1) sample time, (2) replicate, (3) number of viable plates, (4) number
% of infected cells (plaque count > 0), (5) number of lysed cells (plaque count > 1)
exp = csvread(file);

%% Get per time proportion of lysed
times = unique(exp(:,1));
means = [];
for i = 1:numel(times)
    idx = exp(:,1)==times(i);
    means(i) = nanmean(exp(idx,5) ./ exp(idx,4));
end

%% Visual CI
timepoints = 0:0.01:15;
% initialize storage
maxes = zeros(numel(timepoints),4);
mines = zeros(numel(timepoints),4);
MLEs = zeros(numel(timepoints),4);

% Model accounting for plateau correction (0:no, 1:yes)
model_pls = [0,1,0,1];
% Types of CDF distributions 
model_dist = ["gamma","gamma","lognormal","lognormal"];

% For each model 
for i = 1:4

    % Get parameter combinations given model 95% CI 
    Ts = CI(i,1,2):0.1:CI(i,1,3);
    cvs = CI(i,2,2):0.05:CI(i,2,3);
    pls = CI(i,3,2):0.05:CI(i,3,3);
    x_combns = table2array(combinations(Ts,cvs,pls));

    % Get CDF through time for different parameters
    ps = get_CDF(timepoints',x_combns,model_pls(i),model_dist(i));
    
    % Find edges of model CDF CI
    maxes(:,i) = max(ps,[],2);
    mines(:,i) = min(ps,[],2);

    % Make MLE prediction into CDF through time
    MLEs(:,i) = get_CDF(timepoints',CI(i,:,1),model_pls(i),model_dist(i));
end


%% Start figure
figure('Position',[10,10,900,600]);
tl = tiledlayout(2,2);
colors = ["b","b","#9F2B68","#9F2B68"];
markershapes = ["o","square","^","v"];

%%
% for each model 
for i=1:4
    nexttile(i);

    % Plot CI
    fill([timepoints';flipud(timepoints')],[mines(:,i);flipud(maxes(:,i))],"",...
        "FaceColor",colors(i),"EdgeColor","None","FaceAlpha",0.3); hold on;

    % Plot MLE
    plot(timepoints, MLEs(:,i), "Color",colors(i), "LineWidth",1.5);

    % Plot means across experiments 
    scatter(times,means,'Filled','k','SizeData',50);
    
    % Plot experiment data (per replicate)
    replicates = unique(exp(:,2));
    for j = 1:numel(replicates)
        idx = exp(:,2) == replicates(j);
        scatter(exp(idx,1),exp(idx,5)./exp(idx,4),'Filled',...
            'MarkerFaceColor','#999999',...
            'SizeData',50,...
            'Marker',markershapes(j));
    end

    % Plot means across experiments 
    % I plot it twice so the legend appears in the order I want
    scatter(times,means,'Filled','k','SizeData',50);

    % Aesthetics
    box off
    xlim([2,13]);
    ax=gca;
    ax.FontSize=17;
    set(gca,'FontName','Latin Modern Roman')
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);

    xlabel(tl,'Time after infection, hr','FontSize',18,'FontName','Latin Modern Roman');
    ylabel(tl,'Proportion of lysed to infected cells','FontSize',18,'FontName','Latin Modern Roman')
end

% Add legend 
nexttile(2);
legend('CI','Best fit','Replicates mean','Replicates','Direction','reverse');
legend('Location','northeastoutside');
legend('boxoff');

%%
saveas(gcf,'../../Figures/SupplFig3_CDFmodels.svg');