% Code by Marian Dominguez-Mirazo, 2025
% Note that this script uses the combinations function
% which requires MATLAB R2023a and above
% This script generates Figure 2C,D
% It requires the precomputed files: 
% - ~/IntermediateFiles/SingleCellInfections.csv
% - ~/IntermediateFiles/MLE_CDF.mat
% Output: ~/Figures/Figure2CD_SingleCell_CDF_PDF.svg
%% 
clc; clear all; close all;
addpath('../../IntermediateFiles/');
addpath('../Functions/');

%% Load data
% MLE and CI from single cell
load('MLE_CDF.mat');
% MCMC chain from pop-level estimates
file = 'MCMCChain_020325.csv';
tab = readtable(file, 'ReadVariableNames', false);
chain = table2array(tab);
% Experiment data
file = 'SingleCellInfections.csv';
% Column order is: (1) sample time, (2) replicate, (3) number of viable plates, (4) number
% of infected cells (plaque count > 0), (5) number of lysed cells (plaque count > 1)
exp = csvread(file);

%% Get average proportion lysed across replicates
times = unique(exp(:,1));
means = [];
for i = 1:numel(times)
    idx = exp(:,1)==times(i);
    means(i) = nanmean(exp(idx,5) ./ exp(idx,4));
end

%% Visual CDF CI
timepoints = 0:0.01:15;
% Store
maxes = zeros(numel(timepoints),1);
mines = zeros(numel(timepoints),1);
MLEs = zeros(numel(timepoints),1);

model_pls = [0,1,0,1];
model_dist = ["gamma","gamma","lognormal","lognormal"];

% The best model is the second one in the stored file
i=2;
% Get a range of parameter values
Ts = CI(i,1,2):0.1:CI(i,1,3);
cvs = CI(i,2,2):0.05:CI(i,2,3);
pls = CI(i,3,2):0.05:CI(i,3,3);
x_combns = table2array(combinations(Ts,cvs,pls));

% Compute probabilities across values
ps = get_CDF(timepoints',x_combns,model_pls(i),model_dist(i));
% store
maxes = max(ps,[],2);
mines = min(ps,[],2);

% MLE prediction
MLEs(:,i) = get_CDF(timepoints',CI(i,:,1),model_pls(i),model_dist(i));


%% Start figure
figure('Position',[10,10,1100,450]);
tl = tiledlayout(1,2);
prueba = [0.6 0.6 .08 .25];
%% Plot CDF (Figure 2C)

nexttile(1);

colors = ["b","b","#9F2B68","#9F2B68"];
markershapes = ["o","square","^","v"];

i=2;

% Plot CI (shaded area)
fill([timepoints';flipud(timepoints')],[mines;flipud(maxes)],"",...
    "FaceColor",colors(i),"EdgeColor","None","FaceAlpha",0.3); hold on;

% Best single cell prediction
plot(timepoints, MLEs(:,i), "Color",colors(i), "LineWidth",1.5);

% Plot single cell means
scatter(times,means,'Filled','k','SizeData',50);

% Plot single cell experiment data per replicate
replicates = unique(exp(:,2));
for j = 1:numel(replicates)
    idx = exp(:,2) == replicates(j);
    scatter(exp(idx,1),exp(idx,5)./exp(idx,4),'Filled',...
        'MarkerFaceColor','#999999',...
        'SizeData',50,...
        'Marker',markershapes(j));
end

% Plot single cell means
% I plot it twice so the legend appears in the order I want
scatter(times,means,'Filled','k','SizeData',50);

% Aesthetics
box off
xlim([0,12.5]);
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
xlabel('Time after infection, hr','FontSize',18);
ylabel('Proportion of lysed to infected cells','FontSize',18);

legend('CI','MLE','Replicate mean','Replicates');
legend('boxoff');
legend('Location','northwest','Direction','reverse');

%% Plot PDF (Figure 2D)
nexttile(2);
xs = 0:0.1:14;

%% Plot Population-level predicted PDF
% Save lysis rate chain
eta_chain = chain(:,5);
% Save CV chain
cv_chain = chain(:,6);
% Get point estimates
point_eta = mean(eta_chain);
point_cv = mean(cv_chain);
% Calculate gamma shape and scale
[shape_pop,scale_pop] = gamma_fromstats(1/point_eta,point_cv);

% Get 95% CI from chains
% CV
cv_95chain = cv_chain(cv_chain > quantile(cv_chain,0.025) ...
    & cv_chain < quantile(cv_chain,0.975));
% eta, lysis rate
eta_95chain = eta_chain(eta_chain > quantile(eta_chain,0.025) ...
    & eta_chain < quantile(eta_chain,0.975));

% Plot LP distribution of pop-level inference
plot(xs,gampdf(xs,shape_pop,scale_pop),...
    'Color','r','LineWidth',1.5); hold on;
% Plot underlying mean(LP)
xline(1/point_eta,'--','Color','r','LineWidth',1.5);
% Confidence interval for mean latent period
lelimit = [0,2];
fill([min(1./eta_95chain),min(1./eta_95chain),...
     max(1./eta_95chain),max(1./eta_95chain)],...
 [lelimit(1),lelimit(2),lelimit(2),lelimit(1)],'',...
 'FaceColor','r','EdgeColor','None','FaceAlpha',0.1); hold on;

%% Plot single-cell PDF
bestCI = squeeze(CI(2,:,:));
% Get gamma shape and scale
[shape_single,scale_single]= gamma_fromstats(bestCI(1,1),bestCI(2,1));

% Plot LP distribution of single cell
plot(xs,gampdf(xs,shape_single,scale_single),...
    'Color','b','LineWidth',1.5); hold on;
% Plot underlying mean(LP)
xline(bestCI(1,1),'--','Color','b','LineWidth',1.5);
% Plot 95% CI as a geometric figure
% Confidence interval for mean latent period
lelimit = [0,2];
fill([bestCI(1,2),bestCI(1,2),...
     bestCI(1,3),bestCI(1,3)],...
 [lelimit(1),lelimit(2),lelimit(2),lelimit(1)],'k',...
 'FaceColor','b','EdgeColor','None','FaceAlpha',0.1); hold on;

% Aesthetics
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
%xticks(0:xs(end));
yticks(0:0.1:0.8);
box off;
xlabel('Time, hr','FontSize',18);
ylabel('Probability Distribution','FontSize',18);

xlim([xs(1),xs(end)]);
ylim([0,0.4])

% Legend
hL = legend('Pop-level best distribution', 'Pop-level best mean', 'Pop-level CI',...
    'Single-cell best distribution', 'Single-cell best mean', 'Single-cell CI');
legend('boxoff');
legend('Location','northeast')

%% Make a tiny plot for the CV
% create smaller axes in top left, and plot on it
axes('Position',prueba)

% Single-cell CV CI
errorbar(1,bestCI(2,1),bestCI(2,1)-bestCI(2,2),bestCI(2,1)-bestCI(2,3),...
    'LineWidth',1.5,'Color','b'); hold on;
scatter(1,bestCI(2,1),'filled','b');
% Pop-level CV CI
errorbar(2,mean(cv_95chain),...
    mean(cv_95chain)-min(cv_95chain), ...
    mean(cv_95chain)-max(cv_95chain),...
    'LineWidth',1.5,'Color','r'); hold on;
% Pop-level CV best estimate
scatter(2,mean(cv_95chain),'filled','r');

% Aesthetics
xlim([0.5,2.5]);
ylim([0,0.35]);
box on;
ax=gca;
ax.FontSize=13;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
xticklabels('');
ylabel('CV');
title('CV 95% CI')

%%
saveas(gcf,'../../Figures/Figure2CD_SingleCell_CDF_PDF.svg')