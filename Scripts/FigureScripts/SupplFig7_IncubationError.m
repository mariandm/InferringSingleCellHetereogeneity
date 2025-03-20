% Code by Marian Dominguez-Mirazo, 2025
% Caution: This script may take some time to run
% Plot simulates observed means and CVs for different incubation times in 
% the single-cell lysis detection protocol
% Required pre-computed file: 
% - ~/IntermediateFiles/MLE_CDF.mat
% - ~/IntermediateFiles/incubationError.mat
% Output: ~/Figures/Suppl_IncubationError.svg

%%
clear all; close all; clc;
addpath('../Functions');
addpath('../../IntermediateFiles/');

%% Load predicted latent period distribution
load('MLE_CDF.mat');
gamma_CI = squeeze(CI(2,:,:));

%% Load observed latent period mean and CV for varying incubation times
load('incubationError.mat')

%% Plot
figure('Position',[10,10,800,400]);
tiledlayout(1,2);

% Observed average latent period for varying incubation times
nexttile(1);
plot(adsTimes,means,'LineWidth',1.5,'Color','k'); hold on;
yline(gamma_CI(1,1),'LineWidth',1.5,'LineStyle','--','Color','#808080');
scatter([15,15],[means(adsTimes==15),gamma_CI(1,1)],'filled','MarkerFaceColor','r');
ylim([7,8]);
ylabel('Observed latent period mean, hr');

% Observed CV latent period for varying incubation times
nexttile(2);
plot(adsTimes,cvs,'LineWidth',1.5,'Color','k'); hold on;
yline(gamma_CI(2,1),'LineWidth',1.5,'LineStyle','--','Color','#808080');
scatter([15,15],[cvs(adsTimes==15),gamma_CI(2,1)],'filled','MarkerFaceColor','r');
ylim([0.15,0.3]);
ylabel('Observed latent period CV');
legend({'Observed','Underlying','Experimental incubation'});

% Aesthetics
for i = 1:2
    nexttile(i);
    ax=gca;
    ax.FontSize=17;
    xlim([5,30]);
    set(gca,'FontName','Latin Modern Roman')
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    box off;
    xlabel('Incubation time, hr'); 
end
%%
saveas(gcf,'../../Figures/SupplFig7_IncubationError.svg');