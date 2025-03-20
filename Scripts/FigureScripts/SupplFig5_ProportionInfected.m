% Code by Marian Dominguez-Mirazo, 2025
% Plot proportion of infected cells across time
% Required pre-computed file: ~/IntermediateFiles/SingleCellInfections.csv
% Output: ~/Figures/SupplFig5_ProportionInfected.svg
%%
clc; clear all; close all;
addpath('../../IntermediateFiles/');

%% Load data 
file = 'SingleCellInfections.csv';
% Column order is: sample time, replicate, number of viable plates, number
% of infected cells (plaque count > 0), number of lysed cells (plaque count > 1)
exp = csvread(file);

%% Plot proportion of infected across time
replicates = unique(exp(:,2));
for i = 1:numel(replicates)
    this_exp = exp(exp(:,2)==replicates(i),:);
    scatter(this_exp(:,1),this_exp(:,4)./this_exp(:,3),'filled'); hold on;
end

% Aesthetics
ylim([0,1]);
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
xlabel('Time after infection, hr','FontSize',18);
ylabel('Proportion of infected cells','FontSize',18);
box off;
legend('Replicate 1','Replicate 2','Replicate 3', 'Replicate 4','Location','NorthWest');
legend('boxoff');

saveas(gcf,'../../Figures/SupplFig5_ProportionInfected.svg');