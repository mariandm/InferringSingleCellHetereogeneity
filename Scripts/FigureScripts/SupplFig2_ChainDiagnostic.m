% Code by Marian Dominguez-Mirazo, 2025
% Adapted from Dominguez-Mirazo et al., 2024, mBio
clear all; close all; clc;
addpath('../../IntermediateFiles/');
%%
diagnostic = csvread('MCMCdiagnostic.csv');
param_names = ["S_0","V_0","\phi","\beta","\eta","CV"];

%%
figure('Position',[10,10,800,450]);

my_yticks = size(diagnostic,1):-1:1;

nexttile(1);
xline(1,'LineStyle','--','LineWidth',1.5);hold on;
scatter(diagnostic(:,1),my_yticks,'filled','k'); 
xlim([0.95,1.05]);
xlabel('$\hat{R}$','Interpreter','latex');


nexttile(2);
xline(0.1,'LineStyle','--','LineWidth',1.5);hold on;
scatter(diagnostic(:,2),my_yticks,'filled','k'); 
xlim([0,0.2]);
xlabel('ESS');


for i = 1:2
    nexttile(i);
    yticks(flip(my_yticks));
    yticklabels(flip(param_names));
    ylabel('Parameter');
    
    set(gca,'FontName','Latin Modern Roman');
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.025,0.025]);
    set(gca,'LineWidth',0.6);
    ax=gca;
    ax.FontSize=17;
    ylabel('Parameter');
end
%%
saveas(gcf,'../../Figures/SupplFig2_ChainDiagnostic.svg');