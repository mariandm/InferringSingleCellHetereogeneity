% Code by Marian Dominguez-Mirazo, 2024
% Repurposed from Dominguez-Mirazo et al, mBio, 2024
% Relationship between population-level model number of E compartments and
% latent period distribution(Erlang distribution) CV
% Output: ~/Figures/SupplFig1_EtoCV.svg
%%
clear all; close all; clc;
%% Plot E to CV relation
ns = 0:99;
cvs = 1./sqrt(ns+1);
figure('Position',[10,10,400,300]);
plot(ns,cvs,'LineWidth',2,'Color','k');

%% Aesthetics
xlabel('Number of E compartments','FontSize',18);
ylabel('Coefficient of Variation','FontSize',18);
xticks([0,25:25:100]);
box off;
ax=gca;
ax.FontSize=15;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.03,0.03]);
set(gca,'LineWidth',0.6);

%%
saveas(gcf,'../../Figures/SupplFig1_EtoCV.svg');