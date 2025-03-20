% Code by Marian Dominguez-Mirazo, 2025
% This script plots the number of simulated single-cell experiments for
% which the latent period mean and CV was succesfully predicted. 
% The simulated experiments assume a range of underlying distributions.
% Required pre-computed file: ~/IntermediateFiles/CDFAccuracy.mat
% Output: ~/Figures/SupplFig6_CDFAccuracy.svg
%%
clear all; close all; clc;
addpath('../../IntermediateFiles/')
%% Load pre-run data
load('CDFAccuracy.mat');

%% Plot

% Flip matrix for plotting
result2plot=flipud(transpose(result(1:9,:)));

% Plot heatmap
figure('Position',[10,10,550,500]);
imagesc(result2plot);
caxis([0,100]);
colorbar;
xticklabels(mainTs(1:9));
yticklabels(flip(maincvs));
xlabel('Latent Period, hr');
ylabel('CV');
title('Number of experiments successfully predicted');

% Aesthetics
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman')
set(gca,'LineWidth',0.6);

%%
saveas(gcf,'../../Figures/SupplFig6_CDFAccuracy.svg')