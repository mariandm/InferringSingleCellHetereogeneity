% Code by Marian Dominguez-Mirazo, 2026
% This script generates SupplFig11 that showcases RMSE values for a range
% of v (proportion of LP variability explained by eclipse period
% variability), and r (productiion rate), for three values of avergae
% eclipse period assuming a linear model of production period to burst size
% relationship
% It requires the pre-computed files: 
% ~/IntermediateFiles/RMSE_eclipse.mat
% Output: ~/Figures/SupplFig11_SloppyP.svg
%%
clear all; close all; clc;
addpath('../../IntermediateFiles/')
load('../../IntermediateFiles/RMSE_eclipse.mat')

%%
figure('Position',[10,10,400*3,400]);
for i = 1:3
    nexttile(i);
    RSME_matrix_forplot = flipud(transpose(squeeze(RMSE_matrix(i,:,:))));
    
    % Plot heatmap
    
    imagesc(RSME_matrix_forplot);
    % Colorbar
    h = colorbar;        % Add colorbar
    h.Title.String = 'RSME';  % Add label to colorbar
    %h.Label.Rotation = 0;     % horizontal
    %h.Label.VerticalAlignment ='middle';
    %h.Label.Position = [mean(h.Limits), 1.05, 0];

    xrange = 2:2:numel(range_varp);
    xticks(xrange);
    xticklabels(range_varp(xrange));
    yrange = 1:4:numel(range_r);
    yticks(yrange);
    yticklabels(flip(range_r(yrange)));
    
    title(strcat('\mu(\epsilon)= ',num2str(range_d(i)),'hr'));
    
    % Aesthetics
    ax=gca;
    ax.FontSize=15;
    set(gca,'FontName','Latin Modern Sans')
    set(gca,'LineWidth',0.6);
end

nexttile(1);
ylabel('r (hr^{-1}), production rate ');
nexttile(2);
xlabel('v, proportion of latent period variability explained by eclipse period');
%%
saveas(gcf,'../../Figures/SupplFig11_SloppyP.svg')
