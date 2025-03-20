% Code by Marian Dominguez-Mirazo, 2025
% Plot experiment data on phage adhesion to plastic wall of 96-well plate
% output: Figures/Suppl_adhesion.svg
clc; clear all; close all;
%% Load data
tab = readmatrix('../../Data/PlasticAdhesion.csv');
time = tab(1,2:end);
counts = tab(2:end,2:end);
%% Arrange data 
counts_linear = reshape(counts,[numel(counts),1]);
time_linear = sort(repmat(time,1,size(counts,2)));

%% Find best linear fit for exponential decay
lm = fitlm(time_linear,log(counts_linear));

%% Compute fit line for plot
xs = time(1):0.1:time(end);
ys = exp(lm.Coefficients.Estimate(1)) .* exp(lm.Coefficients.Estimate(2) .*xs);

%% Plot
figure('Position',[10,10,400,400]);
plot(xs,ys,'LineWidth',1.5); hold on;
scatter(time,mean(counts,1),'filled','MarkerFaceColor','k');
scatter(time,counts','filled','MarkerFaceColor','#999999'); 
scatter(time,mean(counts,1),'filled','MarkerFaceColor','k');

legend('Exponential fit: $80.4e^{-0.04\,t}$',...
    'Replicate mean', ...
    'Replicates',...
    'Interpreter','latex');

ylabel('PFU per well');
xlabel('Time in well, hr');

box off;
ax=gca;
ax.FontSize=18;
set(gca,'FontName','Latin Modern Roman');
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);

%%
saveas(gcf,'../../Figures/SupplFig8_adhesion.svg')