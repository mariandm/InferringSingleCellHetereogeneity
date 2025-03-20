% Code by Marian Dominguez-Mirazo, 2025
% This script generates Figure 1C,D
% It requires the precomputed file: ~/IntermediateFiles/visualCI_popLevel.mat
% Output: ~/Figures/Figure1CD_PopLevel_Fit_PDF.svg

%% Plot population-level fit
clc; clear all; close all;
addpath('../../IntermediateFiles/');
addpath('../Functions/');
addpath('../../Data/');
%% Load data
% Viral timeseries
file = 'TimeseriesVirus.csv';
virus = readtable(file);
virus = virus{:,:};

%% Set host-only parameters
% Obtained from host parameter MCMC
pars.mu     = 0.04;    % Host growth, hr^-1
pars.K      = 3e8;    % Carrying capacity, CFU/ml

%% Load pre-computed CI
load('../../IntermediateFiles/visualCI_popLevel.mat');

%% Read viral chain

% Viral chain
file = 'MCMCChain_020325.csv';
tab = readtable(file, 'ReadVariableNames', false);
tab = table2array(tab);
H0_chain = tab(:,1);
V0_chain = tab(:,2);
phi_chain = tab(:,3);
beta_chain = tab(:,4);
eta_chain = tab(:,5);
cv_chain = tab(:,6);

%% Run the best fit with viruses
options = odeset('AbsTol',1e-6,'RelTol',1e-6);

pars.phi    = 10^mean(phi_chain);   % Adsorption rate, ml/(CFUxhr)
pars.beta   = mean(beta_chain);    % Burst size
pars.eta = mean(eta_chain); % Lysis rate
pars.cv = mean(cv_chain);
pars.n = round(1/pars.cv^2 -1);
pars.initS  = mean(H0_chain); % Initial host density, CFU/ml
pars.initV  = mean(V0_chain); % Initial viral density, PFU/ml

x0 = zeros(pars.n+3,1);
x0(1) = pars.initS; x0(end) = pars.initV;
t = 0:0.1:virus(end,1);
[tsol,ysol] = ode45(@ODE_SEnIV_sink,t,x0,options,pars);

%% Start figure
figure('Position',[10,10,1100,450]);
tl = tiledlayout(1,2);

%% Plot dynamics
nexttile(1);
% CI
fill([t';flip(t)'],[V_interval(:,1);flip(V_interval(:,2))],...
    'k','FaceColor','r','EdgeColor','None','FaceAlpha',0.2); hold on;
% Best fit 
plot(tsol,ysol(:,end),'r','LineWidth',1.5);
% Experimental mean
scatter(virus(:,1),nanmean(virus(:,2:end),2),'filled','MarkerFaceColor','k',...
    'SizeData',50); 
% Experimental Replicates
for i=2:5
    scatter(virus(:,1),virus(:,i),'filled','MarkerFaceColor','#999999',...
        'SizeData',50);
end
% Experimental mean
scatter(virus(:,1),nanmean(virus(:,2:end),2),'filled','MarkerFaceColor','k',...
    'SizeData',50); 
set(gca, 'YScale', 'log');
legend('CI','Best fit','Replicates mean','Replicates');
legend('Direction','reverse');
legend('Location','SouthEast');
legend('boxoff');

% Aesthetics
ax=gca;
ax.FontSize=17;
set(gca,'FontName','Latin Modern Roman')
set(gca,'TickDir','out');
set(gca,'TickLength',[0.025,0.025]);
set(gca,'LineWidth',0.6);
yticks(0:0.1:0.8);
box off;
xlabel('Time, hr','FontSize',18);
ylabel('PFU/ml','FontSize',18);
ylim([5e4,1e10]);
yticks(10.^(5:10));
xlim([0,24.5]);

%% Plot PDF
nexttile(2);

xs = 0:0.1:14;
% Save lysis rate chain
eta_chain = tab(:,5);
% Save CV chain
cv_chain = tab(:,6);
% Get point estimates
point_eta = mean(eta_chain);
point_cv = mean(cv_chain);

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
 'FaceColor','r','EdgeColor','None','FaceAlpha',0.1);


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

hL = legend('Pop-level best distribution', 'Pop-level best mean', 'Pop-level CI');
legend('boxoff');
legend('Location','northeast')

%%
saveas(gcf,'../../Figures/Figure1CD_PopLevel_Fit_PDF.svg');