% Code by Marian Dominguez-Mirazo, 2025
% Caution: this script takes multiple hours to run
% For a simplified version change to nexp = 5;
% This script simulates single cell experiments for a range of Ts and CVs
% and then predicts the latent perdiod distribution MLE, as we did for the 
% real experimental data. It stores the number of simulate experiments for 
% which the paramaters were succesfully predicted
% Output: ~/IntermediateFiles/CDFAccuracy.mat
%%
tic
clear all; close all; clc;
addpath('../Functions/') 

%% Information on the experimental details
expInfo.sample_size = 30; % Number of plates per timepoint
expInfo.sample_start = 4; % hr, Time of experiment start
expInfo.sample_frequency = 0.5; % h^-1, Sampling frequency
expInfo.sample_end = 12; % hr, Time of experiment end
expInfo.nrep = 4; % Number of replicates
% Parameter info
pars.infected_prob = 0.4; % Probability of cell being infected 

%% Ts and CVs for with which we will simulate the experiments
% i.e. the underlying true parameters
mainTs = 5:0.5:10; % range of T (mean LP) to explore
maincvs = 0.05:0.05:0.4; % range of CV to explore

% Number of experiments to simulate
nexp = 5;

%% Set up for MLE 
nrand = 1e1; % Number of initial conditions for fminsearch 
% x0 range for distribution mean
Trange = [4,10]; 
% x0 range for CV
CVrange = [0.05,0.4];
% x0 range for pl parameter
PLrange = [0.4,1];

% Confidence Interval search ranges
Ts = Trange(1):0.1:Trange(2);
cvs = CVrange(1):0.01:CVrange(2);
pls = PLrange(1):0.01:PLrange(2);
% Combination of parameters
x_combns_pl = table2array(combinations(Ts,cvs,pls));

% Initialize storage
% Stores the number of experiments for which we succesfully predicted T
% within 1 hr of the underlying truth, and CV within a 0.1 distance
result = zeros(numel(mainTs),numel(maincvs));

%%
% How many succesful infections will result in a single plaque
pars.pl = 0.85;

% Simulated T 
for thisT = 1:numel(mainTs)
    disp(thisT);
    pars.T = mainTs(thisT);
    % Simulated CV
    for thisCV = 1:numel(maincvs)
        pars.cv = maincvs(thisCV);

        % Number of experiments to simulate
        for thisexp = 1:nexp

            % Simulate experiment
            exp = simulateSingleCellExp(pars,expInfo);
    
            % Predict distribution
            % Store data
            s_gamma_pl = zeros(nrand,4);

            % Get an x0 sample
            randsamp = lhsu([Trange(1),CVrange(1),PLrange(1)], ...
                [Trange(2),CVrange(2),PLrange(2)], ...
                nrand);

            % For each of the initial conditions for fminsearch
            for i = 1:nrand
                x0 = randsamp(i,:);
                % Gamma with plateau 
                % Search function 
                [this_s_gammapl, this_fval_gammapl] = fminsearch(@(x)Likelihood_CDF(exp,x,1,"gamma"),x0);
                % Save output
                s_gamma_pl(i,1:3) = this_s_gammapl;
                s_gamma_pl(i,4) = this_fval_gammapl;
            end
            % Ensure positive values for predicted T and cv
            idx = logical(prod(s_gamma_pl>0,2));

            % Get MLE
            [A_gamma_pl,I_gamma_pl] = min(s_gamma_pl(:,end));
            MLE_gamma_pl = s_gamma_pl(I_gamma_pl,1:3);

            % Call success or failure of prediction
            if pars.T >= MLE_gamma_pl(1)-1 && pars.T <= MLE_gamma_pl(1)+1 ...
                    && pars.cv >= MLE_gamma_pl(2)-0.1 && pars.cv <= MLE_gamma_pl(2)+0.1
                result(thisT,thisCV) = result(thisT,thisCV) + 1;
            end
        end
    end
end

%%
save('../../IntermediateFiles/CDFAccuracy.mat','result','mainTs','maincvs');

toc