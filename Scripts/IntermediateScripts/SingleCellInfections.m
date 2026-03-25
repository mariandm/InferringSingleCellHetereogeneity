% Code by Marian Dominguez-Mirazo, 2025
% This script takes experimental single cell burst size data 
% from the 'single-cell lysis detection protocol'
% and transforms it to infected cell and lysed cell counts
% Output: ~/IntermediateFiles/SingleCellInfections.csv

%%
clc; clear all; close all;
addpath('../../Data/');

%% Get infection and lysis information from burst size count

file = 'SingleCellData.csv';
tab = readtable(file, 'ReadVariableNames', true);
tab = table2array(tab);
% Get sampling times (they are in the file header)
read_times = table2array(readtable(file));
sample_time = read_times(1,3:end);
% Get replicate ids
replicates = unique(tab(:,1));

exp = [];
% Loop per replicate
for i = 1:numel(replicates)
    this_replicate = tab(tab(:,1)==replicates(i),3:end);
    % Number of succesful plates (sometimes there are NA)
    this_replicate_succs = sum(~isnan(this_replicate),1);
    % index for which we have no plates
    idx_nan = find(this_replicate_succs==0);
    % Number of plates with succesful infections (one plaque or more)
    this_replicate_infected = nansum(this_replicate>0,1);
    % Number of plates with lysed cells (more than one plaque)
    this_replicate_lysed = nansum(this_replicate>1,1);

    % Turn into NaN all timepoints with no data
    this_replicate_succs(idx_nan)=NaN;
    this_replicate_infected(idx_nan)=NaN;
    this_replicate_lysed(idx_nan)=NaN;

    % Store this replicate
    replicate_repetition = repmat(replicates(i),size(sample_time,2),1);
    save_this_replicate = [sample_time',replicate_repetition,...
        this_replicate_succs',...
        this_replicate_infected',...
        this_replicate_lysed'];

    % Store
    exp = [exp; save_this_replicate];
end

%%
writematrix(exp,'../../IntermediateFiles/SingleCellInfections.csv')

