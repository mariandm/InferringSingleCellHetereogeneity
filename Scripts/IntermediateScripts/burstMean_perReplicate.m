% Code by: Marian Dominguez-Mirazo
% Calculate the cumulative effective burst size from original single cell
% data and store in intermediate file
% Output: '~/IntermediateFiles/burstMean_perReplicate.csv'
%%
clc; clear all; close all;
addpath('../../Data/');
addpath('../Functions/');

%% Get effective burst size per replicate
% Load file
file = 'SingleCellData.csv';
tab = table2array(readtable(file, 'ReadVariableNames', true));
% Get sampling times (they are the file header)
read_times = table2array(readtable(file));
sample_time = read_times(1,3:end);

%% Calculate the cumulative effective burst size per replicate
% Get replicate ids
replicates = unique(tab(:,1));
% Initialize storage
efburst_per_replicate = nan(numel(replicates)+1,numel(sample_time));
% Include sample time in first row
efburst_per_replicate(1,:) = sample_time;

% Get mean per time per replicate
for i = 1:numel(replicates)
    % Get specific replicate data
    this_burst = tab(tab(:,1)==replicates(i),3:end);
    % Remove data with one plaque or less
    this_burst(this_burst<2)=NaN;
    % Calculate mean
    efburst_per_replicate(i+1,:) = nanmean(this_burst,1);
end
%%
csvwrite('../../IntermediateFiles/burstMean_perReplicate.csv',efburst_per_replicate);
