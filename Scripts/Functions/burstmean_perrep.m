function [efburst_per_replicate] = burstmean_perrep(data_table)
    % Given experimental data from the 'single-cell lysis detection
    % protocol', calculate average burst size per
    % replicate across timepoints 
    % Code by: Marian Dominguez-Mirazo, 2025
    
    % Read data
    tab = table2array(data_table);
    % Get replicate ids
    replicates = unique(tab(:,1));
    
    % Get mean per time per replicate
    efburst_per_replicate = [];
    for i = 1:numel(replicates)
        % Get specific replicate data
        this_burst = tab(tab(:,1)==replicates(i),3:end);
        % Remove data with less than one plaque
        this_burst(this_burst<2)=NaN;
        % Calculate mean
        efburst_per_replicate(i,:) = nanmean(this_burst,1);
    end
end