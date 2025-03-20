function experiment = simulateSingleCellExp(params,expInfo)
    % Code by Marian Dominguez-Mirazo, 2025
    % Simulate the single-cell lysis detection protocol by:
    % Getting a random number of infected cells
    % Getting a random latent period per individual cell
    % Finding whether the cell lysed by the sampling time point

    % How many plates per time point
    sample_size = expInfo.sample_size;
    % Time points (in hours)
    sample_points = expInfo.sample_start:expInfo.sample_frequency:expInfo.sample_end;
    % Number of replicates
    nrep = expInfo.nrep;
    % Total number of plates
    nplates = sample_size*numel(sample_points)*nrep;

    % Get random latent periods
    % Assuming a Gamma distribution
    [thisshape,thisscale] = gamma_fromstats(params.T,params.cv);
    random_LP = gamrnd(thisshape,thisscale,sample_size,numel(sample_points),nrep);

    % Get the ones that lysed by sampling point
    lysed = random_LP <= sample_points;

    % Get random number of infected cells 
    random_ninfected = binornd(sample_size,params.infected_prob,numel(sample_points),nrep);

    % initialize storage
    experiment = nan(numel(random_ninfected),5);
    % Get proportion of lysed cells per replicate and time
    cnt = 1;
    for i = 1:size(random_ninfected,1)
        for j = 1:size(random_ninfected,2)
            ninf = random_ninfected(i,j);
            nlysed = sum(lysed(1:ninf,i,j));
            % Simulate a percentage of cells having a burst size of 1
            nlysed2 = binornd(nlysed,params.pl);
            % Save in order: timepoint, replicate, nplates, ninfected,
            % nlysed
            experiment(cnt,:) = [sample_points(i),j,sample_size,ninf,nlysed2];
            cnt = cnt +1;
        end
    end
end