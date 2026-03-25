function [beta_eff] = calculate_betaeff(sample_time,CIdist,pars,model)
    % Calculate expected beta effective given a latent period to burst size
    % model and corresponding parameters, and a latent period distribution
    % The code assumes exponential decay of viral particles due to
    % attachment of particles to the well plate plastic surface at a rate
    % of -0.04. See supplementary material.
    % Code by: Marian Dominguez-Mirazo, 2025

    % Calculate gamma shape and scale
    [shape_single,scale_single] = gamma_fromstats(CIdist(1,1),CIdist(2,1));

    % Initialize data storage
    beta_eff = zeros(numel(sample_time),1);
    xs = linspace(0,max(sample_time)+1,1e5);
    
    % Discretize a pdf by substracting cdfs
    this_cdf = gamcdf(xs,shape_single,scale_single);
    this_pdf = this_cdf(2:end) - this_cdf(1:end-1);
    new_xs = (xs(1:end-1) + xs(2:end)) ./ 2;

    % Get the burst distribution for discrete steps
    if model=="linear"
        burst_distr = burst_linear(pars,new_xs);
    elseif model=="kannoly"
        burst_distr = burst_kannoly(pars,new_xs);
    elseif model=="mm"
        burst_distr = burst_mm(pars,new_xs);
    end

    % Multiply burst size by probability
    bs = burst_distr .* this_pdf;

    for i = 1:numel(sample_time)
        % Find probabilities before sample timepoint
        idx = find(new_xs>sample_time(i),1)-1;

        % Include burst size decay
        ed = exp((-0.04) .* (sample_time(i) - new_xs(1:idx)));
        bs_ed = cumsum(ed .* bs(1:idx)); 
        bsp_ed = bs_ed(idx);

        % Normalize to lysis probability by sampling time
        partial_cdf = gamcdf(new_xs(idx),shape_single,scale_single);
        % and store
        beta_eff(i) = bsp_ed/ partial_cdf;
    end
end