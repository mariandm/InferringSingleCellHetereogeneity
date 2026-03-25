function [beta_eff2] = calculate_betaeff_eclipse(sample_time,CIdist_pe, CIdist_e,pars,model)
    % Calculate expected beta effective given a latent period to burst size
    % model and corresponding parameters, and a latent period distribution,
    % accounting for eclipse period variability.
    % The code assumes exponential decay of viral particles due to
    % attachment of particles to the well plate plastic surface at a rate
    % of -0.04. See supplementary material.
    % Code by: Marian Dominguez-Mirazo, 2026

    % Calculate gamma shape and scale
    % Post-eclipse period distribution
    [shape_single,scale_single] = gamma_fromstats(CIdist_pe(1,1),CIdist_pe(2,1));
    % Eclipse period distribution
    [shape_single_e,scale_single_e] = gamma_fromstats(CIdist_e(1,1),CIdist_e(2,1));

    % Initialize data storage
    xs = linspace(0,max(sample_time)+1,1e4);
    
    % Discretize a pdf by substracting cdfs
    this_cdf = gamcdf(xs,shape_single,scale_single);
    this_pdf = this_cdf(2:end) - this_cdf(1:end-1);
    this_cdf_e = gamcdf(xs,shape_single_e,scale_single_e);
    this_pdf_e = this_cdf_e(2:end) - this_cdf_e(1:end-1);
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
    partial_cdf = cumsum(this_pdf);
    partial_cdf_e = cumsum(this_pdf_e);

    % Ahora se viene lo hardcore
    % Iterando sobre epsilons
    beta_eff2=zeros(numel(sample_time),1);
    % loop over the sample time
    for i = 1:numel(sample_time)
        % Last index where the sampling time is larger than the xs vector
        maxj = (sum(new_xs<sample_time(i))-1);
        % First index where the probability of e is larger than 0
        tmp = find(partial_cdf_e>0);
        minj = tmp(1);
        % Initialize
        beta_eff1=zeros(maxj,1);
        % Loop over eclipse time
        for j = minj:maxj
            % Find the relevant index
            % pe + e > t
            idx = find((new_xs+new_xs(j))>sample_time(i),1)-1;
            % virion decay 
            ed = exp((-0.04) .* (sample_time(i) - new_xs(1:idx) - new_xs(j)));
            bs_ed = cumsum(bs(1:idx) .* ed);
            
            beta_eff1(j) = bs_ed(idx) / partial_cdf(idx) * this_pdf_e(j) / partial_cdf_e(maxj);
        end
        beta_eff2(i) = sum(beta_eff1);
    end
end