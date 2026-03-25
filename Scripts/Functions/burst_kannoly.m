function [burst_distr] = burst_kannoly(pars,xs)
    % Latent period (xs) to burst size relationship
    % Kannoly model from Kannoly et al. 2023, Microbiology Spectrum
    % Code by Marian Dominguez-Mirazom 2025
    
    % Set parameters
    pars_d = pars(1);
    pars_r = pars(2);
    pars_kmax = pars(3);
    pars_LT50 = pars(4);

    % Calculate discrete burst distribution
    burst_distr = pars_kmax .* (exp(pars_r .* (xs-pars_d))-1) ./ ...
        (exp(pars_r*(pars_LT50-pars_d)) + exp(pars_r .* (xs-pars_d))-2);
    % Make 0 anything before the eclipse time
    burst_distr(xs<pars_d)=0;
end