function [burst_distr] = burst_mm(pars,xs)
    % Latent period (xs) to burst size relationship
    % Assuming Hill function
    % Code by Marian Dominguez-Mirazo 2025

    % Set parameters
    pars_d = pars(1);
    pars_kmax = pars(2);
    pars_LT50 = pars(3);

    % Calculate discrete burst distribution
    burst_distr = pars_kmax .* (xs-pars_d) ./ ...
        (pars_LT50-pars_d + xs-pars_d);
    % Make 0 anything before the eclipse time
    burst_distr(xs<pars_d)=0;
end