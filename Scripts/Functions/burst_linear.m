function [burst_distr] = burst_linear(pars,xs)
    % Latent period (xs) to burst size relationship
    % Assuming linear function
    % Code by Marian Dominguez-Mirazo 2025

    % Set parameters
    pars_d = pars(1);
    pars_r = pars(2);

    % Calculate discrete burst distribution
    burst_distr = pars_r*(xs-pars_d);
    % Make effective burst 0 for times before the eclipse period
    burst_distr(find(xs < pars_d)) = 0;
end