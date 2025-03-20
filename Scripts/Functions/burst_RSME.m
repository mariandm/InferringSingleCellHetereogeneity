function [RSME] = burst_RSME(data,sample_time,CIdist,pars,model)
    % Calculate the root squared mean error (RSME) given burst size
    % experimental data and a set of parameters
    % Code by: Marian Dominguez-Mirazo, 2025

    % Calculate beta effective for the parameter set
    beta_eff = calculate_betaeff(sample_time,CIdist,pars,model);
    % Get RSME
    RSME = sqrt(nanmean((data - beta_eff').^2,[1,2]));
end