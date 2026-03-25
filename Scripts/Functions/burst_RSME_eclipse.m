function [RSME] = burst_RSME_eclipse(data,sample_time,CIdist,var_p,pars,model)
    % Calculate the root mean squared error (RMSE) given burst size
    % experimental data and a set of parameters assuming eclipse period
    % Code by: Marian Dominguez-Mirazo, 2025
    
    thispars = [0,pars(2:end)];
    
    % Calculate the distribution for the eclipse and post-eclipse
    % periods

    % Prevent eclipse period average to go below 0.1 hr 
    % (to avoid negative values)
    if pars(1)<0.1
        emean = 0.1;
    else
        emean = pars(1);
    end
    [CIdist_e,CIdist_pe] = splitDistr_eclipse(CIdist,var_p,emean);

    % Calculate beta effective for the parameter set
    beta_eff = calculate_betaeff_eclipse(sample_time,CIdist_pe,CIdist_e,thispars,model);
    
    % Get RSME
    RSME = sqrt(nanmean((data - beta_eff').^2,[1,2]));
end