function [CIdist_e,CIdist_pe] = splitDistr_eclipse(CIdist,var_p,emean)
    % Given a Latent Period distribution (CIdist), a proportion of eclipse
    % period variability that explains latent period variability (var_p),
    % and a mean eclipse period, calculate the mean and CV of the eclipse
    % period and post-eclipse period distributions
    % Code by Marian Dominguez-Mirazo, 2026
    
    % Latent period information
    lp_cv = CIdist(2,1);
    lp_mean = CIdist(1,1);
    lp_var = (lp_cv * lp_mean)^2;

    % Eclipse dist
    evar = lp_var * var_p;
    ecv = sqrt(evar)/emean;
    CIdist_e = [emean,1;ecv,1];

    % Post-eclipse dist
    pevar = lp_var * (1-var_p); % post-eclipse variance
    pemean = lp_mean - emean;
    pecv =  sqrt(pevar)/pemean;
    CIdist_pe = [pemean,1;pecv,1];

    disp(evar);
    disp(pevar);
end