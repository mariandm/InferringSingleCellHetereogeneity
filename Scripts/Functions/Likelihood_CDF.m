function ll = Likelihood_CDF(exp,x0,plateau,model)
    % Calculate likelihood for single-cell experimental observations
    % Given a latent period distribution model
    % exp: single-cell experimental counts
    % x0: latent period distribution parameters (T,CV,plateau)
    % plateau: incorporate plateau correction 
    % model: gamma or lognormal 
    % Code by Marian Dominguez-Mirazo, 2025

    % Given specific CDF model and parameters, get the lysis probability 
    % at the experimental timepoints stored in exp(:,1)
    ps = get_CDF(exp(:,1),x0,plateau,model);

    % Total infected cells
    trials1 = exp(:,4);
    % reorder
    trials = repmat(trials1,1,size(ps,2));
    % More than 1 plaque (lysed)
    ks1 = exp(:,5);
    % reorder
    ks = repmat(ks1,1,size(ps,2));
    % Binomial probability
    binoprobs = binopdf(ks,trials,ps);
    % Make probabilities of 0 very small to avoid getting stuck on Inf
    binoprobs(binoprobs==0) = 1e-300;
    % Make negative for minimum finding function
    ll = -sum(log(binoprobs));
end
