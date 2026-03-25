function [beta] = LP_burstdistr_eclipse(xs,CIdist_e,pars,model)
    % Calculate the latent period to burst size relationship for a model
    % that integrates eclipse period variability by convoluting the eclipse
    % period probability distribution and the production time to burst size
    % function
    % Code by marian Dominguez-Mirazo, 2026

    beta = zeros(numel(xs),1);
    % Calculate values for eclipse period distribution
    [shape_e,scale_e] = gamma_fromstats(CIdist_e(1,1),CIdist_e(2,1));
    dt = 0.01; % calculation resolution
    for i=1:numel(xs)
        thisxs=0:dt:xs(i); % post-eclipse time
        es = xs(i) - thisxs; % corresponding eclipse time
        fake_es = es+dt/2; % intermediate es for discrete pdf calculation
        % Probability of each eclipse time given eclipse distribution
        p_es = abs(diff(gamcdf(fake_es,shape_e,scale_e)));
        real_es = es(2:end); % es for which we have a probability
        % Calculate burst for the post-eclipse times
        if model=="linear"
            burstdistr = burst_linear(pars,thisxs);
        elseif model=="mm"
            burstdistr = burst_mm(pars,thisxs);
        elseif model=="kannoly"
            burstdistr = burst_kannoly(pars,thisxs);
        end
        % Get probability of each post-eclipse times
        beta(i) = sum(p_es./sum(p_es).*burstdistr(1:end-1));
    end
end