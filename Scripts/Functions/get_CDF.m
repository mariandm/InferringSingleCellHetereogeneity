function ps = get_CDF(xs,params,plateau,model)
    % Get CDF probabilities for the specified model and parameters
    % Code by Marian Dominguez-Mirazo, 2025

    T = params(:,1)'; % Distribution mean
    cv = params(:,2)'; % Distribution CV

    % Model with plateau correction
    if plateau
        pl = params(:,3)';
    else 
        pl = ones(1,numel(T));%1;
    end

    % Get CDF probabilities
    if model == "gamma"
        % Calculate gamma shape and rate
        [this_shape, this_scale] = gamma_fromstats(T,cv);
        % Get probability and multiply by plateau factor
        ps = zeros(numel(xs),numel(this_shape));
        for i = 1:numel(this_shape)
            this_shape2 = this_shape(i);
            this_scale2 = this_scale(i);
            ps(:,i) = gamcdf(xs,this_shape2,this_scale2) .* pl(i);
        end

    elseif model == "lognormal"
        % Calculate logn mu and sigma
        [this_mu, this_sigma] = lognormal_fromstats(T,cv);
        % Get probability and multiply by plateau factor
        ps = logncdf(xs,this_mu,this_sigma) .* pl;
    end

end
