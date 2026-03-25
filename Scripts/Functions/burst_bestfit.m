function [best_fit] = burst_bestfit(data,sample_time,CI_LP,model,nrand)
    % Randomly samples a parameter space to find the global 
    % minimum best patameter combination that describes the effective burst
    % size data 
    % Code by: Marian Dominguez-Mirazo, 2025

    % Set parameter search range  
    d_range = [0,7]; % Time of first burst (hr)
    kmax_range = [0,500]; % Maximum burst size  
    LT50_range = [0,20]; % Inflexion time or half saturation time (hr)
    
    % Set model-dependent parameter ranges
    % and sample parameter ranges
    if model=="linear"
        nparam = 2;
        r_range = [0,100];
        randsamp = lhsu([d_range(1),r_range(1)],[d_range(2),r_range(2)],nrand);
    
    elseif model=="kannoly"
        nparam = 4;
        r_range = [0,5];
        randsamp = lhsu([d_range(1),r_range(1),kmax_range(1),LT50_range(1)],...
            [d_range(2),r_range(2),kmax_range(2),LT50_range(2)],nrand);

    elseif model=="mm"
        nparam = 3; 
        randsamp = lhsu([d_range(1),kmax_range(1),LT50_range(1)],...
            [d_range(2),kmax_range(2),LT50_range(2)],nrand);
    end

    % store
    s = zeros(nrand,nparam+1);
    
    % For each initial condition look for minimum and store
    for i = 1:nrand
        % initial condition for fminsearch
        x0 = randsamp(i,:);
        % Search function
        [this_s, this_fval,exit_function] = fminsearch(@(pars)burst_RSME(data,sample_time,CI_LP,pars,model),x0);
        
        % Save output only if the search converged
        if exit_function==1
            s(i,1:nparam) = this_s;
            s(i,end) = this_fval;
        else
            s(i,:) = NaN;
        end
    end

    % Remove parameters below 0 (not biologically feasible)
    s2 = s(logical(prod(s>1,2)),:);

    % Find best fit
    [A,I] = min(s2(:,end));
    best_fit = s2(I,:);
end