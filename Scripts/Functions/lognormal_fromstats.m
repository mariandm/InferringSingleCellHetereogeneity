function [mu,sigma] = lognormal_fromstats(mean,cv)
    
    % Calculate lognormal mu and sigma given a mean and 
    % coefficient of variation

    sd = cv .* mean;
    variance = sd .^ 2;
    mu = log(mean .^2 ./ sqrt(variance + mean .^ 2));
    sigma = sqrt(log(variance ./ mean .^2 +1));
end