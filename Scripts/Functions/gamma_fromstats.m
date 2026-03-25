function [shape,scale] = gamma_fromstats(mean,cv)
    % Calculate gamma shape and scale given a mean and 
    % coefficient of variation
    shape = 1 ./ cv .^ 2;
    scale = (cv .^ 2) ./ (1 ./ mean);
end