function [CI] = ML_CI(ll_map,plateau,ML)
    % Calculate confidence intervals based on 
    % likelihood using two (no plateau) or three (with plateau) parameters
    % Code by Marian Dominguez-Mirazo, 2025
    % ll_map is a matrix containing likelihoos value in the last column, calculated 
    % for a range of parameters in columns 1 to 3
    % The plateau parameter states whether the model includes a plateau correction, 
    % and therefore a third parameter. From main text equations y=1, when pleteau =0
    % ML is the maximum likelihood estimate

    % Chi square significance value
    if plateau 
        sig_leveln = chi2inv(0.95,3);
    else
        sig_leveln = chi2inv(0.95,2);
    end

    % Get all matrix index
    idx = abs(ll_map(:,end) - ML) < sig_leveln;
    % Find range for each parameter
    CI_Ts = ll_map(idx,1);
    CI_cvs = ll_map(idx,2);
    
    % If no plateau, make it all 1
    if plateau 
        CI_pls = ll_map(idx,3);
    else
        CI_pls = 1;
    end

    % Store
    CI = zeros(2,2);
    CI(1,1) = min(CI_Ts); CI(1,2) = max(CI_Ts); 
    CI(2,1) = min(CI_cvs); CI(2,2) = max(CI_cvs); 
    CI(3,1) = min(CI_pls); CI(3,2) = max(CI_pls); 
end