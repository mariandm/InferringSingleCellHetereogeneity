function prob = pdf_ads(t,x0,params)
    % Code by Marian Dominguez-Mirazo, 2025
    % Compute an PDF for adsorption times given initial cell and phage
    % densities, adsorption rate, and co-incubation times
    
    S0 = x0(1); % Initial cell density
    V0 = x0(end); % Initial phage density

    phi = params.phi; % adsorption rate
    tf = params.tf; % incubation time

    prob = (phi .* S0 .* V0 .* (S0 - V0) * exp(phi .* t .* (S0 - V0))) ./ ...
        ((V0 - S0 * exp(t .* phi .* (S0 - V0))) .^ 2 .* ...
        (1 / (1 - V0/S0) - 1 / (1 - (V0 / S0) * exp(-phi * (S0 - V0) * tf))));
end