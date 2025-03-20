function rads = random_ads(y,x0,params)
    % Code by Marian Dominguez-Mirazo, 2025
    % Simulate a random adsorption time in single-cell lysis detection
    % protocol

    S0 = x0(1); % Initial cell density
    V0 = x0(end); % Initical phage density

    phi = params.phi; % adsorption rate
    tf = params.tf; % incubation time 

    % get y axis constraint (given tf)
    y_uplimit = pdf_ads(0,x0,params);
    y_lowlimit = pdf_ads(tf,x0,params);
    % transform random number
    y_new = y .* (y_uplimit - y_lowlimit) + y_lowlimit;

    % Like losers get the function and find the closest value for y_new
    % x range
    xt = 0:0.0001:tf*2;
    % get pdf 
    numeric_pdf = pdf_ads(xt,x0,params);

    rads = zeros(numel(y_new),1);
    % If y is many random numbers
    for i = 1:numel(y_new)
        % find the closest pdf value for y_new
        idx = find(numeric_pdf<=y_new(i),1);
        rads(i) = xt(idx);
    end
end