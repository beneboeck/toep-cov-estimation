function K = ExpBounds(eta, P, w, lam)
    % Initialize K with zeros
    K = zeros(P-1, 1);
    
    % Calculate the exponential bounds for the first w elements
    indices = 1:w;
    K(indices) = eta * exp(-lam * indices);
end