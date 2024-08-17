function K = ConstantBounds(eta, P, w)
    % Initialize K with zeros
    K = zeros(P-1, 1);
    
    % Calculate the linear bounds for the first w elements
    indices = 1:w;
    K(indices) = eta;
end