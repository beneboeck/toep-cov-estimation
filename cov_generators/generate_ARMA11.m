function Cov = generate_ARMA11(N, sigma, a, b)
    % Inputs:
    % N     - Size of the covariance matrix
    % sigma - Standard deviation of the white noise
    % a     - AR(1) coefficient
    % b     - MA(1) coefficient
    %
    % Output:
    % Cov   - Generated covariance matrix (N x N)

    % Precompute constants
    a2 = a^2;
    a_b_sum = b + a;
    factor = sigma^2 / (1 - a2);
    
    % Initialize the autocovariance vector r
    r = zeros(N, 1);
    r(1) = sigma^2 * (1 + (a_b_sum^2) / (1 - a2));
    r(2) = factor * (a_b_sum * (1 + a));
    
    % Calculate the rest of the autocovariance values using the AR(1) property
    for k = 3:N
        r(k) = a * r(k-1);
    end
    
    % Generate the covariance matrix
    Cov = toeplitz(r);
    
    % Check if the covariance matrix is positive definite
    if all(eig(Cov) > 0)
        disp('stable')
    end
end