function alpha = PLS(X, sCov, P, N, w, K)
    % Inputs:
    % X     - Data matrix of size (N x P)
    % sCov  - Sample covariance matrix of the data (size: P x P)
    % P     - Dimensionality of the samples
    % N     - Number of samples
    % w     - Order of the autoregressive process
    % K     - Box constraint parameter vector

    % Output:
    % alpha - projected conditional likelihood estimate of the Gohberg-Semencul coefficients

    
    % Preallocate the matrix to store cumulative sums
    S_tilde = zeros(w + 1, w + 1);
    
    % Compute the first row of S_tilde
    for i = 1:w + 1
        S_tilde(1, i) = sum(diag(sCov(1:(P-w), i:i+(P-w-1))));
    end
    
    % Compute the remaining upper triangular part of S_tilde using the previous row's values
    for i = 2:w + 1
        for j = i:w + 1
            S_tilde(i, j) = S_tilde(i-1, j-1) - sCov(i-1, j-1) + sCov(P-w-1 + i, P-w-1 + j);
        end
    end
    
    % Make S_tilde symmetric by adding its transpose and subtracting the diagonal elements
    S_tilde = S_tilde + S_tilde' - diag(diag(S_tilde));
    
    % Extract S_x as the flipped upper-left N_a x N_a submatrix of S_tilde
    S_x = flip(flip(S_tilde(1:w, 1:w), 1), 2);
    
    % Extract c_x as the flipped last row of C_tilde corresponding to the (N_a+1)-th element
    s_x = flip(S_tilde(w + 1, 1:w))';
    
    % Perform Cholesky decomposition and invert C_x to solve for a_hat
    R = chol(S_x);
    Inv = R\(R'\diag(ones(w,1)));
    a_hat = Inv * s_x;
    
    % Compute the residual variance sigma_sq using the computed coefficients
    sigma_sq = 1/(P - w) * (S_tilde(w + 1, w + 1) - s_x' * a_hat);
    
    % Allocate memory for alpha and compute its values based on a_hat and sigma_sq
    alpha = zeros(w + 1, 1);
    alpha(1) = 1 / sigma_sq;
    alpha(2:end) = -a_hat * alpha(1);

    % Project alpha onto the box constraints using the proj function
    alpha = proj(alpha,[0;K]);
    
end


function x = proj(y, K)
    % Inputs:
    % y - Input vector to be projected
    % K - Box constraint parameter vector
    %
    % Output:
    % x - Projected vector onto the box constraints

    % Initialize output vector x
    x = zeros(length(y), 1);
    
    % Apply lower bound constraint on the first element of y
    x(1) = max(y(1), 1e-5);
    
    % Apply box constraints on the remaining elements of y
    for i = 2:length(y)
        if y(i) < -K(i) * x(1)
            x(i) = -K(i) * x(1);
        elseif y(i) > K(i) * x(1)
            x(i) = K(i) * x(1);
        else
            x(i) = y(i);
        end
    end
end