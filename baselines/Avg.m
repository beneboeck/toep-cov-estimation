function Cov_est = Avg(X, sCov, P, N)

    % Inputs:
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples

    % Output:
    % Cov_est   - Estimated covariance matrix

    c = zeros(P,1);
    % computing the means along all diagonals
    for i = [1:P]
        c(i) = mean([diag(sCov,i-1)]);
    end
    Cov_est = toeplitz(c);
end