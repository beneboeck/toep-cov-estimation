function Cov_est = ShU(X, P, N, sCov)

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
    T = toeplitz(c);


    alpha_g = (N)/((N-1) * (N+2)) * ((N+1) * trace(sCov)^2 - 2 * sCov(:).'*sCov(:));
    beta_g = (N)/((N-1) * (N+2)) * (N * sCov(:).'*sCov(:) - trace(sCov)^2);

    mu_g_scaled = zeros(P-1,1);

    for m_p = [1:P-1]
        subS = sCov(1:P-m_p,m_p+1:end);
        tau_m = subS(:).'*subS(:) + trace(sCov * diag(ones(P-m_p,1),m_p) * sCov * diag(ones(P-m_p,1),-m_p));
        mu_g_scaled(m_p) = 2/(P-m_p) * ( (N)/((N-1) * (N+2)) * ((N+1) * trace(sCov * diag(ones(P-m_p,1),m_p))^2 - tau_m) );
    end

    w_g = 1 - (beta_g - 1/P * alpha_g - sum(mu_g_scaled))/(trace((sCov - T) * (sCov - T)));
    w_g = min(1,w_g);
    w_g = max(0,w_g);

    Cov_est = (1 - w_g) * sCov + w_g * T;

end
