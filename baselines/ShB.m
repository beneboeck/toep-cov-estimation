function [Gamma_est,Cov_est] = ShB(X, P, N, sCov)

    % Inputs:
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples

    % Output:
    % Cov_est   - Estimated covariance matrix
    % Gamma_est - Estimated inverse covariance matrix

    H = ones(P,1) * ones(1,P) - diag(ones(P,1));

    T = trace(sCov)/P * diag(ones(P,1)) + (sCov(:).'*H(:))/(P*(P-1)) * H;


    alpha_g = (N)/((N-1) * (N+2)) * ((N+1) * trace(sCov)^2 - 2 * sCov(:).'*sCov(:));
    beta_g = (N)/((N-1) * (N+2)) * (N * sCov(:).'*sCov(:) - trace(sCov)^2);
    c_lambda = (ones(1,P) * sCov * ones(P,1))^2 - 2 * (ones(1,P) * sCov * sCov * ones(P,1)) + sCov(:).'*sCov(:);
    lambda_g = (N)/((N-1) * (N+2)) * ((N+1) * (sCov(:).'*H(:))^2 - 2*c_lambda);


    C_diff = sCov - T;
    w_g = 1 - (beta_g - 1/P * alpha_g - lambda_g/(P*(P-1)))/(C_diff(:).'*C_diff(:));
    w_g = min(1,w_g);
    w_g = max(0,w_g);

    Cov_est = (1 - w_g) * sCov + w_g * T;
    if cond(Cov_est) > 10^5
        Gamma_est = inv(Cov_est + 0.01 * diag(ones(P,1)));
    else
        Gamma_est = inv(Cov_est);
    end
end
