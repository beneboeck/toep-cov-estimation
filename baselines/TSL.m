function [Gamma_est,Cov_est] = TSL(X,P,N,sCov)

    % Inputs:
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples

    % Output:
    % Cov_est   - Estimated covariance matrix
    % Gamma_est - Estimated inverse covariance matrix
    
    Cov_est = zeros(P,P);
    SFnormP1_R = zeros(N,1);
    for i = [1:N]
        if N == 1
            Ri = X * X';
        else
            Ri = X(i,:)' * X(i,:);
        end
        c = zeros(P,1);
        % computing the means along all diagonals
        for j = [1:P]
            c(j) = 1/P * sum([diag(Ri,j-1)]);
        end
        Ri = toeplitz(c);
        SFnormP1_R(i) = sqrt(1 + norm(Ri,'fro')^2);
        Cov_est = Cov_est + Ri./SFnormP1_R(i);
    end
    
    Cov_est = Cov_est / (sum(1./SFnormP1_R));
    Gamma_est = inv(Cov_est);
end