function [Gamma_est,Cov_est] = EM(X, P, N, sCov)

    % Inputs:
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples

    % Output:
    % Cov_est   - Estimated covariance matrix
    % Gamma_est - Estimated inverse covariance matrix

    % some constant parameters
    max_iter = 500;

    % initialize DFT matrix
    F = 1/sqrt(2 * P) * dftmtx(2 * P);
    
    % extract the first P columns
    Fp = F(:,1:P);

    % initialize the covariance estimate with a large value
    C_p_old = 10^(8) * ones(P,P);
    % initialize c vector with ones
    c = ones(2 * P,1);

    % start iterating through the EM algorithm
    for it = [1:max_iter]

        % compute circulant 2P-dimensional covariance
        C_2p = ifft(fft((diag(c))').');
        % extract P-dimensional covariance
        C_p = real(C_2p(1:P,1:P));
        % regularize if ill-conditioned
        if cond(C_p) > 10^5
            C_p = C_p + 0.01 * diag((ones(P,1)));
        end
        
        % invert the P-dimensional covariance
        % for a dimension greater 64 the concatenation of levinson & ICM generation is faster than the naive inv algorithm in general, it is also more stable
        if P > 1024
            c_p =  C_p(:,1);
            [a,e] = levinson(c_p,P-1);
            alpha = [1/e,a(2:end)/e]';
            C_p_inv = gen_Gamma_varA(alpha,length(alpha));
        else
            C_p_inv = inv(C_p);
        end

        % update c vector
        Cwhole = C_p_inv * sCov * C_p_inv';
        Cwhole_big = [Cwhole,zeros(P,P);zeros(P,2*P)];
        C_inv_big = [C_p_inv,zeros(P,P);zeros(P,2*P)];
        four_Cwhole_big = fft(ifft(Cwhole_big').');
        four_C_inv_big = fft(ifft(C_inv_big').');
        C_aux = c .* four_Cwhole_big .* c' + diag(c)  - c .* four_C_inv_big .* c';
        c = real(diag(C_aux));

        % termination criterion
        if 1/(P^2) * norm(C_p_old - C_p,'fro')^2 < 10^(-4)
            break;
        end
        C_p_old = C_p;

    end
    
    % compute final covariance and inverse covariance
    C_2p = ifft(fft((diag(c))').');
    Cov_est = real(C_2p(1:P,1:P));
    if cond(Cov_est) > 10^5
        Cov_est = Cov_est + 0.01 * diag(ones(P,1)); 
    end
    
    if P > 64
        c_p =  Cov_est(:,1);
        [a,e] = levinson(c_p,P-1);
        alpha = [1/e,a(2:end)/e]';
        Gamma_est = gen_Gamma_varA(alpha,length(alpha));
    else
        Gamma_est = inv(Cov_est);
    end
end