function [Gamma_est, Cov_est, memory, la] = hyparaTuningPLS(X, P, N, sCov, w_max)
    % Inputs:
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples
    % w_max     - Maximum memory (order) of the autoregressive process
    % cv_factor - Cross-validation factor (not used in the function)
    %
    % Outputs:
    % Gamma_est - Estimated precision (inverse covariance) matrix
    % Cov_est   - Estimated covariance matrix
    % memory    - Optimal memory selected by BIC
    % la        - Selected exponential parameter

    % Define exponential factors and corresponding constants for K dictionary
    k_exp = [0.6, 1, 1.4, 1.8, 2.2];
    eta_exp = zeros(length(k_exp), w_max);
%     for w = 1:w_max
%         for i = 1:length(k_exp)
%             eta_exp(i,w) = ExpBisection(P,w,k_exp(i));
%         end
%     end

    %%%%%% for our paper, we applied the approximation
    for i = 1:length(k_exp)
        eta_exp(i,:) = ExpBisection(P,15,k_exp(i));
    end
    %%%%%%

    % Preallocate variables
    N_la = length(k_exp);
    BIC = inf(w_max, 1);  % Initialize BIC with large values
    G_list = zeros(w_max, N_la, P, P);
    C_list = zeros(w_max, N_la, P, P);
    likelihood_list = zeros(w_max, N_la);
    idx_l_min = zeros(w_max, 1);

    % Loop over potential autoregressive memory values
    for w = 1:w_max
        % Loop over exponential factors
        for n_exp = 1:N_la
            K_exp = ExpBounds(eta_exp(n_exp,w),P,w,k_exp(n_exp));
            
            alpha = PLS(X, sCov, P, N, w, K_exp);
            Gamma_est = gen_Gamma_varA(alpha, P);
            a_full = [alpha; zeros(P - (w + 1), 1)];  % Extend alpha with zeros
            Cov_est = toeplitz(rlevinson(a_full / a_full(1), 1 / a_full(1)));

            % Compute negative log-likelihood
            likelihood = -logdet(Gamma_est) + Gamma_est(:)' * sCov(:);

            % Save all estimators and likelihoods
            G_list(w, n_exp, :, :) = Gamma_est;
            C_list(w, n_exp, :, :) = Cov_est;
            likelihood_list(w, n_exp) = likelihood;

            % Check for early stopping condition based on likelihood
            if n_exp > 3
                likeli_temp = likelihood_list(w,n_exp-3:n_exp);
                lidx_temp = find(round(likeli_temp,4) == min(round(likeli_temp,4)));
                if length(lidx_temp) ~= 1
                    lidx_temp = lidx_temp(1);
                end
            end
        
            if (n_exp > 3) && (lidx_temp == 1)
                idx_l_min(w) = n_exp - 3;
                break
            end
        end

        % If no early stopping, select the best exponential parameter
        if idx_l_min(w) == 0
            idx_l = find(likelihood_list(w,:) == min(likelihood_list(w,:)));

            % if two estimators led to the same performance
            if length(idx_l) ~= 1
                idx_l_min(w) = idx_l(1);
            else
                idx_l_min(w) = idx_l;
            end
        end

        % Determine the best Gamma estimator based on BIC
        G_est = squeeze(G_list(w, idx_l_min(w), :, :));
        BIC(w) = -logdet(G_est) + G_est(:)' * sCov(:) + (log(N) / N) * (w+1);

        % Early stopping for BIC
        if w > 5
            BIC_temp = BIC(w-5:w);
            idx_temp = find(BIC_temp == min(BIC_temp));
        end
        
        if (w > 5) && (idx_temp == 1)
            BIC = BIC(1:w);
            break
        end
    end

    % Select the best model based on minimum BIC
    idx_min = find(BIC == min(BIC));
    Gamma_est = squeeze(G_list(idx_min, idx_l_min(idx_min), :, :));
    Cov_est = squeeze(C_list(idx_min, idx_l_min(idx_min), :, :));
    memory = idx_min;
    la = k_exp(idx_l_min(idx_min));
end