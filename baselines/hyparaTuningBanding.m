function [Cov_est,win_length] = hyparaTuningBanding(X,P,N,sCov,N_a,cv_factor)
    
    % Banding Estimator for more than one training sample
    % Inputs:
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples
    % N_a       - Maximum considered covariance bandwidth
    % cv_factor - Cross-validation factor (not used in the function)
    %
    % Outputs:
    % Cov_est       - Estimated covariance matrix
    % win_length    - chosen covariance bandwidth


    %% banding with more than one given sample and hyperparameter tuning via
    % cross validation

    % initialize risks
    risks = zeros(N_a,1); % risk for each window length
    batch_size = N/cv_factor; % how many samples are in each cross validation batch

    % loop through different covariance bandwidths
    for wl = [1:N_a]
        r = 0;
        % loop through different cross validation batches
        for cv = [1:cv_factor]
            % compute the corresponding evaluation and training samples
            X_eval = X((cv-1)*batch_size + 1:cv * batch_size,:); % evaluation samples used for the toeplitz sample covariance
            X_train = X(setdiff(1:end,(cv-1)*batch_size + 1:cv * batch_size),:); % training samples used for the banding covariance

            sCov_eval = 1/batch_size * (X_eval' * X_eval); % evaluation sample cov
            sCov_train = 1/(N - batch_size) * (X_train' * X_train); % training sample cov
            
            % compute the first row of the corresponding Avg Estimators
            c_eval = zeros(P,1);
            c_train = zeros(P,1);
            for i = [1:P]
                c_eval(i) = mean(diag(sCov_eval,i-1)); 
                c_train(i) = mean(diag(sCov_train,i-1));
            end
            sC_eval = toeplitz(c_eval); % evaluation toeplitz sample cov

            % estimated MSE of banding estimator and evaluation toeplitz sample cov
            r = r + 1/cv_factor * norm(Banding(wl-1,c_train) - sC_eval,'fro')^2;
            end
            risks(wl) = r;

            % termination criterion checking (if r increased over the last 5 tries, stop)
            if wl > 5
                risks_temp = risks(wl-5:wl);
                idx_temp = find(risks_temp == min(risks_temp));
            end
        
            if (wl > 5) && (idx_temp == 1)
                risks = risks(1:wl);
                break
            end

        end
        idx_m = find(risks == min(risks));
               
        c = zeros(P,1);
            for i = [1:P]
                c(i) = mean(diag(sCov,i-1));
            end

        Cov_est = Banding(idx_m-1,c);
        win_length = idx_m;
end