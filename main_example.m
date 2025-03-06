clear;
clc;

addpath('cov_generators\');
addpath('our_estimators\');
addpath('utils\');
addpath('baselines\');

%% define global variables

P = 64; % dimension of the process
N = 32; % samples considered 
w_max = P-1; % maximal autoregressive order / covariance bandwidth
N_test = 200; % number test samples
cv_factor = 4; % cross validation ratio

%% generate Ground Truth Covariance Matrix 

%C = generate_ARMA11(P,0.8,0.3,0.3);
C = generate_AR(P,0.8,[0.5,0.2,0.05]);
G = inv(C);

nMSEC = zeros(1, N_test);
nMSEG = zeros(1, N_test);

[V,D] = eig(C);

for i = 1:N_test
    X_iid = randn(N,P); % N x P
    X_transpose = V * sqrt(D) * X_iid'; % N x P
    X_data = X_transpose';
    sCov = 1/N * (X_data' * X_data);
    
    %%% COMMENT IN YOUR DESIRED ESTIMATOR %%%%

    %% TEST OUT EITHER THE PROPOSED PLS OR PGD FROM OUR PAPER:
    [G_est, C_est, memory, la] = hyparaTuningPLS(X_data, P, N, sCov, w_max);
    %[G_est, C_est, memory, la] = hyparaTuningPGD(X_data, P, N, sCov, w_max);

    nMSEC(i) = sum((C_est(:) - C(:)).^2) / sum(C(:).^2);
    nMSEG(i) = sum((G_est(:) - G(:)).^2) / sum(G(:).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TEST OUT BASELINES THAT GUARANTEE POSITIVE DEFINITENESS:
    %[G_est, C_est] = Circ(X_data, P, N, sCov);
    %[G_est, C_est] = EM(X_data, P, N, sCov);
    %[G_est, C_est] = TSL(X_data, P, N, sCov);
    %[G_est, C_est] = ShB(X_data, P, N, sCov);

    %nMSEC(i) = sum((C_est(:) - C(:)).^2) / sum(C(:).^2);
    %nMSEG(i) = sum((G_est(:) - G(:)).^2) / sum(G(:).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TEST OUT BASELINES THAT DO NOT ENFORCE POSITIVE DEFINITENESS:
    %C_est = Avg(X_data, sCov, P, N);
    %[C_est,win_length] = hyparaTuningBanding(X_data, P, N, sCov, w_max, cv_factor);
    %[C_est,win_length] = hyparaTuningTapering(X_data, P, N, sCov, w_max, cv_factor);
    %C_est = ShU(X_data, P, N, sCov);

    %nMSEC(i) = sum((C_est(:) - C(:)).^2) / sum(C(:).^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end

fprintf('Average nMSE Covariance: %.4f\n', mean(nMSEC));
fprintf('Average nMSE Inverse Covariance: %.4f\n', mean(nMSEG));
