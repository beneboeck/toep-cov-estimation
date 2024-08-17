function C = generate_AR(N,sigma,a)
    % This script generates the Covariance Matrix of an N-dimensional sample
    % drawn from an autoregressive process. The vector a contains the
    % parameters of the autoregressive process, i.e., X_t = a(1)X_(t-1) +
    % a(2)X_(t-2) + ... + a(end)X_(t-end) + eps, var(eps) = sigma

    %% Stationary Check
    % The matrix B is there to check if the underlying process is
    % stationary (eigenvalue check)
    B = eye(length(a));
    B = B(1:end-1,:);
    B = [a;B];

    if all(abs(eig(B)) < 1)
        disp('stable AR process')
    end

    %% Creating Precision Matrix
    % We utilize the Connection between the Gohberg-Semencul
    % Parameters and the Autoregressive Parameters to create the Precision
    % Matrix

    a0 = 1/(sigma^2);
    ar = -a * a0;
    alpha = [a0;ar'];

    G = gen_Gamma_varA(alpha,N);

    %% Creating Covariance Matrix
    C = inv(G);
end