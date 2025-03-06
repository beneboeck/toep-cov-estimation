function [ob,grad] = gen_Gohberg_obj(alpha,sCov)
    % Inputs:
    % alpha     - alpha vector that paramterizes the inverse covariance
    % sCov      - sample covariance matrix of the data (size: P x P)
    %
    % Outputs:
    % ob        - value of the Gaussian log-likelihood at alpha
    % grad      - gradient of the Gaussian log-likelihood at alpha

    % compute the dimension P
    s_C = size(sCov);
    P = s_C(1);

    % compute the inverse covariance
    Gamma = gen_Gamma_varA(alpha,P);

    % compute the objective
    ob = real(- logdet(Gamma) + Gamma(:).'*sCov(:));

    % compute the covariance
    if P > 4
        alpha_full = [alpha;zeros(P-length(alpha),1)];
        a_AR = alpha_full/(alpha_full(1));
        a_AR(1) = 1;
        G_inv = toeplitz(rlevinson(a_AR,1/alpha_full(1)));
    else
        G_inv = inv(Gamma);
    end

    % generate B and Z of the Gohberg-Semencul parameterization
    n_a = length(alpha);
    alpha2 = [alpha;zeros(P - n_a,1)];
    B = tril(toeplitz(alpha2));
    alpha_prime = [[0];flip(alpha2(2:end))];
    C = tril(toeplitz(alpha_prime));

    % initialize gradient
    grad = zeros(n_a,1);
    
    % compute S - G_inv
    Sdiff = sCov - G_inv;

    % compute gradient (see paper, computational complexity)
    for i = [2:n_a]
        B_diff = zeros(P,P);
        Z_diff = zeros(P,P);
        B_diff(:,i:end) = B(:,1:end-i+1);
        Z_diff(:,end-i+2:end) = C(:,1:i-1);
        grad(i) = 2/alpha2(1) * (B_diff(:).'*Sdiff(:) - Z_diff(:).'*Sdiff(:));
    end

    M_g1 = B + B.' - Gamma;
    grad(1) = 1/alpha2(1) * (M_g1(:).'*Sdiff(:));
end