function alpha = PGD(objfun, X, sCov, P, N, w, K)

    % This function solves the optimization problem
    % minimize -logdet(G) + tr(sCov * G)
    % s.t.     x0 > 0 & -K(i)x0 <= x(i) <= K(i)x0
    % using projected gradient descent

    % Inputs:
    % objfun    - object handle of the objective function
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples
    % w         - Order of the autoregressive process
    % K         - Box constraint parameter vector

    % Output:
    % alpha - projected likelihood estimate of the Gohberg-Semencul coefficients
       
    % some constant parameters
    MAXITER   = 2000; % iterations PGD
    tol       = 1e-6; % termination toleranz for 1) norm of gradient, 2) distance to previous iteration
    s_tol     = 1e-5; % termination toleranz for Armjio stepsize
    c         = 0.5; % backtracking factor for stepsize

    % add a 0 to K to account for alpha_0 not needing any projections
    K = [0;K];

    % initialize starting position within the box constraints
    a_initial = zeros(w+1,1);
    a_initial(1) = 1 + 10 * rand(1);
        for i = [2:w+1]
            a_initial(i) = K(i) * a_initial(1) * 1.9 * (rand(1)-0.5);
        end  
    xk = a_initial;

    % start iterations of projected gradient descent
    for k = 1:MAXITER
        % compute function and gradient
        [f_k,g_k] = objfun(xk);

        % compute search direction (projected steepest descent)
        d_k = pg_grad(-g_k,xk,K);

        % first termination criterion w.r.t. search direction
        if 1/P * norm(d_k,2)^2 < tol
            break;
        end

        % normalize search direction
        d_k = d_k/norm(d_k);

        % initial step size
        s_k = 1;

        %linesearch by armijo backtracking
        yk = xk + s_k*d_k;
        pro_yk = pg(yk,K);
        Gamma = gen_Gamma_varA(pro_yk,P);

        % O(p^2) to compute the logdet term
        f_k1 = - logdet(Gamma) + Gamma(:).'*sCov(:);

        gr_sd = d_k'*g_k;

        % backtracking while loop
        while (f_k1 > f_k + c*s_k*gr_sd) & (s_k > s_tol)
            % backstep s_k
            s_k = 0.5 * s_k;

            % evaluate function at new iterate
            yk = xk + s_k * d_k;
            pro_yk = pg(yk,K);
            Gamma = gen_Gamma_varA(pro_yk,P);
            % cheapest way of computing the objective function
            if P > 512
                % see "On the covariance determinants of moving-average and autoregressive models" Finch (1960)
                % for P -> infty (fastens up the process by a lot)
                b_full = [1;pro_yk(2:end)/pro_yk(1)];
                b_length = length(b_full);
                r = zeros(P,1);
                r(1) = sum(b_full.^2) * 1/pro_yk(1);
                for m = [1:b_length-1]
                    b1 = b_full(1:end - m);
                    b2 = b_full(m+1:end);
                    r(m + 1) = 1/pro_yk(1) * sum(b1 .* b2);
                end
                f_k1 = logdetToepChol(r) + Gamma(:).'*sCov(:);
            else
                f_k1 = - logdet(Gamma) + Gamma(:).'*sCov(:);    
            end
        end

        % update iterate
        xp = xk; % record previous iterate to test convergence
        xk = pg(yk,K);
    
    
        % termination criteria
        if s_k < s_tol
            break;
        end
        if 1/P * norm(xk-xp,2)^2 < tol
            break;
        end
    end

    % final a vector
    alpha = xk;
end


function x = pg(y,K)
% projection on box constraints
    x = zeros(length(y),1);
    if y(1) <= 0.01
        x(1) = 0.01;
    else
        x(1) = y(1);
    end

    for i = [2:length(y)]
        if y(i) < -K(i) * x(1)
            x(i) = -K(i) * x(1);
        elseif y(i) > K(i) * x(1)
            x(i) = K(i) * x(1);
        else
            x(i) = y(i);
        end
    end
end

function g = pg_grad(d,x,K)
% projection of gradient
    N = length(x);
    g = zeros(N,1);
    if x(1) <= 0.01 & d(1) < 0
        g(1) = 0;
    else
        g(1) = d(1);
    end
    for i = [2:N]
        if x(i) <= -K(i) * x(1) + 0.0001 & d(i) < 0
            g(i) = 0;
        elseif x(i) >= K(i) * x(1) - 0.0001 & d(i) > 0
            g(i) = 0;
        else
            g(i) = d(i);
        end
    end
end