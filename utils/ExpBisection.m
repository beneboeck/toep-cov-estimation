function eta = ExpBisection(P, w, lam)
    % Initialize parameters for bisection
    eta_start = 1e-5;
    eta_end = 10;
    tol = 1e-4;
    max_iter = 100;

    % Compute initial values for B_start and B_end
    B_start = computeB(ExpBounds(eta_start, P, w, lam));
    B_end = computeB(ExpBounds(eta_end, P, w, lam));

    % Bisection loop
    for n = 1:max_iter
        % Calculate the midpoint eta
        eta = (eta_end + eta_start) / 2;
        
        % Compute B_eta at the current midpoint eta
        B_eta = computeB(ExpBounds(eta, P, w, lam));

        % Check convergence condition
        if (1 - B_eta < tol) && (1 - B_eta > 0)
            return;
        end
        
        % Update the interval based on the sign of (1 - B_eta)
        if sign(1 - B_eta) == sign(1 - B_start)
            eta_start = eta;
            B_start = B_eta;  % Update B_start to B_eta
        else
            eta_end = eta;
            B_end = B_eta;  % Update B_end to B_eta
        end
    end
end