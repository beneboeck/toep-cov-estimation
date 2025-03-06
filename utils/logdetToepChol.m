function logdetT = logdetToepChol(c)
    % from Paper "On the stability of the Bareiss and realted Toeplitz factorization algorithms"
    P = length(c);
    % generators of the toeplitz matrix (utilizes the displacement rank 2 property of toeplitz matrices)

    break_c = false;
    uk = c/sqrt(c(1));
    vk = uk;
    vk(1) = 0;

    diagL = zeros(P,1);
    for k = [1:P-1]
        ek = zeros(P,1);
        ek(k) = 1;
        ekp1 = zeros(P,1);
        ekp1(k+1) = 1;

        sinth = (ekp1'*vk)/(ek'*uk);
        costh = sqrt(1-sinth^2);
        
        Zuk = [0;uk(1:end-1)];

        ukp1 = (Zuk - sinth*vk)/costh;
        vkp1 = -sinth*ukp1 + costh * vk;

        diagL(k) = uk(k);
        if (k > 1) & (diagL(k-1) == diagL(k))
            break_c = true;
            diagL(k+1:end) = diagL(k);
            break;
        end
        uk = ukp1;
        vk = vkp1;

    end
    if break_c == false
        diagL(P) = uk(P);
    end
    logdetT = 2 * sum(log(diagL));
end