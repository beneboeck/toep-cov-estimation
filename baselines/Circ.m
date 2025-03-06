function [Gamma_est,Cov_est] = Circ(X,P,N,sCov)

    % Inputs:
    % X         - Data matrix of size (N x P)
    % sCov      - Sample covariance matrix of the data (size: P x P)
    % P         - Dimensionality of the samples
    % N         - Number of samples

    % Output:
    % Cov_est   - Estimated covariance matrix
    % Gamma_est - Estimated inverse covariance matrix

    Cov_est = ifft(fft(diag(diag(fft(ifft(sCov').')))').');

    if P > 1024
        c_p =  Cov_est(:,1);
        [a,e] = levinson(c_p,P-1);
        alpha = [1/e,a(2:end)/e]';
        Gamma_est = gen_Gamma_varA(alpha,length(alpha));     
    else
        Gamma_est = inv(Cov_est);
    end
end