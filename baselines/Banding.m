function C = Banding(k,c)

    % Inputs:
    % k         - covariance bandwidth for the banding
    % c         - first row of the Toeplitz covariance that is supposed to banded
    %
    % Output:
    % C         - Estimated covariance matrix

    w = zeros(length(c),1);
    w(1:k+1) = 1;
    C = toeplitz(w .*c);
end