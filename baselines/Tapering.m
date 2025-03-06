function C = Tapering(k,c,P)
    
    % Inputs:
    % k         - covariance bandwidth for the banding
    % c         - first row of the Toeplitz covariance that is supposed to banded
    % P         - sample dimension
    %
    % Output:
    % C         - Estimated covariance matrix

    s = length(c);

    % compute the window (1, linearly decreasing, 0)
    w = zeros(s(1),1);
    w(1:k+1) = 1;
    for i = [k+2:min(2*k+1,P)]
        w(i) = (2 - (i-1)/(k));
    end
    % Computing the tapered Covariance Matrix
    C = toeplitz(c .* w);
end