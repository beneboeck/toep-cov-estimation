function B = computeB(K)
    % Calculate the dimension P as length(K) + 1
    P = length(K) + 1;
    
    % Initialize the generalized Fibonacci sequence fib with zeros
    fib = zeros(P, 1);
    fib(1) = 1;  % Set the first element of fib to 1 (fib starts with index 0 in the paper)
    
    % Calculate the generalized Fibonacci sequence
    for i = 1:length(K)
        fib(i+1) = sum(K(1:i) .* flip(fib(1:i)));
    end
    
    % Initialize g with zeros
    g = zeros(length(K), 1);
    
    % Calculate the g sequence
    K_rev = flip(K);  % Reverse K
    for k = 1:length(K)
        g(k) = sum(K_rev(1:k) .* flip(fib(1:k)));
    end
    
    % Calculate the final B
    B = sum((P-1:-1:1)' .* g.^2);
end