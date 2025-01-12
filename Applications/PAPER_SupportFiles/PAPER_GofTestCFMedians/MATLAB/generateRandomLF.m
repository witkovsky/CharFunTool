function randomMatrix = generateRandomLF(n, N, a)
    % Generate random numbers using the inverse transform method
    u = rand(n, N);  % Generate uniform random numbers between 0 and 1
    
    % Invert the CDF to get the random numbers
    randomMatrix = sqrt(-2 * log(1 - u) / a);
end