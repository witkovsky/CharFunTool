function random_numbers = generate_random_numbers_S3(k, n, N)
    % Generisanje slučajnih brojeva sa zadatom CDF
    u = rand(n, N);
    random_numbers=zeros(n,N);
    
    % Primena inverzne transformacije na svaki slučajni broj
    for i = 1:n
        for j=1:N
        if u(i,j) < 0.5
            random_numbers(i,j) =0.5 - ((0.5-u(i,j)) / 2^(k-1))^(1/k);
        else
            random_numbers(i,j) = 0.5 + ((u(i,j) -0.5) / 2^(k-1))^(1/k);
        end
        end 
    end
end