function random_numbers_S1 = generate_random_numbers_S1(k, n,N)
    % Funkcija za generisanje slučajnih brojeva sa zadatom funkcijom gustine
    % f(x) = k(1-x)^(k-1), 0 < x < 1
    
   % Generisanje slučajnih brojeva sa uniformnom raspodelom
    u = rand(n, N);
    random_numbers_S1=zeros(n,N);
    % Inverzna transformacija
   
     for i=1:n
         for j=1:N
    
    random_numbers_S1(i,j) = 1-(1-u(i,j)).^(1/k);
         end
     end
end
   