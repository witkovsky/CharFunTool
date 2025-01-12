function cf = cf_MedianExponential(t, lambda, N)
%  cf_MedianExponential calculates the CF of the sample median from a
%  random sample of size N from the exponential distribution EXP(lambda),
%  with lambda >  0.  
%
%  The algorithm is a realization of equation (16) in Theorem 6 of the
%  PMW manuscript.
%
% Syntax:
%   cf = cf_MedianExponential(t, lambda, N)
%
% % EXAMPLE
%  lambda = 0.1;                % Replace with your specific values
%  N = 4;                       % Replace with your specific values
%  t = linspace(-2, 2, 101);    % Range of t values
%  % Calculate CF for the entire range of t
%  cf = cf_MedianExponential(t, lambda, N);
%  % Plot the CF
%  plot(t, real(cf), t, imag(cf));
%  xlabel('t');
%  ylabel('CF \phi_{Median}(t)');
%  title('CF of the sample median from Exponential distribution');
%  legend('Real(CF)','Imag(CF)','Location','ne')
%  grid on;
%
% % EXAMPLE
%  lambda = 0.1;                % Replace with your specific values
%  N    = 4;                    % Replace with your specific values
%  cf   = @(t) cf_MedianExponential(t, lambda, N);
%  x    = linspace(0,30,101)';
%  prob = [0.9, 0.95, 0.975, 0.99];
%  clear options
%  options.xMin = 0;
%  result = cf2DistGP(cf,x,prob,options);
%  pdf = pdf_MedianExponential(x, lambda, N); 
%  figure
%  plot(x,pdf,'.-',x,result.pdf,'--')
%  xlabel('x');
%  ylabel('PDF_{Median}(x)');
%  title(['Exact PDF vs. PDF from inverted CF for \lambda = ', num2str(lambda), ', N = ', num2str(N)]);
%  legend('Exact','Inverted from CF','Location','ne')
%  grid on;

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver. 20-Jan-2024 14:45:47

%% Algorithm

cf = zeros(size(t));  % Initialize the output vector

if mod(N, 2) == 0
    for k = 0:(N/2 - 1)
        terms = nchoosek(N/2 - 1, k) .* ...
            ((-1)^(k) ./ ( (lambda*(N/2+1+k) - 1i*t) .* (lambda*N - 1i*t)) );
        cf = cf + terms;
    end
    % Multiply by precomputed constants
    cf = nchoosek(N, N/2) * (N^2/2) * lambda^2 * cf;
else
    for k = 0:(N-1)/2
        terms = nchoosek((N-1)/2, k) .* ...
            ((-1)^(k) ./ (lambda*((N+1)/2+k) - 1i*t) );
        cf = cf + terms;
    end
    % Multiply by precomputed constants
    cf = nchoosek(N, (N+1)/2) * (N+1)/2 * lambda * cf;
end

end
