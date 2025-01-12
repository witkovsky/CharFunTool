function pdf = pdf_MedianExponential(x, lambda, N)
%  pdf_MedianExponential for evaluating the PDF of medians from a random
%    sample of size N from the exponential distribution EXP(lambda). The
%    PDF is evaluated at specified median values x >= 0.
%
%  The algorithm is a realization of equation (11) and (12) in Theorem 5 of
%  the PMW manuscript. 
% 
% Syntax:
%  pdf = pdf_MedianExponential(m, lambda, N)
%
% % EXAMPLE
%  lambda = 0.1;            % Replace with your specific values
%  N = 4;                   % Replace with your specific values
%  x= linspace(0, 25, 100); % Range of x values
%  % Calculate PDF for the entire range of x
%  pdf = pdf_MedianExponential(x, lambda, N);
%  % Plot the PDF
%  figure
%  plot(x, pdf,'LineWidth',2,'LineStyle','-')
%  xlabel('x');
%  ylabel('PDF f_M(x)');
%  title(['PDF for \lambda = ', num2str(lambda), ', N = ', num2str(N)]);
%  grid on;
%
% % EXAMPLE (Check the integral over PDF)
%  lambda = 0.1;            % Replace with your specific values
%  N = 4;                   % Replace with your specific values
%  upperLimit = 100;
%  I = integral(@(x)pdf_MedianExponential(x, lambda, N),0,upperLimit)
%
% % EXAMPLE (Exact vs Simulated PDF)
%  lambda = 0.1;            % Replace with your specific values
%  N = 4;                   % Replace with your specific values
%  M = 10000;
%  R = exprnd(1/lambda,M, N); % Generated random values from Exponential 
%  Medians = median(R,2);
%  x = linspace(0, 25, 100);  % Range of x values
%  pdf = pdf_MedianExponential(x, lambda, N);
%  figure
%  plot(x, pdf,'LineWidth',2,'LineStyle','-.')
%  hold on
%  histogram(Medians,'Normalization','pdf');
%  xlabel('x');
%  ylabel('PDF_{Median}(x)');
%  title(['Exact vs. Empirical PDF for \lambda = ', num2str(lambda), ', N = ', num2str(N)]);
%  grid on;

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver. 20-Jan-2024 14:45:47

%% Algorithm
pdf = zeros(size(x));  % Initialize the output vector

% Find indices where m is within the defined range
id = (x >= 0);

% Calculate the PDF for valid indices
x1 = x(id);

% Calculate the PDF for each element in x

if mod(N, 2) == 0
    for k = 0:(N/2 - 2)
        terms = (-1)^(k) * nchoosek(N/2 - 1, k)  ...
            .* (exp(lambda * x1 * (N/2 -1 -k)) - 1) ./ (lambda *(N/2 -1 -k));
        pdf(id) = pdf(id) + terms ;
    end

    %pdf(id) = nchoosek(n, n/2) * (n^2/2) * lambda .* exp(-2*lambda*x1);
    term1 = (N^2/2) * nchoosek(N, N/2) * lambda^2 .* exp(-lambda * x1 *N);

    % Complete the calculation for valid indices
    pdf(id) = term1.*( pdf(id) + (-1)^(N/2 -1).* x1);
else
    pdf(id) = ((N+1)/2) * nchoosek(N, (N+1)/2) * lambda * ...
        (1 - exp(-lambda * x1)).^((N-1)/2) .* exp(-lambda * x1 *(N+1)/2);
end
    
end