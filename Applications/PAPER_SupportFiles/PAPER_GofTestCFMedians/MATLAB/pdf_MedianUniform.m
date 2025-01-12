function pdf = pdf_MedianUniform(x, a, b, N)
%  pdf_MedianUniform for evaluating the PDF of medians from a random
%    sample of size N from the uniform distribution U(a,b). The
%    PDF is evaluated at specified median values a <= x <= b.
%
%  The algorithm is a realization of equation (4) and (5) in Theorem 2 of
%  the PMW manuscript. 
% 
% Syntax:
%  pdf = pdf_MedianUniform(x, a, b, N)
%
% % EXAMPLE
%  a = -3;                      % Replace with your specific values
%  b = 3;                       % Replace with your specific values
%  N = 16;                      % Replace with your specific values
%  x = linspace(a, b, 1000);    % Range of median values
%  % Calculate PDF for the entire range of x
%  pdf_values = pdf_MedianUniform(x, a, b, N);
%  % Plot the PDF
%  figure
%  plot(x, pdf_values);
%  xlabel('x');
%  ylabel('PDF_{Median}(x)');
%  title(['PDF for a = ', num2str(a), ', b = ', num2str(b), ', N = ', num2str(N)]);
%  grid on;
%
% % EXAMPLE
%  a = -4;                      % Replace with your specific values
%  b = 1;                       % Replace with your specific values
%  N = 16;                      % Replace with your specific values
%  x = linspace(a, b, 1000);    % Range of m values
%  % Calculate PDF for the entire range of x
%  pdf_values = pdf_MedianUniform(x, a, b, N);
%  % Plot the PDF
%  figure
%  plot(x, pdf_values);
%  xlabel('x');
%  ylabel('PDF_{Median}(x)');
%  title(['PDF for a = ', num2str(a), ', b = ', num2str(b), ', N = ', num2str(N)]);
%  grid on;

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver. 20-Jan-2024 14:45:47

%% Algorithm

szm = size(x);
x   = x(:);
x0  = 2/(b-a) * (x - (b+a)/2);
pdf = zeros(size(x));  % Initialize the output vector

if mod(N, 2) == 0
    % Find indices where m is within the defined range
    id1 = (x0 >= -1) & (x0 < 0);
    id2 = (x0 >= 0) & (x0 <= 1);

    pdf(id1) = hfun(x0(id1),N);
    pdf(id2) = hfun(-x0(id2),N);

else
    pdf = ((N+1)/2^(N+1)) * nchoosek(N, (N+1)/2) .* (1 - x0.^2).^((N-1)/2);
end

pdf = pdf * 2/(b-a);
pdf = reshape(pdf,szm);

end
%% function hfun
function f = hfun(m,n)
% Auxiliary function

f = 0;
for k = 0:(n/2 - 1)
    for j = 0:(n/2 - 1)
        f = f + nchoosek(n/2-1, k) * nchoosek(n/2-1, j) ...
            .* (1 - 2*m).^(n/2-j-1) ...
            .* (m.^(k+j+1) - (-1)^(k+j+1))/(k+j+1);
    end
end

f = nchoosek(n, n/2) * (n^2/2^(n+1)) * f;
end