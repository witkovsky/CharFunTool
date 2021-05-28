function coeffs = ChebCoefficients(values)
% ChebCoefficients converts function values, say values = f(x), evaluated
%  at Chebyshev points of the 2nd kind, say n-dimensional x, to
%  n-dimensional Chebyshev coefficients.  
%  
%  C = ChebCoefficients(D) returns the n-dimensional vector of
%  coefficients such that f(x) = C(1)*T_{0}(x) + C(2)*T_{1}(x)+ ... +
%  C(N)*T_{n-1}(x) (where T_{k}(x) denotes the k-th 1st-kind Chebyshev
%  polynomial) interpolates the data [D(1); ... ; D(n)] at Chebyshev
%  points x = [x(1); ... ; x(n)] of the 2nd kind. 
% 
%  If D is an (n x m)-matrix, then C = ChebCoefficients(D) returns the
%  (n x m)-matrix of coefficients C.
% 
%  The ChebCoefficients algorithm is based on the ChebFun vals2coeffs
%  algorithm. For more details see http://www.chebfun.org/.  
%  
% SYNTAX:
%   coeffs = ChebCoefficients(values)
%
% EXAMPLE1 (Chebyshev coefficients of the Sine function on (-pi,pi))
%   n      = 2^5+1;
%   domain = [-pi,pi];
%   x      = ChebPoints(n,domain);
%   f      = sin(x);
%   coeffs = ChebCoefficients(f);
%   disp([x coeffs])
%   x      = linspace(-pi,pi)';
%   pval   = ChebPolyValues(coeffs,x,domain);
%   figure; plot(x,pval,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev polynomial')
%   title('Chebyshev Polynomial Values Specified by its Coefficients')
%
% EXAMPLE2 (Chebyshev coefficients of the Sine and the Cosine on (-pi,pi))
%   n      = 2^5+1;
%   domain = [-pi,pi];
%   x      = ChebPoints(n,[-pi,pi]);
%   f      = [sin(x) cos(x)];
%   coeffs = ChebCoefficients(f);
%   disp([x coeffs])
%   x      = linspace(-pi,pi)';
%   pval   = ChebPolyValues(coeffs,x,domain);
%   figure; plot(x,pval,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev polynomial')
%   title('Chebyshev Polynomial Values Specified by its Coefficients')

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-May-2021 14:28:24
% Revisions: 24-Jul-2017 10:06:48
%
% Based on the algorithm VALS2COEFFS of ChebFun.
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% FUNCTION
%  coeffs = ChebCoefficients(values)

%% ALGORITHM

n = size(values, 1);

if ( n <= 1 )
    coeffs = values; 
    return
end

% Check symmetry
isEven = max(abs(values-flipud(values)),[],1) == 0;
isOdd  = max(abs(values+flipud(values)),[],1) == 0;

% DCT using an FFT
tmp = [values(n:-1:2,:) ; values(1:n-1,:)];

if isreal(values)               % Real-valued case
    coeffs = ifft(tmp);
    coeffs = real(coeffs);
elseif isreal(1i*values)        % Imaginary-valued case
    coeffs = ifft(imag(tmp));
    coeffs = 1i*real(coeffs);
else                            % General case
    coeffs = ifft(tmp);
end

% Truncate, scale the interior coefficients, and adjust for symmetry
coeffs                 = coeffs(1:n,:);
coeffs(2:n-1,:)        = 2*coeffs(2:n-1,:);
coeffs(2:2:end,isEven) = 0;
coeffs(1:2:end,isOdd)  = 0;

end