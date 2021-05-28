function values = ChebValues(coeffs)
% ChebValues converts the n-dimensional Chebyshev coefficients to
%  n-dimensional function values, say f(x), evaluated at the n-dimensional
%  Chebyshev points x of the 2nd kind. 
%  
%  V = ChebValues(C) returns the n-dimensional vector of (polynomial)
%  values evaluated at the Chebyshev points x such that V(i) = f(x(i)) =
%  C(1)*T_{0}(x(i)) + C(2)*T_{1}(x(i)) + ... + C(N)*T_{N-1}(x(i)) (where T_k(x)
%  denotes the k-th 1st-kind Chebyshev polynomial, and x(i) are the
%  2nd-kind Chebyshev nodes. 
%
%  If the input C is an (n x m)-matrix then V = ChebValues(C) returns the
%  (n x m)-matrix of values V such that V(i,j) = P_j(x_i) =
%  C(1,j)*T_{0}(x(i)) + C(2,j)*T_{1}(x(i)) + ... + C(N,j)*T_{N-1}(x(i)).
% 
%  The ChebValues algorithm is based on the ChebFun VALS2COEFFS algorithm.
%  For more details see http://www.chebfun.org/.  
%  
% SYNTAX:
%   values = ChebValues(coeffs)
%
% EXAMPLE1 (Values of Sine function evaluated Chebyshev points on (-pi,pi))
%   n      = 2^5+1
%   domain = [-pi,pi];
%   x      = ChebPoints(n,domain);
%   f      = sin(x);
%   coeffs = ChebCoefficients(f);
%   V      = ChebValues(coeffs);
%   disp([x coeffs f V])
%
% EXAMPLE2 (Chebyshev values of the Sine and the Cosine on (-pi,pi))
%   n      = 2^5+1
%   domain = [-pi,pi];
%   x      = ChebPoints(n,domain);
%   f      = [sin(x) cos(x)];
%   coeffs = ChebCoefficients(f)
%   V      = ChebValues(coeffs);
%   disp([x coeffs f V])

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-May-2021 14:28:24
% Revisions: 24-Jul-2017 10:06:48
%
% Based on the algorithm VALS2COEFFS of ChebFun.
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%% FUNCTION
%  coeffs = ChebValues(values)

%% ALGORITHM

n = size(coeffs, 1);

if ( n <= 1 )
    values = coeffs; 
    return
end

% Check symmetry
isEven = max(abs(coeffs(2:2:end,:)),[],1) == 0;
isOdd  = max(abs(coeffs(1:2:end,:)),[],1) == 0;

% Scaling
coeffs(2:n-1,:) = coeffs(2:n-1,:)/2;

% DCT using an FFT
tmp = [ coeffs ; coeffs(n-1:-1:2,:) ];

if isreal(coeffs)                       % Real-valued case
    values = real(fft(tmp));
elseif isreal(1i*coeffs)                % Imaginary-valued case
    values = 1i*real(fft(imag(tmp)));
else                                    % General case
    values = fft(tmp);
end

% Flip and truncate:
values = values(n:-1:1,:);

% Symmetry
values(:,isEven) = (values(:,isEven)+ flipud(values(:,isEven)))/2;
values(:,isOdd)  = (values(:,isOdd) - flipud(values(:,isOdd)))/2;

end