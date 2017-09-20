function [coefs,nPts] = LegendreSeries(fun,n)
% LegendreSeries calculates the coefficients of the Fourier-Legendre Series
%  expansion of the function FUN over interval [0,Inf].

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Sep-2017 22:46:25

%% ALGORITHM CALL
%[coefs,nPts] = LegendreSeries(fun,nPts)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, n = []; end

if isempty(n)
    n = 100;
end

L = 10.^(-10:10);
L = [1e-12 L 1e+308];
nL = length(L);

%% ALGORITHM

% Nodes and weights of the n-th order Gauss-Legendre quadrature on [-1,1]
[x,w]  = LegendrePts(n+1);
P      = zeros(n+1);
P(1,:) = 1;

% The Legendre polynomials of orders 0:nPts evaluated at the nodes x
p1 = 1;
p2 = 0;
for k = 1:n
    p3 = p2;
    p2 = p1;
    p1 = ((2*k-1) .* x .* p2 - (k-1) * p3) / k;
    P(k+1,:) = p1;
end

% Function fun (weighted, shifted and scaled) evaluated at x
coefs = 0;
for i = 2:nL
    scale = (L(i) - L(i-1)) / 2;
    shift = (L(i) + L(i-1)) / 2;
    F     = fun(scale*x + shift);
    coefs = coefs + F.*w;
end

% Coefficients (greater than 1e-14) of the Legendre Series
coefs = sum(bsxfun(@times,P,coefs),2);
coefs = coefs(abs(coefs)>1e-14);
nPts  = length(coefs);

end