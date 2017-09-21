function [coefs,scale,shift,x,w,P] = LegendreSeries(fun,A,B,nMax,tol,x,w,P)
% LegendreSeries calculates the coefficients of the n-th order (less or
%  equal than given nMax) FOURIER-LEGENDRE SERIES expansion of the function
%  FUN over the interval [A,B]. The calculated coefficients are truncated
%  such that abs(coefs) > tol. 
%
%  The output of LegendreSeries can be further used, e.g., by the algorithm
%  LegendreSeriesFourierIntegral for computing the Fourier Integral of
%  function FUN approximated by the n-th order FOURIER-LEGENDRE SERIES
%  expansion over [A,B]. 
%
% SYNTAX:
%   [coefs,scale,shift,x,w,P] = LegendreSeries(fun,A,B,nMax,tol,x,w,P)
%
% EXAMPLE: (PDF by inverting CF of ChiSquared distribution)
%  x   = linspace(1e-12,15,501);
%  cf  = @(t)cf_ChiSquare(t,3);
%  fun = @(t) real(cf(t));
%  A   = 1e-12;
%  B   = 50;
%  [coefs,scale,shift] = LegendreSeries(fun,A,B);
%  FI  = LegendreSeriesFourierIntegral(x,coefs,scale,shift);
%  PDF = 2*real(FI)/pi;
%  plot(x,PDF,'.-');grid

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Sep-2017 09:33:32

%% ALGORITHM CALL
% [coefs,scale,shift,x,w,P] = LegendreSeries(fun,A,B,nMax,tol,x,w,P)

%% CHECK THE INPUT PARAMETERS
narginchk(1,8);
if nargin < 8, P    = []; end
if nargin < 7, w    = []; end
if nargin < 6, x    = []; end
if nargin < 5, tol  = []; end
if nargin < 4, nMax = []; end
if nargin < 3, B    = []; end
if nargin < 2, A    = []; end

if isempty(nMax), nMax = 100; end
if isempty(B), B   = 100; end
if isempty(A), A   = 0; end

if isempty(tol)
    tol = 1e-8; 
end

% Nodes and weights of the n-th order Gauss-Legendre quadrature on [-1,1]
if isempty(x) || isempty(w)
    [x,w]  = LegendrePoints(nMax+1);
end
% The Legendre polynomials of orders 0:nPts evaluated at the nodes x
if isempty(P)
    P      = zeros(nMax+1);
    P(1,:) = 1;
    p1     = 1;
    p2     = 0;
    for k  = 1:nMax
        p3 = p2;
        p2 = p1;
        p1 = ((2*k-1) .* x .* p2 - (k-1) * p3) / k;
        P(k+1,:) = p1;
    end
end

%% ALGORITHM

% Function fun (weighted, shifted and scaled) evaluated at x
scale = (B - A) / 2;
shift = (B + A) / 2;
F     = fun(scale*x + shift);
coefs = F.*w;

% Coefficients (greater than 1e-14) of the Legendre Series
coefs = sum(bsxfun(@times,P,coefs),2);
coefs = coefs(abs(coefs)>tol);

end