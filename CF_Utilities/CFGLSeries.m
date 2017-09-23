function [nCoefs,pdfCoefs,cdfCoefs,scale,shift,x,w,P] = ...
         CFGLSeries(cf,A,B,nMax,tol,x,w,P)
% CFGLSeries (Gauss-Legendre series expansion)calculates the PDF and CDF
%  coefficients of the n-th order (less or equal than given nMax)
%  GAUSS-LEGENDRE SERIES EXPANSION of the characteristic function CF over
%  the interval [A,B]. The calculated coefficients are truncated such that
%  abs(coefs) > tol.  
%
%  The output of CFGLSeries can be further used, e.g., by the algorithm
%  CFGLSeriesFourierIntegral for computing the Fourier Integral of
%  function FUN approximated by the n-th order GAUSS-LEGENDRE SERIES
%  EXPANSION over [A,B]. 
%
% SYNTAX:
%   [nCoefs,pdfCoefs,cdfCoefs,scale,shift,x,w,P] = ...
%                                      CFGLSeries(cf,A,B,nMax,tol,x,w,P)
%
% EXAMPLE: (PDF by inverting CF of ChiSquared distribution)
%  x   = linspace(1e-12,15,501);
%  cf  = @(t)cf_ChiSquare(t,3);
%  A   = 1e-12;
%  B   = 50;
%  [nCoefs,pdfCoefs,cdfCoefs,scale,shift] = CFGLSeries(cf,A,B)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 22-Sep-2017 11:11:11

%% ALGORITHM CALL
% [nCoefs,pdfCoefs,cdfCoefs,scale,shift,x,w,P] = ...
%                                     CFGLSeries(cf,A,B,nMax,tol,x,w,P)

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
    tol = 1e-15; 
end

% Nodes and weights of the n-th order Gauss-Legendre quadrature on [-1,1]
if isempty(x) || isempty(w)
    [x,w]  = LegendrePoints(nMax+1);
end
% The Legendre polynomials of orders 0:nCoefs evaluated at the nodes x
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
scale    = (B - A) / 2;
shift    = (B + A) / 2;
t        = scale * x + shift;
F        = cf(t);
pdfCoefs = real(F) .* w;
cdfCoefs = imag(F) .* w ./ t;
pdfCoefs = sum(bsxfun(@times,P,pdfCoefs),2);
cdfCoefs = sum(bsxfun(@times,P,cdfCoefs),2);

pdfCoefs = pdfCoefs(abs(pdfCoefs)>tol);
cdfCoefs = cdfCoefs(abs(cdfCoefs)>tol);

% Coefficients (greater than tol = 1e-15) of the Legendre Series
if ~isempty(pdfCoefs) && ~isempty(cdfCoefs)
    nCoefs = max(length(pdfCoefs),length(cdfCoefs));
elseif ~isempty(pdfCoefs) && isempty(cdfCoefs)
    nCoefs = length(pdfCoefs);
elseif isempty(pdfCoefs) && ~isempty(cdfCoefs)
    nCoefs = length(cdfCoefs);
else
    nCoefs = 0;
end

if ~isempty(pdfCoefs) && length(pdfCoefs) < nCoefs
    coef0 = pdfCoefs;
    pdfCoefs = zeros(nCoefs,1);
    pdfCoefs(1:length(coef0)) = coef0;
end
if ~isempty(cdfCoefs) && length(cdfCoefs) < nCoefs
    coef0 = cdfCoefs;
    cdfCoefs = zeros(nCoefs,1);
    cdfCoefs(1:length(coef0)) = coef0;
end

end