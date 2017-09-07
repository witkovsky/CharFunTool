function [pval,pcoef] = ChebPoly(n,x)
% ChebPoly Evaluates the Chebyshev polynomial of the first kind T_n(x).
%   See https://en.wikipedia.org/wiki/Chebyshev_polynomials.
%
% SYNTAX:
%   [pval pcoef] = ChebPoly(n,x)
%
% EXAMPLE 1
%   n    = 0:5;
%   x    = ChebPoints(32,[-1,1])';
%   pval = ChebPoly(n,x);
%   plot(x,pval,'.-')
%   title('Chebyshev polynomials T_n(x) evaluated at Chebyshev points')
%
% MOTIVATED BY:
%   ChebyshevPoly.m by David Terr, Raytheon, 5-10-04

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Sep-2017 16:20:17

%% FUNCTION
%  [pval pcoef] = ChebPoly(n,x)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, x = []; end

if isempty(x)
    x = linspace(-1,1,101);
end

szx = size(x);
x   = x(:)';
Nx  = length(x);

%% ALGORITHM

Nn    = length(n);
pval  = zeros(Nn,Nx);
pcoef = cell(1,Nn);
for k = 1:Nn
    nk = n(k);
    if nk == 0
        coef = 1;
    elseif nk == 1
        coef = [1 0]';
    else
        aux1      = zeros(nk+1,1);
        aux1(nk)   = 1;
        aux2      = zeros(nk+1,1);
        aux2(nk+1) = 1;
        for i = 2:nk
            coef = zeros(nk+1,1);
            for j = nk-i+1:2:nk
                coef(j) = 2*aux1(j+1) - aux2(j);
            end
            if mod(i,2) == 0
                coef(nk+1) = (-1)^(i/2);
            end
            if i < nk
                aux2 = aux1;
                aux1 = coef;
            end
        end
    end
    pcoef{k}  = coef;
    pval(k,:) = polyval(coef,x);
end

if Nn == 1
    pval  = reshape(pval,szx);
    pcoef = coef;
end