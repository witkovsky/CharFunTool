function pval = ChebPoly(n,x)
% ChebPoly Evaluates the Chebyshev nth order polynomial of the first kind
%   T_n(x), where T_n(x) = cos(n*acos(x)). For more details see also
%   https://en.wikipedia.org/wiki/Chebyshev_polynomials.  
%
% SYNTAX:
%   pval = ChebPoly(n,x)
%
% EXAMPLE 1
%   n    = 0:5;
%   x    = ChebPoints(32,[-1,1])';
%   pval = ChebPoly(n,x);
%   plot(x,pval,'.-')
%   title('Chebyshev polynomials T_n(x) evaluated at Chebyshev points')

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Sep-2017 16:20:17

%% FUNCTION
%  pval = ChebPoly(n,x)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, x = []; end

if isempty(x)
    x = linspace(-1,1,101);
end

szx = size(x);
x   = x(:)';

%% ALGORITHM
pval  = cos(n(:)*acos(x));

if length(n) == 1
    pval  = reshape(pval,szx);
end
end