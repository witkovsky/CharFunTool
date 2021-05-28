function [pval,x,domain] = ChebPolyValues(coeffs,x,domain)
% ChebPolyValues Evaluates the values of the Chebyshev polynomial specified
% by its coefficients at the specified grid of points x, where x = xMin +
% xx*(xMax-xMin), with domain = [xMin, xMax] and xx being the grid of
% points from the interval [-1,1].   
%
% SYNTAX:
%   [pval,x,domain] = ChebPolyValues(coeffs,x,domain)
%
% EXAMPLE1 (Chebyshev polynomial values specified by coefficients of Sine)
%   n      = 2^5+1;
%   domain = [-pi,pi];
%   nodes  = ChebPoints(n,domain);
%   f      = sin(nodes);
%   coeffs = ChebCoefficients(f);
%   x      = linspace(-pi,pi)';
%   pval   = ChebPolyValues(coeffs,x);
%   figure; plot(x,pval,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev polynomial')
%   title('Chebyshev Polynomial Specified by its Coefficients')
%
% EXAMPLE2 (Chebyshev polynomial values of the Sine and Cosine functions)
%   n      = 2^5+1;
%   domain = [-pi,pi];
%   nodes  = ChebPoints(n,domain);
%   f      = [sin(nodes) cos(nodes) sin(nodes).*cos(nodes)];
%   coeffs = ChebCoefficients(f);
%   [pval,x] = ChebPolyValues(coeffs,[],domain);
%   figure; plot(x,pval,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev polynomial')
%   title('Chebyshev Polynomial Specified by its Coefficients')
%
% EXAMPLE3 (Chebyshev polynomial values of the Normal PDF and CDF)
%   n      = 2^7+1;
%   nodes  = ChebPoints(n,[-8,8]);
%   f      = [normpdf(nodes) normcdf(nodes)];
%   coeffs = ChebCoefficients(f);
%   x      = linspace(-4,4,101);
%   domain = [-8, 8];
%   [pval,x] = ChebPolyValues(coeffs,x,domain);
%   figure; plot(x,pval,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev polynomial')
%   title('Chebyshev Polynomial Specified by its Coefficients')

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-May-2021 14:28:24

%% FUNCTION
%  [pval,x,domain] = ChebPolyValues(coeffs,x,domain)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, domain = []; end
if nargin < 2, x = []; end

szx = size(x);
x   = x(:);

% Set x and the domain such that x in [-1,1] and domain = [xMin, xMax]
if isempty(domain)
    if isempty(x)
        domain = [-1,1];
        x = linspace(-1,1,101)';
    else
        domain = [min(x),max(x)];
        x = 2*(x(:) - domain(1))/ (domain(2)-domain(1)) - 1;
    end
else
    if isempty(x)
        x = linspace(-1,1,101)';
    else
        if (min(x) < domain(1) || max(x) > domain(2))
            warning('Some values are outside the stated domain')
        end
        x = 2*(x(:) - domain(1))/ (domain(2)-domain(1)) - 1;
    end
end

n = 0:length(coeffs(:,1))-1;

%% ALGORITHM

pval  = cos(n.*acos(x))*coeffs;
x     = domain(2)*(x + 1)/2 + domain(1)*(1 - x)/2;

if all(szx) && size(pval,2) == 1
    pval = reshape(pval,szx);
end

end