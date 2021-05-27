function [fun,x,domain] = InterpChebValues(values,x,domain)
% InterpChebValues Evaluates the function fun by using Chebyshev
%   polynomial approximation specified by the function values evaluated at
%   Chebyshev points with specified domain.   
%
% SYNTAX:
%   [fun,x,domain] = InterpChebValues(values,x,domain)
%
% EXAMPLE1 (Interpolation by Chebyshev polynomial approximation of Sine)
%   n      = 2^5+1;
%   nodes  = ChebPoints(n,[-pi,pi]);
%   values = sin(nodes);
%   x      = linspace(-pi,pi)';
%   domain = [-pi,pi];
%   fun    = InterpChebValues(values,x,domain)
%   figure; plot(x,fun,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev Interpolant')
%   title('Chebyshev Polynomial Approximation')
%
% EXAMPLE2 (Interpolation by Chebyshev polynomial approximation)
%   n      = 2^5+1;
%   nodes  = ChebPoints(n,[-pi,pi]);
%   values = [sin(nodes) cos(nodes) sin(nodes).*cos(nodes)];
%   x      = linspace(-pi,pi)';
%   domain = [-pi,pi];
%   fun    = InterpChebValues(values,x,domain)
%   figure; plot(x,fun,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev Interpolant')
%   title('Chebyshev Polynomial Approximation')
%
% EXAMPLE3 (CDF/PDF interpolation by Chebyshev polynomial approximation)
%   n      = 2^7+1;
%   nodes  = ChebPoints(n,[-8,8]);
%   values = [normpdf(nodes) normcdf(nodes)];
%   x      = linspace(-4,4,101);
%   domain = [-8, 8];
%   fun    = InterpChebValues(values,x,domain);
%   figure; plot(x,fun,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev Interpolant')
%   title('Chebyshev Polynomial Approximation')

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 27-May-2021 17:56:16

%% FUNCTION
%  [fun,x,domain] = InterpChebValues(values,x,domain)

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
            error('Dimension mismatch')
        else
            x = 2*(x(:) - domain(1))/ (domain(2)-domain(1)) - 1;
        end
    end
end

coeffs = ChebCoefficients(values);
n = 0:length(coeffs(:,1))-1;

%% ALGORITHM

fun  = cos(n.*acos(x))*coeffs;
x     = domain(2)*(x + 1)/2 + domain(1)*(1 - x)/2;

if size(fun,2) == 1
    fun = reshape(fun,szx);
end

end