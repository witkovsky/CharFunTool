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
%   domain = [-pi,pi];
%   nodes  = ChebPoints(n,domain);
%   values = sin(nodes);
%   x      = linspace(domain(1),domain(2))';
%   fun    = InterpChebValues(values,x,domain);
%   figure; plot(x,fun,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev Interpolant')
%   title('Chebyshev Polynomial Approximation')
%
% EXAMPLE2 (Interpolation by Chebyshev polynomial approximation)
%   n      = 2^5+1;
%   domain = [-pi,pi];
%   nodes  = ChebPoints(n,domain);
%   values = [sin(nodes) cos(nodes) sin(nodes).*cos(nodes)];
%   x      = linspace(domain(1),domain(2))';
%   fun    = InterpChebValues(values,x,domain);
%   figure; plot(x,fun,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev Interpolant')
%   title('Chebyshev Polynomial Approximation')
%
% EXAMPLE3 (CDF/PDF interpolation by Chebyshev polynomial approximation)
%   n      = 2^7+1;
%   domain = [-8,8];
%   nodes  = ChebPoints(n,domain);
%   values = [normpdf(nodes) normcdf(nodes)];
%   x      = linspace(-4,4,101);
%   fun    = InterpChebValues(values,x,domain);
%   figure; plot(x,fun,'.-'); grid
%   xlabel('x')
%   ylabel('Chebyshev Interpolant')
%   title('Chebyshev Polynomial Approximation')

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-May-2021 14:28:24

%% FUNCTION
%  [fun,x,domain] = InterpChebValues(values,x,domain)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, domain = []; end
if nargin < 2, x = []; end

%% ALGORITHM

coeffs = ChebCoefficients(values);
[fun,x,domain] = ChebPolyValues(coeffs,x,domain);

end