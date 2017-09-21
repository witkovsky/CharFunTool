function [x,w] = LegendrePoints(n,A,B)
% LegendrePoints Evaluates the n-th order nodes (zeros of Legendre
%   polynomial Pn(x)) and the weights for the Gauss-Legendre quadrature
%   rule on the interval [A,B]. 
%
% SYNTAX:
%   [x,w] = LegendrePoints(n,A,B)
%
% SEE also: https://en.wikipedia.org/wiki/Gaussian_quadrature. 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Sep-2017 16:20:17

%% FUNCTION
%  [x,w] = LegendrePoints(n,A,B)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, B = []; end
if nargin < 2, A = []; end

if isempty(A), A = -1; end
if isempty(B), B = 1; end

%% ALGORITHM

x = zeros(1, n);
w = zeros(1, n);
m = (n+1)/2;
d = B-A;

for i  = 1:m
    zr0 = cos(pi * (i-0.25) / (n+0.5));
    zr1 = zr0+1;
    while abs(zr0-zr1) > eps
        p1 = 1;
        p2 = 0;
        for jj = 1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj-1) * zr0 * p2 - (jj-1) * p3) / jj;
        end
        pp  = n * (zr0*p1-p2) / (zr0^2-1);
        zr1 = zr0;
        zr0 = zr1 - p1 / pp;
    end
    x(i) = zr0;
    x(n+1-i) = -zr0;
    w(i) = d / ((1-zr0^2) * pp^2);
    w(n+1-i) = w(i);
end

if A ~= -1 || B ~= 1
    x = (x+1)*(d/2) + A;
end
end