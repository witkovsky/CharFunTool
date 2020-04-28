function cf = cfS_Rectangular(t,coef,niid)
%% cfS_Rectangular
%  Characteristic function of the zero-mean symmetric RECTANGULAR
%  distribution defined on the interval (-1,1). 
%
%  cfS_Rectangular is an ALIAS of the more general function
%  cf_RectangularSymmetric, used to evaluate the characteristic function of
%  a linear combination of independent RECTANGULAR distributed random
%  variables.
%
%  The characteristic function of X ~ RectangularSymmetric is defined by
%   cf(t) = sinc(t) = sin(t)/t;
%
% SYNTAX:
%  cf = cfS_Rectangular(t);
%  cf = cfS_Rectangular(t,coef,niid);
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1. 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
%
% EXAMPLE 1:
% % CF of the Rectangular distribution on (-1,1)
%   t = linspace(-50,50,501);
%   cf = cfS_Rectangular(t);
%   figure; plot(t,cf),grid
%   title('CF of the Rectangular distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the Rectangular distribution on (-1,1)
%   cf = @(t) cfS_Rectangular(t);
%   x = linspace(-2,2,101);
%   xRange = 2;
%   clear options
%   options.N = 2^5;
%   options.dt = 2*pi/xRange;
%   result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
%   WITKOVSKY V. (2016). Numerical inversion of a characteristic
%   function: An alternative tool to form the probability distribution of
%   output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 3);
if nargin < 3, niid = []; end
if nargin < 2, coef = []; end

cf = cf_RectangularSymmetric(t,coef,niid);

end