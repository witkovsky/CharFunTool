function cf = cfS_Trapezoidal(t,lambda,coef,niid)
%% cfS_Trapezoidal
%  Characteristic function of the zero-mean symmetric TRAPEZOIDAL
%  distribution defined on the interval (-1,1). 
%
%  cfS_Trapezoidal is an ALIAS of the more  general function
%  cf_TrapezoidalSymmetric, used to evaluate the characteristic function of
%  a linear combination of independent TRAPEZOIDAL distributed random
%  variables.
%
%  The characteristic function of X ~ TrapezoidalSymmetric(lambda), where
%  0<= lambda <=1 is the offset parameter is defined by
%   cf(t) = (sin(w*t)./(w*t)).*(sin((1-w)*t)./((1-w)*t)), 
%  where w  = (1+lambda)/2.
%
% SYNTAX:
%  cf = cfS_Trapezoidal(t,lambda);
%  cf = cfS_Trapezoidal(t,lambda,coef,niid);
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  lambda - parameter of the offset, 0 <= lambda <=1. If empty, default
%           value is lambda = 0.
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1. 
%
% EXAMPLE 1:
% % CF of the symmetric Trapezoidal distribution with lambda = 0.5
%   lambda = 0.5;
%   t = linspace(-50,50,501);
%   cf = cfS_Trapezoidal(t,lambda);
%   figure; plot(t,cf),grid
%   title('CF of the symmetric Trapezoidal distribution with lambda = 0.5')
%
% EXAMPLE 2:
% % PDF/CDF of the symmetric Trapezoidal distribution, lambda = 0.5
%   lambda = 0.5;
%   cf = @(t) cfS_Trapezoidal(t,lambda)
%   x = linspace(-1,1,101);
%   xRange = 2;
%   clear options
%   options.N = 2^10;
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
narginchk(1, 4);
if nargin < 4, niid  = []; end
if nargin < 3, coef  = []; end
if nargin < 2, lambda = []; end

cf = cf_TrapezoidalSymmetric(t,lambda,coef,niid);

end