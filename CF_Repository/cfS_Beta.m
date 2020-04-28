function cf = cfS_Beta(t,theta,coef,niid)
%% cfS_Beta
%  Characteristic function of the zero-mean symmetric BETA distribution
%  defined on the interval (-1,1). 
%
%  cfS_Beta is an ALIAS of the more general function cf_BetaSymmetric, used
%  to evaluate the characteristic function of a linear combination of
%  independent BETA distributed random variables.
%
%  The characteristic function of X ~ BetaSymmetric(theta) is defined by
%   cf(t) = cf_BetaSymmetric(t,theta) 
%          = gamma(1/2+theta) * (t/2)^(1/2-theta) * besselj(theta-1/2,t).
%
% SYNTAX:
%  cf = cfS_Beta(t,theta);
%  cf = cfS_Beta(t,theta,coef,niid);
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  theta - the 'shape' parameter theta > 0. If empty, default
%          value is theta = 1. 
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
%  https://en.wikipedia.org/wiki/Beta_distribution
%
% EXAMPLE 1:
% % CF of the symmetric Beta distribution with theta = 3/2 on (-1,1)
%   theta = 3/2;
%   t = linspace(-50,50,501);
%   cf = cfS_Beta(t,theta);
%   figure; plot(t,cf),grid
%   title('CF of the symmetric Beta distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the the symmetric Beta distribution on (-1,1)
%   theta = 3/2;
%   cf = @(t) cfS_Beta(t,theta);
%   x = linspace(-1,1,101);
%   xRange = 2;
%   clear options
%   options.N = 2^8;
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
if nargin < 2, theta = []; end

cf = cf_BetaSymmetric(t,theta,coef,niid);

end