function cf = cfS_Gaussian(t,mu,sigma,coef,niid)
%% cfS_Gaussian
%  Characteristic function of a GAUSSIAN (standard normal) distribution.
%
%  cfS_Gaussian is an ALIAS of the more general function cf_Normal, used to
%  evaluate the characteristic function of a linear combination of
%  independent normally distributed random variables.
%
%  The characteristic function of the standard normally distributed random
%  variable, X ~ N(0,1), is defined by 
%   cf(t) = cfS_Gaussian(t) = exp(-t^2/2);
%
% SYNTAX
%  cf = cfS_Gaussian(t)
%  cf = cfS_Gaussian(t,mu,sigma,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the 'location' parameters mu in R. If empty, default
%          value is mu = 0.  
%  sigma - vector of the 'scale' parameters sigma > 0. If empty, default
%          value is sigma = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          Normal random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
% 
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Normal_distribution.
%
% EXAMPLE 1:
% % CF of the Gaussian distribution N(0,1)
%   t = linspace(-5,5,501);
%   cf = cfS_Gaussian(t);
%   figure; plot(t,cf),grid
%   title('CF of the Gaussian distribution N(0,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the Gaussian distribution N(0,1)
%   cf   = @(t) cfS_Gaussian(t);
%   x    = linspace(-4,4,101);
%   prob = [0.9 0.95 0.975 0.99];
%   clear options
%   options.N = 2^5;
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,x,prob,options)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

cf = cf_Normal(t,mu,sigma,coef,niid);

end