function cf = cfX_Exponential(t,lambda,coef,niid)
%% cfX_Exponential
%  Characteristic function of the EXPONENTIAL distribution with the rate
%  parameter lambda > 0.
%
%  cfX_Exponential is an ALIAS NAME of the more general function
%  cf_Exponential, used to evaluate the characteristic function of a linear
%  combination of independent EXPONENTIAL distributed random variables.
%
%  The characteristic function of the EXPONENTIAL distribution is
%   cf(t) = lambda / (lambda - 1i*t)
%
% SYNTAX:
%  cf = cfX_Exponential(t,lambda)
%  cf = cfX_Exponential(t,lambda,coef,niid)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  lambda - vector of the 'rate' parameters lambda > 0. If empty, default
%           value is lambda = 1.   
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
%  https://en.wikipedia.org/wiki/Exponential_distribution.
%
% EXAMPLE 1:
%  % CF of the Exponential distribution with lambda = 5
%  lambda = 5;  
%  t = linspace(-50,50,501);
%  cf = cfX_Exponential(t,lambda);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Exponential distribution with lambda = 5')
%
% EXAMPLE 2:
%  % PDF/CDF of the Exponential distribution with lambda = 5
%  lambda = 5;  
%  cf = @(t) cfX_Exponential(t,lambda);
%  x  = linspace(0,1.5,101);
%  clear options
%  options.xMin = 0;
%  options.SixSigmaRule = 8;
%  result = cf2DistGP(cf,x,[],options)
%
% EXAMPLE 3:
% % PDF/CDF of the compound Binomial-Exponential distribution
%  n = 25;  
%  p = 0.3;
%  lambda = 5;
%  cfX  = @(t) cfX_Exponential(t,lambda);
%  cf   = @(t) cfN_Binomial(t,n,p,cfX);
%  x    = linspace(0,5,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jun-2017 10:07:43
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 4);
if nargin < 4, niid = []; end
if nargin < 3, coef = []; end
if nargin < 2, lambda = []; end

cf = cf_Exponential(t,lambda,coef,niid);

end