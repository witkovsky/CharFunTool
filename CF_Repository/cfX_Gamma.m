function cf = cfX_Gamma(t,alpha,beta,coef,niid)
%% cfX_Gamma 
%  Characteristic function of the GAMMA distribution with the shape
%  parameter alpha > 0 and the rate parameter beta > 0.
%
%  cfX_Gamma is an ALIAS NAME of the more general function cf_Gamma,
%  used to evaluate the characteristic function of a linear combination of
%  independent GAMMA distributed random variables.
%
%  The characteristic function of the GAMMA distribution is defined by
%   cf(t) = (1 - i*t/beta)^(-alpha)
%
% SYNTAX:
%  cf = cfX_Gamma(t,alpha,beta)
%  cf = cfX_Gamma(t,alpha,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - the shape parameter alpha > 0. If empty, default value is alpha
%          = 1.   
%  beta  - the rate (1/scale) parameter beta > 0. If empty, default value
%          is beta = 1.   
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
%   https://en.wikipedia.org/wiki/Gamma_distribution.
%
% EXAMPLE 1: 
%  % CF of the Gamma distribution with alpha = 2, beta = 2
%  alpha = 2;
%  beta = 2;
%  t = linspace(-20,20,501);
%  cf = cfX_Gamma(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Gamma distribution with alpha = 2, beta = 2')
%
% EXAMPLE 2: 
%  % PDF/CDF of the Gamma distribution with alpha = 2, beta = 2
%  alpha = 2;
%  beta = 2;
%  x = linspace(0,5,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^14;
%  cf = @(t) cfX_Gamma(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3: 
%  % PDF/CDF of the compound Binomial-Gamma distribution
%  n = 25;  
%  p = 0.3;
%  alpha = 2;
%  beta = 2;
%  cfX = @(t)  cfX_Gamma(t,alpha,beta);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,25,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jun-2017 10:07:43
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end

cf = cf_Gamma(t,alpha,beta,coef,niid);

end