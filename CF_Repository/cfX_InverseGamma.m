function cf = cfX_InverseGamma(t,alpha,beta,coef,niid)
%% cfX_InverseGamma
%  Characteristic function of the INVERSE GAMMA distribution with the shape
%  parameter alpha > 0 and the rate parameter beta > 0.
%
%  cfX_InverseGamma is an ALIAS NAME of the more general function
%  cf_InverseGamma, used to evaluate the characteristic function of a
%  linear combination of independent INVERSE GAMMA distributed random
%  variables.
%
%  The characteristic function of the GAMMA distribution is defined by
%   cf(t) =  2 / gamma(alpha) * (-1i*beta*t).^(alpha/2) ...
%                 * besselk(alpha,sqrt(-4i*beta*t)).
%
% SYNTAX:
%  cf = cfX_InverseGamma(t,alpha,beta)
%  cf = cfX_InverseGamma(t,alpha,beta,coef,niid)
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
%  https://en.wikipedia.org/wiki/Inverse-gamma_distribution
%
% EXAMPLE 1:
%  % CF of the Inverse Gamma distribution with alpha = 2, beta = 2
%  alpha = 2;
%  beta = 2;
%  t = linspace(-20,20,501);
%  cf = cfX_InverseGamma(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the InverseGamma distribution with alpha = 2, beta = 2')
%
% EXAMPLE 2:
%  % PDF/CDF of the Inverse Gamma distribution with alpha = 2, beta = 2
%  alpha = 2;
%  beta = 2;
%  x = linspace(0,15,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.SixSigmaRule = 10;
%  options.N = 2^14;
%  cf = @(t) cfX_InverseGamma(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3:
%  % PDF/CDF of the compound Binomial-InverseGamma distribution
%  n = 25;  
%  p = 0.3;
%  alpha = 2;
%  beta = 2;
%  cfX = @(t)  cfX_InverseGamma(t,alpha,beta);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,70,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
%  WITKOVSKY, V.: Computing the distribution of a linear combination of
%  inverted gamma variables, Kybernetika 37 (2001), 79-90.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end

cf = cf_InverseGamma(t,alpha,beta,coef,niid);

end