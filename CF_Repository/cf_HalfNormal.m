function cf = cf_HalfNormal(t,sigma,coef,niid)
%cf_HalfNormal 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Half-Normal random variables.     
%   
%  The Half-Normal distribution is a continuous probability distribution
%  for positive-valued random variables which is proportional to the
%  chi-distribution with one degree of freedom. 
% 
%  In particular, if Z ~ N(0,1) then X = sigma*|Z| has Half-Normal
%  distribution with scale parameter sigma, X ~ HN(sigma). That is, the
%  random variable X is multiple of the random variable |Z| which has
%  chi-distribution with one degree of freedom.  
%
%  Note that MATHEMATICA uses different parametrization of the Half-Normal
%  distribution: X ~ HalfNormalDistribution[theta], where theta =
%  sqrt(pi/2)/sigma, and otherwise, sigma = sqrt(pi/2)/theta. 
%
%  cf_HalfNormal evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Half-Normal(sigma_i) are
%  inedependent Half-Normal RVs with the scale parameters sigma_i > 0, for
%  i  = 1,...,N. 
%
%  The characteristic function of X ~ HN(sigma) is defined by
%   cf_HalfNormal(t) = cf_Chi(sigma*t,df=1), 
%  where cf_Chi(t,df) denotes the characteristic function of the Chi
%  distribution with df degrees of freedom. Hence, the characteristic
%  function of Y is 
%   cf(t) = Prod ( cf_HalfNormal(coef_i*t,sigma_i) )
%         = Prod ( cf_HalfNormal(coef_i*sigma_i*t,1) )
%
%  Alternatively, the characteristic function of X ~ HN(sigma) can  be
%  expressed by 
%   cf_HalfNormal(t) = exp(-sigma^2*t^2/2) * (1-1i*erfi(sigma*t/sqrt(2))),
%          = exp(-sigma^2*t^2/2) + (2i/sqrt(pi)) * Dawson(sigma*t/sqrt(2))
%  where erf(z) denotes the complex error function and Dawson(z) denotes
%  the Dawson function. Note that the second expression is more suitable
%  for practical purposesis as it is numerically more stable for large t.
%  See also erfiZX.m and Dawson.m.  
%
% SYNTAX:
%  cf = cf_HalfNormal(t,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Half-Normal random
%          variables. If sigma is scalar, it is assumed that all scale
%          parameters are equal. If empty, default value is sigma = 1. 
%  coef  - vector of the coefficients of the linear combination of the
%          Half-Normal random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Half-normal_distribution
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% NOTES (from Wikipedia)
%  - The Half-Normal distribution is a special case of the Folded-Normal
%    distribution with its location parameter mu = 0. 
%  - The Half-Normal distribution is a special case of the generalized
%    gamma distribution with d = 1, p = 2, a = sqrt(2)*sigma.  
%  - The Half-Normal distribution coincides with a zero-mean normal
%    distribution truncated from below at zero (see truncated normal
%    distribution).
%  - If Y has a Half-Normal distribution, then (Y/sigma)^2 has a chi-square
%    distribution with 1 degree of freedom, i.e. Y/sigma has a
%    chi-distribution with 1 degree of freedom.  
%
% EXAMPLE 1:
% % CF of the distribution of Half-Normal RV with sigma = 3
%   sigma = 3;
%   t     = linspace(-5,5,501);
%   cf    = cf_HalfNormal(t,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the Half-Normal RV with sigma = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent Half-Normal RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   t     = linspace(-1,1,501);
%   cf    = cf_HalfNormal(t,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Half-Normal RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Half-Normal RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   cf    = @(t) cf_HalfNormal(t,sigma,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(5,40,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_Chi, erfZX, Dawson

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Sep-2019 12:33:09
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_HalfNormal(t,sigma,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid  = []; end
if nargin < 3, coef  = []; end
if nargin < 2, sigma = []; end

if isempty(sigma)
    sigma = 1;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and sigma 
[errorcode,coef,sigma] = distchck(2,coef(:),sigma(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

% CF of the linear combination of the Half-Normal RVs (expressed by using Chi)
df = 1;
cf = cf_Chi(t,df,sigma.*coef,niid);

end