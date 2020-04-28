function cf = cf_Rice(t,distance,sigma,coef,niid,tol)
%cf_Rice 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Rice distributed random variables, with the location
%  parameter distance >= 0 (which can be interpreted as the distance
%  between the reference point and the center of the bivariate
%  distribution), and the scale parameter sigma > 0.
%  
%  The Rice distribution is a continuous probability distribution for
%  positive-valued random variables. It is a special case of the
%  non-central chi distribution with two degrees of freedom. The Rice
%  distribution is also related to the non-central Rayleigh distribution.
%
%  cf_Rice evaluates the characteristic function cf(t) of Y = sum_{i=1}^N
%  coef_i * X_i, where X_i~ RiceNC(distance_i,sigma_i) are inedependent
%  Rice distributed RVs with the location parameters distance_i >= 0, and
%  the scale parameters sigma_i > 0, for i  = 1,...,N.
%
%  The characteristic function of X ~ RiceNC(distance,sigma) is defined by
%   cf_Rice(t) = cf_ChiNC(sigma*t,df=2,delta=distance/sigma), 
%  where by cf_ChiNC(t,df,delta) we denote the characteristic function of
%  the noncentral chi distribution with df degrees of freedom and the
%  non-centrality parameter delta. Hence, the characteristic function of Y
%  is  
%   cf(t) = Prod ( cf_Rice(coef_i*t,distance_i,sigma_i) )
%
% SYNTAX:
%  cf = cf_Rice(t,distance,sigma,coef,niid,tol)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Rice distributed
%          random variables, sigma > 0.  If sigma is scalar, it is assumed
%          that all scale parameters are equal. If empty, default value is
%          sigma = 1.
%  distance - vector of the location parameters, distance >= 0. If empty, default
%          value is distance = 0. Notice that each component of the parameter
%          distance can be interpreted as a distance from the origin, i.e. as
%          a square root of the sum of squared means, distance_i =
%          sqrt(mu_{i,1}^2 + mu_{i,2}^2), of the associated generating
%          input variables X_i = sqrt(Z_{i,1}^2 + Z_{i,2}^2), where Z_{i,j}
%          ~ N(mu_{i,j},sigma_{i,j}^2), j = 1,2.
%  coef  - vector of the coefficients of the linear combination of the
%          Rice distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1. 
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%  https://en.wikipedia.org/wiki/Rice_distribution  
%
% REMARK:
%  The probability density function of the Rice distribution with the
%  location parameter distance >= 0, and the scale parameter sigma > 0 is 
%  pdf(x) = x\sigma^2 * exp((-(x^2+distance^2)/(2*sigma^2)) ...
%           * besseli(0,(x*distance)/sigma^2),
%  where besseli(0,z) is the modified Bessel function of the first kind
%  with order zero. 
%
% EXAMPLE 1:
% % CF of the distribution of Rice RV with sigma = 3
%   sigma = 1;
%   distance = 5;
%   t     = linspace(-5,5,501);
%   cf    = cf_Rice(t,distance,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the non-central Rice RV')
%
% EXAMPLE 2:
% % PDF/CDF of the distribution of Rice RV with sigma = 3
%   sigma = 1;
%   distance = 5;
%   cf    = @(t) cf_Rice(t,distance,sigma);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a linear combination of independent Rice RVs
%   sigma = [1 2 3];
%   distance = [1 1 1];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_Rice(t,distance,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Rice RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent Rice RVs
%   sigma = [1 2 3];
%   distance = [1 1 1];
%   coef  = [1 1 1];
%   cf    = @(t) cf_Rice(t,distance,sigma,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,20,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_ChiNC, cf_RayleighNC

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 20:29:23
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_Rice(t,distance,sigma,coef,niid,tol)
%  This is equivalent with the CF of non-central Rayleigh distribution
%  cf = cf_RayleighNC(t,sigma,distance./sigma,coef,niid,tol);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, tol   = []; end
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 3, sigma = []; end
if nargin < 2, distance = []; end

df    = 2;
delta = distance./sigma;
coef  = sigma.*coef;
cf    = cf_ChiNC(t,df,delta,coef,niid,tol);

end