function cf = cf_HalfNormalNC(t,sigma,delta,coef,niid,tol)
%cf_HalfNormalNC 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Half-Normal distributed random variables, with
%  the scale parameters sigma_i > 0, and the non-centrality parameters
%  delta_i >= 0, for i  = 1,...,N. 
%  
%  The non-central Half-Normal distribution is a continuous probability
%  distribution for positive-valued random variables. It is a special case
%  of the non-central chi distribution with one degree of freedom.
% 
%  When delta = mu, the non-central Half-Normal distribution is known
%  as the Folded-Normal distribution with parameters mu and sigma.
%
%  cf_HalfNormalNC evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i~ HalfNormalNC(sigma_i,delta_i) are
%  independent non-central Half-Normal distributed RVs with the scale
%  parameters sigma_i > 0, and the non-centrality parameters delta_i >= 0,
%  for i = 1,...,N.
%
%  The characteristic function of X ~ HalfNormalNC(sigma,delta) is defined
%  by 
%   cf_HalfNormalNC(t) = cf_ChiNC(sigma*t,df=1,delta/sigma), 
%  where by cf_ChiNC(t,df,delta) we denote the characteristic function of
%  the noncentral chi distribution with df degrees of freedom and the
%  non-centrality parameter delta. Hence, the characteristic function of Y
%  is  
%   cf(t) = Prod ( cf_HalfNormalNC(t,sigma_i,delta_i) )
%
% SYNTAX:
%  cf = cf_HalfNormalNC(t,sigma,delta,coef,niid,tol)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Half-Normal distributed
%          random variables.  If sigma is scalar, it is assumed that all
%          scale parameters are equal. If empty, default value is sigma =
%          1.
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0.  Notice that each component of the
%          non-centrality parameter delta can be interpreted as a square
%          root of the sum of squared standardized means, delta_i =
%          sqrt((mu_{i,1}/sigma_{i,1})^2 + (mu_{i,2}/sigma_{i,2})^2), of
%          the associated generating input variables X_i = sqrt(Z_{i,1}^2 +
%          Z_{i,2}^2), where Z_{i,j} ~ N(mu_{i,j},sigma_{i,j}^2), j = 1,2.
%  coef  - vector of the coefficients of the linear combination of the
%          Half-Normal distributed random variables. If coef is scalar, it
%          is assumed that all coefficients are equal. If empty, default
%          value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%  https://en.wikipedia.org/wiki/Half-normal_distribution
%
% NOTES 
%  Mathematically, the Half-Normal distribution with sigma = 1 is the
%  chi-distribution with one degree of freedom.     
%
% EXAMPLE 1:
% % CF of the distribution of non-central Half-Normal RV with sigma = 3
%   sigma = 1;
%   delta = 5;
%   t     = linspace(-5,5,501);
%   cf    = cf_HalfNormalNC(t,sigma,delta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the non-central Half-Normal RV')
%
% EXAMPLE 2:
% % PDF/CDF of the non-central Half-Normal RV with sigma = 3
%   sigma = 1;
%   delta = 5;
%   cf    = @(t) cf_HalfNormalNC(t,sigma,delta);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a linear combination of independent non-central Half-Normal RVs
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_HalfNormalNC(t,sigma,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Half-Normal RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent NC Half-Normal RVs
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   cf    = @(t) cf_HalfNormalNC(t,sigma,delta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,20,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 5:
% % PDF/CDF of a linear combination of independent NC Half-Normal RVs
%   sigma = [1 2 3];
%   delta = [1 1 1]./sigma;
%   coef  = [1 1 1];
%   cf    = @(t) cf_HalfNormalNC(t,sigma,delta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,20,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_ChiNC, cf_HalfNormal

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 22-Sep-2019 23:08:23
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_HalfNormalNC(t,sigma,delta,coef,niid,tol)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, tol   = []; end
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 3, delta = []; end
if nargin < 2, sigma = []; end

if isempty(sigma)
    sigma = 1;
end

if isempty(delta)
    delta = 0;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and sigma 
[errorcode,coef,sigma,delta] = distchck(3,coef(:),sigma(:),delta(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

%% CF OF THE LINEAR COMBINATION OF THE NON-CENTRAL HALF-NORMAL RVs
%  Here, we assume delta = sqrt(sum(mu_i^2/sigma_i^2))
%  Alternatively, if delta = sqrt(sum(mu_i^2)), set delta = delta./sigma;

df   = 1;
coef = sigma.*coef;
cf   = cf_ChiNC(t,df,delta,coef,niid,tol);

end