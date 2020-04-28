function cf = cf_FoldedNormal(t,mu,sigma,coef,niid)
%cf_FoldedNormal 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Folded-Normal random variables.     
%   
%  The Folded-Normal distribution is a continuous probability distribution
%  for positive-valued random variables which is proportional to the
%  non-central chi-distribution with one degree of freedom. 
% 
%  In particular, if U ~ N(mu/sigma,1) then X = sigma*|U| has Folded-Normal
%  distribution with the parameters mu and sigma, X ~ FN(mu,sigma). That
%  is, the random variable X is multiple of the random variable |U| which
%  has non-central chi-distribution with one degree of freedom and the
%  noncentrality parameter delta = mu/sigma.   
%
%  cf_FoldedNormal evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Folded-Normal(mi_i,sigma_i) are
%  inedependent Folded-Normal RVs with the location parameters mu_i and the
%  scale parameters sigma_i > 0, for i  = 1,...,N. 
%
%  The characteristic function of X ~ FN(mu,sigma) is defined by
%   cf_FoldedNormal(t) = cf_ChiNC(sigma*t,df=1,delta=mu/sigma), 
%  where cf_ChiNC(t,df,delta) denotes the characteristic function of the
%  non-central Chi-distribution with df degrees of freedom and the
%  non-centrality parameter delta. Hence, the characteristic 
%  function of Y is 
%   cf(t) = Prod ( cf_FoldedNormal(coef_i*t,sigma_i) )
%         = Prod ( cf_FoldedNormal(coef_i*sigma_i*t,1) )
%
%  Alternatively, see Tsagris et al. (2014), the characteristic function of
%  X ~ FN(mu,sigma) can be expressed by 
%   cf_FoldedNormal(t) = exp(-sigma^2*t^2/2 + 1i*mu*t) 
%                        * (1/2)*(1-erf((-mu/sigma+1i*sigma*t)/sqrt(2)))
%                      + exp(-sigma^2*t^2/2 - 1i*mu*t) 
%                        * (1/2)*(1-erf((mu/sigma+1i*sigma*t)/sqrt(2)))
%  where erf(z) denotes the complex error function. See also erfZX.m
%
% SYNTAX:
%  cf = cf_FoldedNormal(t,mu,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the location parameters of the Folded-Normal random
%          variables. If mu is scalar, it is assumed that all location
%          parameters are equal. If empty, default value is mu = 0. 
%  sigma - vector of the scale parameters of the Folded-Normal random
%          variables. If sigma is scalar, it is assumed that all scale
%          parameters are equal. If empty, default value is sigma = 1. 
%  coef  - vector of the coefficients of the linear combination of the
%          Folded-Normal random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Folded_normal_distribution
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%
% NOTES (from Wikipedia)
%  The folded normal distribution is a probability distribution related to
%  the normal distribution. Given a normally distributed random variable X
%  with mean ? and variance ?2, the random variable Y = |X| has a folded
%  normal distribution. Such a case may be encountered if only the
%  magnitude of some variable is recorded, but not its sign. The
%  distribution is called "folded" because probability mass to the left of
%  the x = 0 is folded over by taking the absolute value. In the physics of
%  heat conduction, the folded normal distribution is a fundamental
%  solution of the heat equation on the upper plane (i.e. a heat kernel).
%
%  - When mu = 0, the distribution X ~ FN(0,sigma) is a half-normal
%    distribution, i.e. X ~ HN(sigma).
%  - The random variable (X/?)^2 has a non-central chi-squared distribution
%    with 1 degree of freedom and noncentrality equal to (mu/?)^2.  
%  - The folded normal distribution can also be seen as the limit of the
%    folded non-standardized t-distribution as the degrees of freedom go to
%    infinity.   
%
% EXAMPLE 1:
% % CF of the distribution of Folded-Normal RV with mu = 1 and sigma = 3
%   mu    = 1;
%   sigma = 3;
%   t     = linspace(-5,5,501);
%   cf    = cf_FoldedNormal(t,mu,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the Folded-Normal RV with mu = 1 and sigma = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent Folded-Normal RVs
%   mu    = [1 1 1 1 1];
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_FoldedNormal(t,mu,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Folded-Normal RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Folded-Normal RVs
%   mu    = [1 1 1 1 1];
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   cf    = @(t) cf_FoldedNormal(t,mu,sigma,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,30,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_ChiNC, erfZX
%
% REFERENCES:
%  [1] Tsagris M., Beneki C. and Hassani H. (2014). On the folded normal
%  distribution. Mathematics, 2(1), pp.12-28. 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Sep-2019 12:33:09
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_FoldedNormal(t,mu,sigma,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

if isempty(mu)
    mu = 0;
end

if isempty(sigma)
    sigma = 1;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and sigma 
[errorcode,coef,sigma,mu] = distchck(3,coef(:),sigma(:),mu(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

% CF of the linear combination of the Folded-Normal RVs (by using ChiNC) 
df = 1;
delta = mu./sigma;
cf = cf_ChiNC(t,df,delta,sigma.*coef,niid);

end