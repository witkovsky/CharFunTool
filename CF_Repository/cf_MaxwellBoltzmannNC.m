function cf = cf_MaxwellBoltzmannNC(t,sigma,delta,coef,niid)
%cf_MaxwellBoltzmannNC 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Maxwell-Boltzmann distributed random variables,
%  with the scale parameters sigma_i > 0, and the non-centrality parameters
%  delta_i >= 0, for i  = 1,...,N.      
%  
%  The non-central Maxwell-Boltzmann distribution is a continuous
%  probability distribution for positive-valued random variables. It is a
%  special case of the non-central chi distribution in three degrees of
%  freedom.  
%
%  cf_MaxwellBoltzmannNC evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i~ MaxwellBoltzmannNC(sigma_i,delta_i)
%  are inedependent non-central Maxwell-Boltzmann distributed RVs with the
%  scale parameters sigma_i > 0,  and the non-centrality parameters
%  delta_i >= 0, for i  = 1,...,N.      
%
%  The characteristic function of X ~ MaxwellBoltzmannNC(sigma,delta) is
%  defined by 
%   cf_MaxwellBoltzmannNC(t) = cf_ChiNC(sigma*t,df=3,delta/sigma), 
%  where by cf_ChiNC(t,df,delta) we denote the characteristic function of
%  the noncentral chi distribution with df degrees of freedom and the
%  non-centrality parameter delta. Hence, the characteristic function of Y
%  is  
%   cf(t) = Prod ( cf_MaxwellBoltzmannNC(t,sigma_i,delta_i) )
%
% SYNTAX:
%  cf = cf_MaxwellBoltzmannNC(t,sigma,delta,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Maxwell-Boltzmann
%          distributed random variables.  If sigma is scalar, it is assumed
%          that all scale parameters are equal. If empty, default value is
%          sigma = 1.
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0.  Notice that each component of the
%          non-centrality parameter delta can be interpreted as a square
%          root of the sum of squared standardized means, delta_i =
%          sqrt((mu_{i,1}/sigma_{i,1})^2 + (mu_{i,2}/sigma_{i,2})^2 +
%          (mu_{i,3}/sigma_{i,3})^2), of the associated generating input
%          variables X_i = sqrt(Z_{i,1}^2 + Z_{i,2}^2 + Z_{i,3}^2), where
%          Z_{i,j} ~ N(mu_{i,j},sigma_{i,j}^2), j = 1,2,3.
%  coef  - vector of the coefficients of the linear combination of the
%          Maxwell-Boltzmann distributed random variables. If coef is
%          scalar, it is assumed that all coefficients are equal. If empty,
%          default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%  https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
%
% NOTES (from Wikipedia)
%  Mathematically, the Maxwell–Boltzmann distribution with sigma = 1 is the
%  chi distribution with three degrees of freedom (the components of the
%  velocity vector in Euclidean space).     
%
% EXAMPLE 1:
% % CF of the distribution of Maxwell-Boltzmann RV with sigma = 3
%   sigma = 1;
%   delta = 5;
%   t     = linspace(-5,5,501);
%   cf    = cf_MaxwellBoltzmannNC(t,sigma,delta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the non-central Maxwell-Boltzmann RV')
%
% EXAMPLE 2:
% % PDF/CDF of the distribution of Maxwell-Boltzmann RV with sigma = 3
%   sigma = 1;
%   delta = 5;
%   cf    = @(t) cf_MaxwellBoltzmannNC(t,sigma,delta);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a linear combination of independent Maxwell-Boltzmann RVs
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_MaxwellBoltzmannNC(t,sigma,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Maxwell-Boltzmann RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent Maxwell-Boltzmann RVs
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   cf    = @(t) cf_MaxwellBoltzmannNC(t,sigma,delta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,25,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_ChiNC, cf_MaxwellBoltzmann 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 5-Oct-2018 17:08:51

%% ALGORITHM
%  cf = cf_MaxwellBoltzmannNC(t,sigma,delta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
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

%% CF OF THE LINEAR COMBINATION OF THE NON-CENTRAL MAXWELL-BOLTZMANN RVs

%  REMARK
%  Notice that there are two possibilities how to set the noncentrality
%  parameter delta:
%  1. Here, by default we assume that delta = sqrt(sum(mu_i^2/sigma_i^2)),
%     which is normalized with respect to the scale parameter sigma.
%  2. Alternatively, we can assume that delta represents fixed distance
%     from the origin, which does not depend on the scale parameter sigma,
%     i.e. delta = sqrt(sum(mu_i^2)). In this case use delta./sigma as the
%     non-centrality parameter input:
%     cf = cf_MaxwellBoltzmannNC(t,sigma,delta./sigma,coef,niid)

df   = 3;
coef = sigma.*coef;
cf   = cf_ChiNC(t,df,delta,coef,niid);


end