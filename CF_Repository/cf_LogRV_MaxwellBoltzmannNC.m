function cf = cf_LogRV_MaxwellBoltzmannNC(t,sigma,delta,coef,niid,tol)
%cf_LogRV_MaxwellBoltzmannNC 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED NON-CENTRAL MAXWELL-BOLTZMANN distributed
%  random variables, with the scale parameters sigma_i > 0, and the
%  non-centrality parameters delta_i >= 0, for i  = 1,...,N.
%  
%  The Maxwell-Boltzmann distribution is a continuous probability
%  distribution for positive-valued random variables. It is a scaled chi
%  distribution with three degrees of freedom.
%
%  cf_LogRV_MaxwellBoltzmannNC evaluates the characteristic function cf(t)
%  of Y = sum_{i=1}^N coef_i * log(X_i), where X_i~
%  MaxwellBoltzmannNC(sigma_i,delta_i) are inedependent non-central
%  Maxwell-Boltzmann distributed RVs with the scale parameters sigma_i > 0,
%  and the non-centrality parameters delta_i >= 0, for i  = 1,...,N.
%
%  The characteristic function of log(X) with X ~ MaxwellBoltzmannNC(sigma,
%  delta) is 
%   cf_LogRV_MaxwellBoltzmannNC(t) = exp(1i*t*log(sigma)) ...
%                                    * cf_LogRV_ChiNC(t,df,delta)                          
%  where cf_LogRV_ChiNC(t,df,delta) denotes the characteristic function of
%  the log-transformed non-central Chi distribution, here with df = 3
%  degrees of freedom and the non-centrality parameter delta.
%
%  Hence, the characteristic function of Y is 
%   cf(t) = Prod ( cf_LogRV_MaxwellBoltzmannNC(coef_i*t,sigma_i,delta_i) )
%
% SYNTAX:
%  cf = cf_LogRV_MaxwellBoltzmannNC(t,sigma,delta,coef,niid,tol)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Maxwell-Boltzmann
%          distributed random variables.  If sigma is scalar, it is assumed
%          that all scale parameters are equal. If empty, default value is
%          sigma = 1.
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0. Notice that the noncentrality
%          parameter delta can be interpreted as a square root of the sum
%          of standardized squared means, delta = sqrt(sum_{i=1}^3
%          mu_i^2/sigma^2_i). 
%  coef  - vector of the coefficients of the linear combination of the
%          Maxwell-Boltzmann distributed random variables. If coef is
%          scalar, it is assumed that all coefficients are equal. If empty,
%          default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1. 
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution       
%
% EXAMPLE 1:
% % CF of the log-transformed non-central Maxwell-Boltzmann RV
%   sigma = 1;
%   delta = 5;
%   t     = linspace(-20,20,501);
%   cf    = cf_LogRV_MaxwellBoltzmannNC(t,sigma,delta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the log-transformed non-central Maxwell-Boltzmann RV')
%
% EXAMPLE 2:
% % PDF/CDF of the log-transformed non-central Maxwell-Boltzmann RV
%   sigma = 1;
%   delta = 5;
%   cf    = @(t) cf_LogRV_MaxwellBoltzmannNC(t,sigma,delta);
%   clear options
%   options.N = 2^10;
%   x = linspace(0.5,2.5,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a combinantion of log-transformed non-central MaxwellBoltzmann RV
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   t     = linspace(-5,5,501);
%   cf    = cf_LogRV_MaxwellBoltzmannNC(t,sigma,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a combinantion of log-transformed non-central M-B RV')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of log-transformed non-central M-B RV
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   cf    = @(t) cf_LogRV_MaxwellBoltzmannNC(t,sigma,delta,coef);
%   clear options
%   options.N = 2^10;
%   x = linspace(0,6,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% SEE ALSO: cf_LogRV_ChiNC, cf_LogRV_RayleighNC

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-Oct-2018 23:48:44

%% ALGORITHM
%  cf = cf_LogRV_MaxwellBoltzmannNC(t,sigma,delta,coef,niid,tol)

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

%% CF OF THE LOG-TRANSFORMED NON-CENTRAL Maxwell-Boltzmann RVs
%  (expressed by using cf_LogRV_ChiNC) 
%  Here, we assume delta = sqrt(sum(mu_i^2/sigma_i^2))
%  Alternatively, if delta = sqrt(sum(mu_i^2)), set delta = delta./sigma;

df    = 3;
shift = sum(coef.*log(sigma));
cf    = exp(1i*t*shift) .* cf_LogRV_ChiNC(t,df,delta,coef,niid,tol);

end