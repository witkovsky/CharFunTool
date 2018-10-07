function cf = cf_LogRV_RayleighNC(t,sigma,delta,coef,niid)
%cf_LogRV_RayleighNC 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED NON-CENTAL RAYLEIGH distributed random
%  variables, with the scale parameters sigma_i > 0, and the non-centrality
%  parameters delta_i >= 0, for i  = 1,...,N.
%  
%  The non-central Rayleigh distribution is a continuous probability
%  distribution for positive-valued random variables. It is a scaled
%  non-central chi distribution in three degrees of freedom. 
%
%  cf_LogRV_RayleighNC evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * log(X_i), where X_i~ RayleighNC(sigma_i,delta_i)
%  are inedependent non-central Rayleigh distributed RVs with the scale
%  parameters sigma_i > 0,  and the non-centrality parameters delta_i >= 0,
%  for i  = 1,...,N.
%
%  The characteristic function of log(X) with X ~ RayleighNC(sigma,delta)
%  is 
%   cf_LogRV_RayleighNC(t) = cf_LogRV_ChiNC(sigma*t,delta,df=2) ...
%                          * exp(1i*t*log(sigma))
%  where cf_LogRV_ChiNC(t,df,delta) denotes the characteristic function of
%  the log-transformed non-central Chi distribution with df degrees of
%  freedom and the noncentrality parameter delta.
%  
%  Hence, the characteristic function of Y is 
%   cf(t) = Prod ( cf_LogRV_RayleighNC(coef_i*t,sigma_i,delta_i) )
%
% SYNTAX:
%  cf = cf_LogRV_RayleighNC(t,sigma,delta,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Rayleigh distributed
%          random variables.  If sigma is scalar, it is assumed that all
%          scale parameters are equal. If empty, default value is sigma =
%          1.
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0. Notice that the noncentrality
%          parameter delta can be interpreted as a square root of the sum
%          of squared means, delta = sqrt(sum_{i=1}^2 mu_i^2).
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed non-central Rayleigh distributed random
%          variables. If coef is scalar, it is assumed that all
%          coefficients are equal. If empty, default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Rayleigh_distribution
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution  
%
% EXAMPLE 1:
% % CF of the log-transformed non-central Rayleigh RV
%   sigma = 1;
%   delta = 5;
%   t     = linspace(-5,5,501);
%   cf    = cf_LogRV_RayleighNC(t,sigma,delta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the log-transformed non-central Rayleigh RV')
%
% EXAMPLE 2:
% % PDF/CDF of the log-transformed non-central Rayleigh RV
%   sigma = 1;
%   delta = 5;
%   cf    = @(t) cf_LogRV_RayleighNC(t,sigma,delta);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a combination of the log-transformed non-central Rayleigh RV
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_LogRV_RayleighNC(t,sigma,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Rayleigh RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a combination of the log-transformed non-central Rayleigh RV
%   sigma = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   cf    = @(t) cf_LogRV_RayleighNC(t,sigma,delta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,20,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% SEE ALSO: cf_LogRV_ChiNC, cf_LogRV_Rayleigh

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 18:36:25

%% ALGORITHM
%  cf = cf_LogRV_RayleighNC(t,sigma,delta,coef,niid)

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

%% CF OF THE LINEAR COMBINATION OF THE NON-CENTRAL RAYLEIGH RVs
%  Here, we assume delta = sqrt(sum(mu_i^2/sigma_i^2))
%  Alternatively, if delta = sqrt(sum(mu_i^2)), set delta = delta./sigma;

df   = 2;
coef = sigma.*coef;
cf   = cf_ChiNC(t,df,delta,coef,niid);

end