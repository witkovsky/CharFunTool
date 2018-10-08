function cf = cf_LogRV_MaxwellBoltzmann(t,sigma,coef,niid)
%cf_LogRV_MaxwellBoltzmann 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED MAXWELL-BOLTZMANN distributed random
%  variables.      
%  
%  The Maxwell-Boltzmann distribution is a continuous probability
%  distribution for positive-valued random variables. It is a scaled chi
%  distribution with three degrees of freedom.
%
%  cf_LogRV_MaxwellBoltzmann evaluates the characteristic function cf(t) of
%  Y = sum_{i=1}^N coef_i * log(X_i), where X_i ~ MaxwellBoltzmann(sigma_i)
%  are inedependent non-central Maxwell-Boltzmann distributed RVs with the scale
%  parameters sigma_i > 0, for i  = 1,...,N.
%
%  The characteristic function of log(X) with X ~ MaxwellBoltzmann(sigma)
%  is 
%   cf_LogRV_MaxwellBoltzmann(t) = cf_LogRV_Chi(sigma*t,df=3) ...
%                                 * exp(1i*t*log(sigma))
%  where cf_LogRV_Chi(t,df) denotes the characteristic function of the
%  log-transformed Chi distribution with df degrees of freedom.
%
%  Hence, the characteristic function of Y is 
%   cf(t) = Prod ( cf_LogRV_MaxwellBoltzmann(coef_i*t,sigma_i) )
%
% SYNTAX:
%  cf = cf_LogRV_MaxwellBoltzmann(t,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the log-transformed
%          Maxwell-Boltzmann distributed random variables.  If sigma is
%          scalar, it is assumed that all scale parameters are equal. If
%          empty, default value is sigma = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed Maxwell-Boltzmann distributed random variables.
%          If coef is scalar, it is assumed that all coefficients are
%          equal. If empty, default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
%  https://en.wikipedia.org/wiki/Chi_distribution   
%
% EXAMPLE 1:
% % CF of the distribution of Maxwell-Boltzmann RV with sigma = 3
%   sigma = 3;
%   t     = linspace(-10,10,501);
%   cf    = cf_LogRV_MaxwellBoltzmann(t,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the Maxwell-Boltzmann RV with sigma = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent Maxwell-Boltzmann RVs
%   sigma = [1 2 3];
%   coef  = [1 1 1];
%   t     = linspace(-5,5,501);
%   cf    = cf_LogRV_MaxwellBoltzmann(t,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Maxwell-Boltzmann RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Maxwell-Boltzmann RVs
%   sigma = [1 2 3];
%   coef  = [1 1 1];
%   cf    = @(t) cf_LogRV_MaxwellBoltzmann(t,sigma,coef);
%   clear options
%   options.N = 2^10;
%   x = linspace(-2,6,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% SEE ALSO: cf_MaxwellBoltzmann, cf_LogRV_Chi, cf_LogRV_Rayleigh 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 19:10:00

%% ALGORITHM
%  cf = cf_LogRV_MaxwellBoltzmann(t,sigma,coef,niid)

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

% CF of the linear combination of the log-transformed Maxwell-Boltzmann RVs
% (expressed by using cf_LogRV_Chi) 
df    = 3;
shift = sum(coef.*log(sigma));
cf    = exp(1i*t*shift) .* cf_LogRV_Chi(t,df,coef,niid);

end