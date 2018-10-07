function cf = cf_LogRV_Rayleigh(t,sigma,coef,niid)
%cf_LogRV_Rayleigh 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED RAYLEIGH random variables.     
%  
%  The Rayleigh distribution is a continuous probability distribution for
%  positive-valued random variables. It is a scaled chi distribution with
%  two degrees of freedom.
%
%  cf_LogRV_Rayleigh evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * log(X_i), where X_i ~ Rayleigh(sigma_i) are
%  independent Rayleigh RVs with the scale parameters sigma_i > 0, for i
%  = 1,...,N. 
%
%  The characteristic function of log(X) with X ~ Rayleigh(sigma) is
%   cf_LogRV_Rayleigh(t) = exp(1i*t*log(sigma)) * cf_LogRV_Chi(t,df)
%  where cf_LogRV_Chi(t,df) denotes the characteristic function of the
%  log-transformed Chi distribution with df=2 degrees of freedom. 
%  
%  Hence, the characteristic function of Y is 
%   cf(t) = Prod ( cf_LogRV_Rayleigh(coef_i*t,sigma_i) )
%
% SYNTAX:
%  cf = cf_LogRV_Rayleigh(t,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the log-transformed Rayleigh
%          random variables.  If sigma is scalar, it is assumed that all
%          scale parameters are equal. If empty, default value is sigma =
%          1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed Rayleigh random variables. If coef is scalar, it
%          is assumed that all coefficients are equal. If empty, default
%          value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Rayleigh_distribution
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% EXAMPLE 1:
% % CF of the distribution of log-transformed Rayleigh RV with sigma = 3
%   sigma = 3;
%   t     = linspace(-10,10,501);
%   cf    = cf_LogRV_Rayleigh(t,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the log-transformed Rayleigh RV with sigma = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent log-transformed Rayleigh RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   t     = linspace(-3,3,501);
%   cf    = cf_LogRV_Rayleigh(t,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a combination of log-transformed Rayleigh RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Rayleigh RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   cf    = @(t) cf_LogRV_Rayleigh(t,sigma,coef);
%   clear options
%   options.N = 2^10;
%   x = linspace(-2,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_Rayleigh, cf_LogRV_Chi, cf_LogRV_MaxwellBoltzmann

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 18:07:49

%% ALGORITHM
%  cf = cf_LogRV_Rayleigh(t,sigma,coef,niid)

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

% CF of the linear combination of the log-transformed Rayleigh RVs
% (expressed by using cf_LogRV_Chi) 
df    = 2;
shift = sum(coef.*log(sigma));
cf    = exp(1i*t*shift) .* cf_LogRV_Chi(t,df,coef,niid);

end