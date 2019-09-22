function cf = cf_LogRV_HalfNormal(t,sigma,coef,niid)
%cf_LogRV_HalfNormal 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED HALF-NORMAL random variables.     
%  
%  The Half-Normal distribution is a continuous probability distribution
%  for positive-valued random variables. It is a scaled chi distribution
%  with one degree of freedom.
%
%  cf_LogRV_HalfNormal evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * log(X_i), where X_i ~ HN(sigma_i) are independent
%  Half-Normal RVs with the scale parameters sigma_i > 0, for i = 1,...,N.
%
%  The characteristic function of log(X) with X ~ HN(sigma) is
%   cf_LogRV_HalfNormal(t) = exp(1i*t*log(sigma)) * cf_LogRV_Chi(t,df)
%  where cf_LogRV_Chi(t,df) denotes the characteristic function of the
%  log-transformed Chi distribution with df=1 degree of freedom. 
%  
%  Hence, the characteristic function of Y is 
%   cf(t) = Prod ( cf_LogRV_HalfNormal(coef_i*t,sigma_i) )
%
% SYNTAX:
%  cf = cf_LogRV_HalfNormal(t,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the log-transformed Half-Normal
%          random variables.  If sigma is scalar, it is assumed that all
%          scale parameters are equal. If empty, default value is sigma =
%          1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed Half-Normal random variables. If coef is scalar, it
%          is assumed that all coefficients are equal. If empty, default
%          value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Half-Normal_distribution
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% EXAMPLE 1:
% % CF of the distribution of log-transformed Half-Normal RV with sigma = 3
%   sigma = 3;
%   t     = linspace(-10,10,501);
%   cf    = cf_LogRV_HalfNormal(t,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the log-transformed Half-Normal RV with sigma = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent log-transformed Half-Normal RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   t     = linspace(-3,3,501);
%   cf    = cf_LogRV_HalfNormal(t,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a combination of log-transformed Half-Normal RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Half-Normal RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   cf    = @(t) cf_LogRV_HalfNormal(t,sigma,coef);
%   clear options
%   options.N = 2^10;
%   x = linspace(-10,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_HalfNormal, cf_LogRV_Chi

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 22-Sep-2019 23:24:55

%% ALGORITHM
%  cf = cf_LogRV_HalfNormal(t,sigma,coef,niid)

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

% CF of the linear combination of the log-transformed Half-Normal RVs
% (expressed by using cf_LogRV_Chi) 
df    = 1;
shift = sum(coef.*log(sigma));
cf    = exp(1i*t*shift) .* cf_LogRV_Chi(t,df,coef,niid);

end