function cf = cf_LogRV_Weibull(t,alpha,beta,coef,niid)
%%cf_LogRV_Weibull 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED WEIBULL random variables (RVs) log(X), where
%  X ~ Weibull(alpha,beta), and alpha > 0 and beta > 0 represent the scale
%  and the shape parameters of the WEIBULL distribution.
%  
%  That is, cf_LogRV_Weibull evaluates the characteristic function cf(t) of
%  Y = coef_i*log(X_1) +...+ coef_N*log(X_N), where X_i ~
%  Weibull(alpha_i,beta_i), with the parameters alpha_i > 0 and beta_i > 0,
%  for i = 1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ Weibull(alpha,beta)
%  is defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) =
%  E(X^(1i*t)). That is, the characteristic function can be derived from
%  expression for the r-th moment of X, E(X^r) by using (1i*t) instead of
%  r. In particular, the characteristic function of Y = log(X) is
%   cf_Y(t) = alpha^(1i*t) .* gamma(1 + 1i*t/beta).
%
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters alpha(i) and beta(i).
%
% SYNTAX
%  cf = cf_LogRV_Weibull(t,alpha,beta,coef,niid)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  alpha  - vector of the 'scale' parameters alpha > 0. If empty, default
%           value is alpha = 1.  
%  beta   - vector of the 'shape' parameters beta > 0. If empty, default
%           value is beta = 1.  
%  coef   - vector of the coefficients of the linear combination of the
%           LOG-TRANSFORMED WEIBULL random variables. If coef is scalar, it
%           is assumed that all coefficients are equal. If empty, default
%           value is coef = 1.
%  niid   - scalar convolution coeficient niid, such that Z = Y +...+ Y is
%           sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%           coef(i) * log(X_i) is independently and identically distributed
%           random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Weibull_distribution.
%
% REMARK:
%  The probability density function of a Weibull random variable is
%    pdf(x) = beta/alpha * (x/alpha)^{beta-1} * exp(-(x/alpha)^beta),
%  for x >=0, where beta > 0 is the shape parameter and alpha > 0 is the
%  scale parameter of the distribution. The characteristic function of
%  log(X) is given by  
%   cf(t) = alpha^(1i*t) .* gamma(1 + 1i*t/beta).
% 
% EXAMPLE 1:
% % CF of a linear combination of independent log-transformed Weibull RVs
%   coef   = [1 2 3];
%   alpha  = 5/2;
%   beta   = 3/2;
%   t      = linspace(-2,2,1001);
%   cf     = cf_LogRV_Weibull(t,alpha,beta,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a combination of independent log-transformed Weibull RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a linear combination of independent log-Gamma RVs
%   coef   = [1 2 3];
%   alpha  = 5/2;
%   beta   = 3/2;
%   cf     = @(t) cf_LogRV_Weibull(t,alpha,beta,coef);
%   x      = linspace(-10,15);
%   prob   = [.9,.95,.99];
%   clear options
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 7-Oct-2018 14:31:53

%% ALGORITHM
% cf = cf_LogRV_Weibull(t,alpha,beta,coef,n)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end

%%
if isempty(beta)
    beta = 1;
end

if isempty(alpha) 
    alpha = 1;
end

if isempty(coef) 
    coef = 1;
end

if isempty(niid)
    niid = 1;
end

%% Check size of the parameters
[errorcode,coef,alpha,beta] = distchck(3,coef(:)',alpha(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of linear combination 
szt = size(t);
t   = t(:);
cf  = prod(exp(GammaLog(1+1i*t*(coef./beta)) + 1i*t*(coef.*log(alpha))),2);
cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end