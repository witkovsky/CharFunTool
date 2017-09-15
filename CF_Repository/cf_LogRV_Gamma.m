function cf = cf_LogRV_Gamma(t,alpha,beta,coef,niid)
%%cf_LogRV_Gamma 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED GAMMA random variables (RVs) log(X), where X
%  ~ Gamma(alpha,beta), and alpha > 0 and beta > 0 represent the shape and
%  the rate parameters of the GAMMA distribution.
%  
%  That is, cf_LogRV_Gamma evaluates the characteristic function cf(t) of Y
%  = coef_i*log(X_1) +...+ coef_N*log(X_N), where X_i ~
%  Gamma(alpha_i,beta_i), with the parameters alpha_i > 0 and beta_i > 0,
%  for i = 1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ Gamma(alpha,beta) is
%  defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)). 
%  That is, the characteristic function can be derived from expression for
%  the r-th moment of X, E(X^r) by using (1i*t) instead of r. In
%  particular, the characteristic function of Y = log(X) is
%   cf_Y(t) = (1/beta^(1i*t)) .* (gamma(alpha + 1i*t) / gamma(alpha)).
%
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters alpha(i) and beta(i).
%
% SYNTAX
%  cf = cf_LogRV_Gamma(t,alpha,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'rate' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          logGamma random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y +...+ Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Gamma_distribution
%
% EXAMPLE 1:
% % CF of a weighted linear combination of independent log-Gamma RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   alpha  = 5/2;
%   beta   = 3/2;
%   t      = linspace(-20,20,1001);
%   cf     = cf_LogRV_Gamma(t,alpha,beta,-weight);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a linear combination of minus log-Gamma RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a linear combination of independent log-Gamma RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   alpha  = 5/2;
%   beta   = 3/2;
%   cf     = @(t) cf_LogRV_Gamma(t,alpha,beta,weight);
%   clear options
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,[],[],options);
%   disp(result)
%
% EXAMPLE 3:
% % PDF/CDF of a log-Gamma (log-ChiSquared) RV with alpha = 5/2, beta = 1/2 
%   alpha = 5/2;
%   beta = 1/2;
%   cf = @(t) cf_LogRV_Gamma(t,alpha,beta);
%   clear options
%   options.N = 2^12;
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%   qfGamma = exp(result.qf)
%
% REFERENCES:
%   MATHAI, A.M. (1973). A review of the different techniques used for
%   deriving the exact distributions of multivariate test criteria.
%   Sankhya: The Indian Journal of Statistics, Series A, 39-60.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24

%% ALGORITHM
% cf = cf_LogRV_Gamma(t,alpha,beta,coef,n)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end

%%
if isempty(beta) && ~isempty(alpha)
    beta = 1;
elseif isempty(beta) && ~isempty(coef)
    beta = 1;
end

if isempty(alpha) && ~isempty(coef)
    alpha = 1;
elseif isempty(alpha) && ~isempty(beta)
    alpha = 1;
end

if isempty(coef) && ~isempty(beta)
    coef = 1;
elseif isempty(coef) && ~isempty(alpha)
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
aux = 1i*t*coef;
aux = GammaLog(bsxfun(@plus,aux,alpha))-ones(length(t),1)*GammaLog(alpha);
aux = aux - 1i*t*log(beta);
cf  = prod(exp(aux),2);
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