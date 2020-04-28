function cf = cf_GeneralizedExponential(t,alpha,lambda,mu,coef,niid)
%% cf_GeneralizedExponential 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent GENERALIZED EXPONENTIAL random variables (RVs) where X ~
%  GE(alpha,lambda,mu). Here alpha > 0 is a shape parameter, lambda > 0 is
%  a scale parameter and mu is a location parameter.  
%  
%  We say that the random variable X has a generalized exponential
%  distribution, X ~ GE(alpha,lambda,mu), if X has the cumulative
%  distribution function  
%   F(x;alpha,lambda,mu) = (1 - e^{-(x-mu)/lambda})^alpha 
%  for x > mu, alpha > 0, lambda > 0, with the corresponding probability
%  density function  
%   f(x;alpha,lambda,mu) = (alpha/lambda) * ...
%                          (1 - e^{-(x-mu)/lambda})^(alpha-1) * ...
%                          e^{-(x-mu)/lambda}.
%  In particular if X ~ GE(alpha) = GE(alpha,1,0) then mu + lambda * X ~
%  GE(alpha,lambda,mu). 
%  The characteristic function of X ~  GE(alpha,lambda,mu) is given by
%   cf_X(t) = exp(1i*t*mu) * gamma(alpha+1) * gamma(1 - 1i*t*lambda) ...
%             / gamma(alpha - 1i*t*lambda +1)
%  For more details see e.g. Gupta and Kundu (1999).
%
%  That is, cf_GeneralizedExponential evaluates the characteristic function
%  cf(t) of  Y =  coef_i*X_1 +...+ coef_N*X_N, where X_i ~
%  GE(alpha_i,lambda_i,mu_i) are inedependent RVs, with the parameters
%  alpha_i > 0, lambda_i > 0, and real mu_i, for i = 1,...,N. Hence,the
%  characteristic function of Y  = coef(1)*X1 + ... + coef(N)*XN  is
%  cf_Y(t) =  cf_X1(coef(1)*t) * ... * cf_XN(coef(N)*t), where cf_Xi(t) 
%  is evaluated with the parameters alpha(i), lambda(i) and mu(i).
%
% SYNTAX
%  cf = cf_GeneralizedExponential(t,alpha,lambda,mu,coef,niid)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  alpha  - vector of the 'shape' parameters alpha > 0. If empty, default
%           value is alpha = 1.  
%  lambda - vector of the 'scale' parameters lambda > 0. If empty, default
%           value is alpha = 1. 
%  mu     - vector of the 'location' parameters, mu is real. If empty,
%           default value is mu = 0.  
%  coef   - vector of the coefficients of the linear combination of the
%           random variables. If coef is scalar, it is assumed that all
%           coefficients are equal. If empty, default value is coef = 1.
%  niid   - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%           sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%           coef(i) * X_i is independently and identically distributed
%           random variable. If empty, default value is niid = 1.  
%
% EXAMPLE 1:
% % CF of a linear combination of generalized exponential RVs
%   alpha  = [1 2 3 4 5];
%   lambda = [1 1 2 2 3];
%   mu     = [0 0 1 1 1];
%   coef   = [1 2 3 4 5]/15;
%   t      = linspace(-5,5,201);
%   cf     = cf_GeneralizedExponential(t,alpha,lambda,mu,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('Characteristic function of a linear combination of GE RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a linear combination of GE RVs from the CF by cf2DistGP
%   alpha  = [1 2 3 4 5];
%   lambda = [1 1 2 2 3];
%   mu     = [0 0 1 1 1];
%   coef   = [1 2 3 4 5]/15;
%   cf     = @(t) cf_GeneralizedExponential(t,alpha,lambda,mu,coef);
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options);
%
% REFERENCES:
%   [1] Gupta, R.D. and Kundu, D., 1999. Generalized exponential
%       distributions. Australian & New Zealand Journal of Statistics,
%       41(2), 173-188.   
%   [2] Gupta, R.D. and Kundu, D., 2007. Generalized exponential
%       distribution: Existing results and some recent developments.
%       Journal of Statistical Planning and Inference, 137(11), 3537-3547.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 14-Sep-2018 13:36:38
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf_GeneralizedExponential(t,alpha,lambda,mu,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid   = []; end
if nargin < 5, coef   = []; end
if nargin < 4, mu     = []; end
if nargin < 3, lambda = []; end
if nargin < 2, alpha  = []; end

%%
if isempty(alpha)
    alpha = 1;
end

if isempty(lambda)
    lambda = 1;
end

if isempty(mu)
    mu = 0;
end

if isempty(coef)
    coef = 1;
end

%% Check size of the parameters
[errorcode,coef,alpha,lambda,mu] = ...
    distchck(4,coef(:)',alpha(:)',lambda(:)',mu(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);
cf  = GammaLog(1+(ones(size(t))*alpha));
cf  = cf + GammaLog(1-1i*bsxfun(@times,t,coef.*lambda));
cf  = cf - GammaLog(1+(ones(size(t))*alpha)-1i*bsxfun(@times,t,coef.*lambda));
cf  = cf + 1i*bsxfun(@times,t,coef.*mu);
cf  = prod(exp(cf),2);
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