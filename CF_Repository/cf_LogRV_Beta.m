function cf = cf_LogRV_Beta(t,alpha,beta,coef,niid)
%% cf_LogRV_Beta 
%  Characteristic function of a linear combination (resp.
%  convolution) of independent LOG-TRANSFORMED BETA random variables (RVs)
%  log(X), where X ~ BETA(alpha,beta), and alpha > 0 and beta > 0
%  represent the shape parameters of the BETA distribution. 
%  
%  That is, cf_Beta evaluates the characteristic function cf(t) of  Y =
%  coef_i*log(X_1) +...+ coef_N*log(X_N), where X_i ~ Beta(alpha_i,beta_i)
%  are inedependent RVs, with the shape parameters alpha_i > 0 and beta_i
%  >0, for i = 1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ Beta(alpha,beta) is
%  defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)). 
%  That is, the characteristic function can be derived from expression for
%  the r-th moment of X, E(X^r) by using (1i*t) instead of r. In
%  particular, the characteristic function of Y = log(X), with X ~
%  Beta(alpha,beta) is defined by  
%   cf_Y(t) = gamma(alpha + 1i*t) / gamma(alpha) .* ...
%             gamma(alpha + beta) / gamma(alpha + beta + 1i*t).
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters alpha(i) and beta(i).
%
% SYNTAX
%  cf = cf_LogRV_Beta(t,alpha,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'shape' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
%
%  WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Beta_distribution
%
% EXAMPLE 1:
% % CF of a linear combination of K=50 independent log-Beta RVs
%   coef = 1./(((1:50) - 0.5)*pi).^2;
%   figure; plot(coef,'.-'); grid on;
%   title('Coefficients of the linear combination of log-Beta RVs')
%   alpha = 5/2;
%   beta  = 3/2;
%   t     = linspace(-100,100,201);
%   cf    = cf_LogRV_Beta(t,alpha,beta,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('Characteristic function of a linear combination of log-Beta RVs')
%
% EXAMPLE 2:
% % PDF/CDF from the CF by cf2DistGP
%   alpha = 5/2;
%   beta  = 3/2;
%   coef  = 1./(((1:50) - 0.5)*pi).^2;
%   cf    = @(t) cf_LogRV_Beta(t,alpha,beta,coef);
%   clear options
%   options.xMax = 0;
%   result = cf2DistGP(cf,[],[],options);
%
% EXAMPLE 3: 
% % Distribution of log(R), where R=geometric/arithmetic mean of Gamma RVs
% % Let X_1,...,X_n are iid RVs, X_j ~ Gamma(A,B), where A > 0 is the known
% % shape parameter and B > 0 is the (unknown, common) rate parameter.  
% % Let R = geometricmean(X)/mean(X). According to Glaser (JASA 1976) 
% % log(R) ~ (1/n) * sum_{j=1}^{n-1} log(Y_j), Y_j ~ Beta(alpha,beta_j),
% % where alpha = A and beta_j = j/n for j = 1,...,n-1. That is, log(R) is
% % distributed as linear combination of independent logBeta random
% % variables log(Y_j). 
% % Here we evaluate the PDF/CDF of W = -log(R) (i.e. minus of log(R))
%   n = 10;
%   alpha = 1; % A = 1, i.e. X_j are from exponential distribution
%   beta  = (1:n-1)'/n;
%   coef  = -1/n;
%   cf = @(t) cf_LogRV_Beta(t,alpha,beta,coef);
%   t = linspace(-25,25,201);
%   figure; plot(t,real(cf(t)),t,imag(cf(t))); grid on;
%   prob = [ 0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%
% EXAMPLE 4 (Distribution of Wilks Lambda statistic) 
% % If E ~ Wp(m,Sigma) and H ~ Wp (n,Sigma) with m >= p, then the
% % distribution of Lambda = L = det(E)/det(E + H) is Wilks Lambda
% % distribution, denoted by L ~ Lambda(p,m,n), with L in (0,1). 
% % It holds that L ~ Prod_{i=1}^p B_i, where B_i follow independent
% % Beta distributions with Bi ~ B{(m + 1 - i)/2, n/2)}, for i = 1,...,p. 
% % Let W = -log(L) and cdf_W(x) = Prob(W <= x) = Prob(-sum(log(Bi)) <= x). 
% % Then, cdf_L(u) = Prob(L <= u) = Prob(W > -log(u)) = 
% % 1 - Prob(W <= -log(u)) = 1 - cdf_W(x), where x = -log(u) and u =
% % exp(-x). Moreover, pdf_Lambda(u) = pdf_W(x)/exp(-x) = pdf_W(-log(u))/u. 
% % The Lambda statistic is used to test null (significance) hypothesis
% % expressed by the matrix H. The null hypothesis is rejected for small
% % values of the observed statistic Lambda, or large value -log(Lambda).
%   p = 10;
%   m = 20;
%   n = 7;
%   i  = 1:p;
%   alpha = (m+1-i)/2;
%   beta  = n/2;
%   coef  = -1;
%   cf = @(t) cf_LogRV_Beta(t,alpha,beta,coef);
%   t = linspace(-5,5,201);
%   figure; plot(t,real(cf(t)),t,imag(cf(t))); grid on;
%   prob = [ 0.99 0.95 0.9];
%   x = linspace(2,8);
%   clear options
%   options.xMin = 0; 
%   % Distribution of -log(Lambda)
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%   % PDF of Lambda
%   % If W = -log(Lambda), then pdf_Lambda(u) = pdf_W(x)/exp(-x), 
%   % where x = -log(u) and u = exp(-x)
%   figure; plot(exp(-result.x),result.pdf./exp(-result.x))
%   xlabel('\Lambda')
%   ylabel('PDF')
%   title('PDF of Wilks \Lambda distribution with p=10, m = 20, n=7')
%   % CDF of Lambda
%   % If W = -log(Lambda), then cdf_Lambda(u) = 1-cdf_W(x), 
%   % where x = -log(u) and u = exp(-x)
%   figure; plot(exp(-result.x),1-result.cdf)
%   xlabel('\Lambda')
%   ylabel('CDF')
%   title('CDF of Wilks \Lambda distribution with p=10, m = 20, n=7')
%   disp('Quantiles of Wilks Lambda distribution')
%   disp(['alpha    = ',num2str(1-prob)])
%   disp(['quantile = ',num2str(exp(-result.qf))])
%
% REFERENCES:
%   GLASER, R.E. (1976). The ratio of the geometric mean to the arithmetic
%   mean for a random sample from a gamma distribution. Journal of the
%   American Statistical Association, 71(354), 480-487.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24

%% ALGORITHM
% cf = cf_LogRV_Beta(t,alpha,beta,coef,niid)

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

%% Check size of the parameters
[errorcode,coef,alpha,beta] = distchck(3,coef(:)',alpha(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);
aux = 1i*bsxfun(@times,t,coef);
aux = GammaLog(bsxfun(@plus,aux,alpha)) - ...
      GammaLog(bsxfun(@plus,aux,alpha+beta));
aux = bsxfun(@plus,aux,ones(length(t),1) * ...
      (GammaLog(alpha+beta)-GammaLog(alpha)));
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