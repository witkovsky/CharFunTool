function cf = cf_LogWilksLambdaRV(t,p,m,n,coef,niid)
%% cf_LogWilksLambdaRV 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED WILK's LAMBDA random variables.
%  
%  That is, cf_LogWilksLambdaRV evaluates the characteristic function of a
%  random variable Y  = coef_1*W_1 +...+ coef_N*W_N, such that cf_Y(t) =
%  cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t), where cf_W_i(t) is CF of W_i =
%  log(Lambda_i), and each Lambda_i has  WILK's LAMBDA distribution,
%  Lambda_i ~ Lambda(p_i,m_i,n_i), for i = 1,...,N. 
%
%  In particular Lambda_i = det(E_i)/det(E_i + H_i), with  E_i and H _i 
%  being independent random matrices with Wishart distributions E_i ~
%  Wp_i(m_i,Sigma_i) and H_i ~ Wp_i(n_i,Sigma_i) with m_i >= p_i and
%  Sigma_i > 0 (an unknown positive definite symmetric covariance matrix),
%  for all i = 1,...,N. Each particular Lambda_i statistic can be used to
%  test specific null hypothesis (by measuring its effect expressed by the
%  matrix H_i). In this case, the null hypothesis is rejected for small
%  values of the observed statistic Lambda_i, or large -log(Lambda_i).   
%
%  The Wilks’ Lambda distribution of Lambda_i ~ Lambda(p_i,m_i,n_i), with
%  Lambda_i in (0,1), is defined by Lambda_i ~ Prod_{j=1}^p_i B_{i,j},
%  where B_{i,j} ~ Beta{(m_i+1-j)/2, n_i/2)}, i.e. B_{i,j} follow
%  independent Beta distributions for all i = 1,...,N and j = 1,...,p_i.   
%
% SYNTAX
%  cf = cf_LogWilksLambdaRV(t,p,m,n,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  p     - vector of the dimension parameters p = (p_1,...,p_N). If empty,
%          default value is p = (1,...,1).  
%  m     - vector of sample size parameters m = (m_1,...,m_N). If empty,
%          default value is m = (1,...,1). 
%  n     - vector of sample size parameters n = (n_1,...,n_N). If empty,
%          default value is n = (1,...,1). 
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically
%          distributed random variable. If empty, default value is n = 1.  
%
% EXAMPLE 1:
% % CF of log Wilks Lambda RV distribution Lambda(p,m,n)
%   p  = 5;
%   m  = 10;
%   n  = 3;
%   t  = linspace(-10,10,201);
%   cf = cf_LogWilksLambdaRV(t,p,m,n);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of log Wilks Lambda RV with Lambda(p,m,n), p=5, m=10, n=3')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of minus log Wilks Lambda RVs 
%   p    = [5 5 5];
%   m    = [10 15 20];
%   n    = [3 2 1];
%   coef = -[10 15 20]/45;
%   t    = linspace(-20,20,201);
%   cf   = cf_LogWilksLambdaRV(t,p,m,n,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a weighted linear combination of -log Wilks Lambda RVs')
%
% EXAMPLE 3:
% % PDF/CDF of minus log Wilks Lambda RV (p=5, m=10, n=3) from its CF
%   p    = 5;
%   m    = 10;
%   n    = 3;
%   coef = -1;
%   cf   = @(t) cf_LogWilksLambdaRV(t,p,m,n,coef);
%   x    = linspace(0,5)';
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 4:
% % Compare the the exact distribution with the Bartlett's approximation
% % The Bartlett's approximation (see e.g. Wikipedia) is given by:
% % ((p-n+1)/2 - m)*log(Lambda(p,m,n)) ~ chi^2_{n*p}
%   p    = 15;
%   m    = 30;
%   n    = 3;
%   coef = (p-n+1)/2 - m;
%   cf   = @(t) cf_LogWilksLambdaRV(t,p,m,n,coef);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%   x = result.x;
%   figure;plot(x,result.cdf,x,chi2cdf(x,n*p));grid
%   title('Exact CDF vs. the Bartlett approximation')
%   xlabel('transformed \Lambda(p,m,n)')
%   ylabel('CDF')
%   disp(prob)
%   disp(result.qf)
%   disp(chi2inv(prob,n*p))
%
%  WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Wilks%27s_lambda_distribution

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Jun-2017 17:18:39

%% ALGORITHM
% cf = cf_LogWilksLambdaRV(t,p,m,n,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, n = []; end
if nargin < 3, m = []; end
if nargin < 2, p = []; end

%%
if isempty(p), p = 1; end
if isempty(m), m = 1; end
if isempty(n), n = 1; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Check size of the parameters
[errorcode,p,m,n,coef] = distchck(4,p(:)',m(:)',n(:)',coef(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:length(coef)
    alpha = (m(i)+1-(1:p(i)))/2;
    beta  = n(i)/2;
    cf = cf .* cf_LogBetaRV(coef(i)*t,alpha,beta);
end
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