function cf = cf_LogRV_WilksLambdaNC(t,p,m,n,delta,coef,niid,MAX)
%% cf_LogRV_WilksLambdaNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED WILK's LAMBDA random variables.
%
%  That is, cf_LogRV_WilksLambdaNC evaluates the characteristic function of a
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
%  where B_{i,j} ~ Beta{(m_i+1-j)/2,n_i/2)}, i.e. B_{i,j} follow
%  independent Beta distributions for all i = 1,...,N and j = 1,...,p_i.
%
% SYNTAX
%  cf = cf_LogRV_WilksLambdaNC(t,p,q,n,delta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  p     - vector of the dimension parameters p = (p_1,...,p_N). If empty,
%          default value is p = (1,...,1).
%  m     - vector of degrees of freedom of the Wishart matrices E_i, m =
%          (m_1,...,m_N). If empty, default value is m = (1,...,1).
%  n     - vector of degrees of freedom of the Wishart matrices H_i, n =
%          (n_1,...,n_N). If empty, default value is n = (1,...,1).
%  delta - p.s.d. matrix or vector of nonnegative eigenvalues or cell array
%          of matrices or vectors of eigenvalues of the non-centrality
%          parameters, delta = {delta_1,...,delta_N). Default value is
%          delta = [].
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically
%          distributed random variable. If empty, default value is n = 1.
%  MAX   - maximum number of partitions used for computing the
%          hypergeometric 1F1 function with matrix argument, for more
%          details see HypergeomMatFun.m. If empty MAX, default value is
%          MAX = 10;
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Wilks%27s_lambda_distribution
%
% EXAMPLE 1:
% % CF of log of noncentral Wilks Lambda RV distribution Lambda(p,m,n,delta)
%   p  = 5;
%   m  = 10; % elsewhere it is denoted as n (d.f. of within SS&P)
%   n  = 3;  % elsewhere it is denoted as q (d.f. of between SS&P)
%   delta = sort(rand(1,p));
%   t  = linspace(-10,10,201);
%   cf = cf_LogRV_WilksLambdaNC(t,p,m,n,delta);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of log of noncentral Wilks Lambda RV')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of minus log Wilks Lambda RVs 
%   p     = [5 5 5];
%   m     = [10 15 20];
%   n     = [3 2 1];
%   delta = {sort(rand(1,p(1))),sort(rand(1,p(2))),sort(rand(1,p(3)))};
%   coef  = -[10 15 20]/45;
%   t     = linspace(-20,20,201);
%   cf    = cf_LogRV_WilksLambdaNC(t,p,m,n,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a weighted linear combination of -log Wilks Lambda RVs')
%
% EXAMPLE 3:
% % PDF/CDF of minus log Wilks Lambda RV (p=5, m=10, n=3) from its CF
%   p     = 5;
%   m     = 10;
%   n     = 3;
%   delta = sort(rand(1,p));
%   coef  = -1;
%   cf    = @(t) cf_LogRV_WilksLambdaNC(t,p,m,n,delta,coef);
%   x     = linspace(0,5)';
%   prob  = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 4:
% % Compare the exact distribution with the Bartlett's approximation
% % The Bartlett's approximation (see e.g. Wikipedia) is given by:
% % ((p-n+1)/2 - m)*log(Lambda(p,m,n)) ~ chi^2_{n*p}
%   p     = 15;
%   m     = 30;
%   n     = 3;
%   delta = [];
%   coef  = (p-n+1)/2 - m;
%   cf    = @(t) cf_LogRV_WilksLambdaNC(t,p,m,n,delta,coef);
%   prob  = [0.9 0.95 0.99];
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
% REMARK:
%  cf_LogRV_WilksLambdaNC is FRAGILE and could lead to nonstable results.
%  Computing CF of the LOG-TRANSFORMED NON-CENTARL WILK's LAMBDA random
%  variable depends on computing the generalized hypergeometric function of
%  matrix argument. Current version of cf_LogRV_WilksLambdaNC uses modified
%  implementation of the algorithm suggested in Koev and Edelman (2006). 
%
% REFERENCES:
% [1] Koev, P. and Edelman, A., 2006. The efficient evaluation of the
%     hypergeometric function of a matrix argument. Mathematics of
%     Computation, 75(254), 833-846.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Oct-2017 12:44:48

%% ALGORITHM
% cf = cf_LogRV_WilksLambdaNC(t,p,m,n,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 8);
if nargin < 8, MAX   = []; end
if nargin < 7, niid  = []; end
if nargin < 6, coef  = []; end
if nargin < 5, delta = []; end
if nargin < 4, n = []; end
if nargin < 3, m = []; end
if nargin < 2, p = []; end

%%
if isempty(p), p = 1; end
if isempty(m), m = 1; end
if isempty(n), n = 1; end
if isempty(delta), delta = cell(1); end
if isempty(coef),  coef  = 1; end
if isempty(niid),  niid  = 1; end

if ~iscell(delta)
    [p1,p2] =size(delta);
    if p1 == p2
        delta{1} = eig((delta+delta')/2);
    elseif min(p1,p2) == 1
        delta = mat2cell(delta,1);
    elseif min(p1,p2) ~= 1
        error(message('InputSizeMismatch'));
    end
else
    for i = 1:length(delta)
        [p1,p2] =size(delta{i});
        if p1 == p2 && p1 > 0
            delta{i} = eig((delta{i}+delta{i}')/2);
        elseif min(p1,p2) > 1
            error(message('InputSizeMismatch'));
        end
    end
end

%% Check size of the parameters
[errorcode,p,m,n,coef] = distchck(4,p(:)',m(:)',n(:)',coef(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

if length(delta) == 1
    delta = repmat(delta,length(coef));
end

%% Characteristic function of a linear combination
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:length(coef)
    alpha = (m(i)+1-(1:p(i)))/2;
    beta  = n(i)/2;
    cf = cf .* cf_LogRV_Beta(coef(i)*t,alpha,beta);
    if ~isempty(delta{i})
        cf = cf .* Hypergeom1F1Mat(1i*coef(i)*t, ...
            1i*coef(i)*t + (m(i)+n(i))/2,-delta{i},MAX);
    end
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
%% FUNCTION Hypergeom1F1Mat
function f = Hypergeom1F1Mat(a,b,x,MAX)
% Hypergeom1F1Mat Computes truncated hypergeometric function
% 1F1^alpha(a;b;x) with parameters a and b and matrix argument x. the
% maximum number of partitions in the truncated hypergeometric series is
% equal to given parameter MAX (default value is MAX = 10).

%% ALGORITHM
if nargin < 4, MAX = 10; end

alpha = 1;
y     = [];
nt    = length(a);
f     = zeros(nt,1);

for i = 1:nt
    f(i) = HypergeompFqMat(a(i),b(i),x,y,alpha,MAX);
end
end