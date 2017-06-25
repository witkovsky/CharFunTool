function cf = cf_LogRV_MeansRatioW(t,n,alpha,weight,coef,niid)
%% cf_LogRV_MeansRatioW 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED WEIGHTED MEANS-RATIO random variables (RVs)
%  W_i = log(R_i), for i = 1,...,N, where each R_i = G_i/A_i is a ratio of
%  the weighted geometric mean G_i and the (unweighted) arithmetic mean
%  A_i of independent RVs X_{i,1},...,X_{i,n_i} with GAMMA distributions,
%  with the shape parameters alpha_{i,j} and the (common) rate parameter
%  beta_i for i = 1,...,N and j = 1,...,n_i. 
%  
%  That is, cf_LogRV_MeansRatioW evaluates the characteristic function of a
%  random variable Y = coef_1*W_1 +...+ coef_N*W_N, such that cf_Y(t) =
%  cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t), where cf_W_i(t) is CF of W_i
%  = log(R_i), where R_i is the ratio statistic of the weighted geometric
%  mean and the arithmetic mean of independent gamma distributed RVs, which
%  distribution depends on the shape parameter alpha_{i,j} and the vector 
%  of weights w_{i,j}, for i = 1,...,N and j = 1,...,n_i.
%
%  Here, for each fixed i = 1,...,N, the weighted geometric mean is defined
%  by G_i = (X_{i,1}^w_{i,1} *...* X_{i,n_i}^w_{i,n_i}), where the weights
%  w_{i,j}, j = 1,...,n_i, are such that w_{i,1} +...+ w_{i,n_i} = 1, and
%  the arithmetic mean is defined by A_i = (X_{i,1} +...+ X_{i,n_i})/n_i.
%  The random variables  X_{i,j} are mutually independent with X_{i,j} ~
%  Gamma(alpha_{i,j},beta_i) for all i = 1,...,N and j = 1,...,n_i, where
%  alpha_{i,j} are the shape parameters and beta_i is the (common) rate
%  parameter of the gamma distributions. 
%
%  Note that the weighted ratio random variables R_i are scale invariant,
%  so their distribution does not depend on the common rate (or scale)
%  parameter beta_i for each i = 1,...,N.
%
%  The distribution of the logarithm of the means ratio log(R_i) is defined
%  by its characteristic function, see e.g. Chao and Glaser (JASA 1978),
%  which is
%   cf_{log(R_i)}(t) = (n_i)^(1i*t) * ...
%   Gamma(sum(alpha_{i,j}))/Gamma(sum(alpha_{i,j})+1i*t)) * ...
%   Prod_{j=1}^n_i Gamma(alpha_{i,j}+1i*w_{i,j}*t)/Gamma(alpha_{i,j}),
%  for each i = 1,...,N.
%
% SYNTAX
%  cf = cf_LogRV_MeansRatioW(t,n,alpha,weights,coef,niid)
%
% INPUTS
%  t       - vector or array of real values, where the CF is evaluated.
%  n       - vector of sample size parameters n = (n_1,...,n_N). 
%  alpha   - an N-by-1 cell array of weights vectors with shape parameters
%            alpha{i} = [alpha_{i,1},...,alpha_{i,n_i}].  If empty,
%            default value is alpha{i} = [1,...,1].
%  weights - an N-by-1 cell array of weights vectors with weights{i} =
%            [w_{i,1},...,w_{i,n_i}], for i = 1,...,N. If empty, default
%            value is weights{i}  = [1/k_i,...,1/k_i].
%  coef    - vector of the coefficients of the linear combination of the
%            log-transformed random variables. If coef is scalar, it is
%            assumed that all coefficients are equal. If empty, default
%            value is coef = 1.
%  niid    - scalar convolution coeficient niid, such that Z = Y +...+ Y
%            is sum of niid random variables Y, where each Y = sum_{i=1}^N
%            coef_i * log(X_i) is independently and identically distributed
%            random variable. If empty, default value is niid = 1.
%
%  WIKIPEDIA
%   https://en.wikipedia.org/wiki/Bartlett%27s_test
%
% EXAMPLE 1:
% % CF of log MeansRatio RV with and n = 5 and alpha = 7/2
%   n     = 5;
%   alpha = 7/2;
%   t     = linspace(-100,100,201);
%   cf    = cf_LogRV_MeansRatioW(t,n,alpha);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of log MeansRatio RV with and n = 5 and alpha = 7/2')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of minus log MeansRatio RVs 
%   clear
%   n        = [5 7 10];
%   alpha{1} = 7/2;
%   alpha{2} = 10/2;
%   alpha{3} = 3/2;
%   weight   = [];
%   coef     = -1/3;
%   t        = linspace(-100,100,201);
%   cf       = cf_LogRV_MeansRatioW(t,n,alpha,weight,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a weighted linear combination of minus log MeansRatio RVs')
%
% EXAMPLE 3:
% % PDF/CDF of minus log MeansRatio RV, n = 5 and alpha = 7/2, from its CF
%   clear
%   n     = 5;
%   alpha = 7/2;
%   coef  = -1;
%   cf    = @(t) cf_LogRV_MeansRatioW(t,n,alpha,[],coef);
%   x = linspace(0,0.6)';
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 4:
% % PDF/CDF of minus log of the (weighted) means ratio RV, from its CF
%   clear
%   n         = 5;
%   alpha{1}  = [3 5 7 10 3]/2;
%   weight{1} = alpha{1}/sum(alpha{1});
%   coef      = -1;
%   cf        = @(t) cf_LogRV_MeansRatioW(t,n,alpha,weight,coef);
%   x = linspace(-0.15,1.15,200)';
%   prob = [0.9 0.95 0.99];
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%
% EXAMPLE 5:
% % Compare the exact distribution with the Bartlett's approximation
%   clear
%   k         = 15; % k normal populations with unequal sample sizes
%   df        = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]; % degrees of freedom
%   DF        = sum(df);
%   alpha{1}  = df/2;
%   weight{1} = alpha{1}/sum(alpha{1});
%   C         = DF * log(k * prod(weight{1}.^weight{1}));
%   C_B       = 1 + 1/(3*(k-1))*(sum(1./df) - 1/DF);
%   shift     = C/C_B;
%   coef      = -DF/C_B;
%   cf_R      = @(t) cf_LogRV_MeansRatioW(t,k,alpha,weight,coef);
%   cf        = @(t) exp(1i*t*shift) .* cf_R(t)
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%   x = result.x;
%   figure;plot(x,result.cdf,x,chi2cdf(x,k-1));grid
%   title('Exact CDF vs. the Bartlett approximation')
%   xlabel('corrected test statistic')
%   ylabel('CDF')
%   disp(prob)
%   disp(result.qf)
%   disp(chi2inv(prob,k-1))
%
%  REFERENCES
%  Glaser, R. E. (1976a). The ratio of the geometric mean to the arithmetic
%  mean for a random sample from a gamma distribution. Journal of the
%  American Statistical Association, 71(354), 480-487.   
%
%  Glaser, R. E. (1976b). Exact critical values for Bartlett's test for
%  homogeneity of variances. Journal of the American Statistical
%  Association, 71(354), 488-490.  
%
%  Chao, M. T., & Glaser, R. E. (1978). The exact distribution of
%  Bartlett's test statistic for homogeneity of variances with unequal
%  sample sizes. Journal of the American Statistical Association, 73(362),
%  422-426.   

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Jun-2017 17:18:39

%% ALGORITHM
% cf = cf_LogRV_MeansRatioW(t,k,alpha,weight,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 6);
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, weight = []; end
if nargin < 3, alpha = []; end

%% Check size of the parameters
N = length(n);

if isempty(coef)
    coef = ones(1,N); 
end

[errorcode,n,coef] = distchck(2,n,coef);
if errorcode > 0
    error(message('InputSizeMismatch'));
end

if isempty(niid), niid = 1; end

if isempty(alpha) 
    alpha = cell(N,1); 
    for i = 1:N
        alpha{i} = ones(1,n(i));
    end
end

if isempty(weight) 
    weight = cell(N,1); 
    for i = 1:N
        weight{i} = ones(1,n(i))/n(i);
    end
end

if ~iscell(alpha)
    alpha = num2cell(alpha);
end

if ~iscell(weight)
    weight = num2cell(weight);
end

for i = 1:N
    aux = ones(1,n(i));
    [errorcode,alpha{i},weight{i},aux] = distchck(3,alpha{i},weight{i},aux);
    if errorcode > 0
        error(message('InputSizeMismatch'));
    end
end

if isempty(niid), niid = 1; end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:N
    cf = cf .* cf_LogR(coef(i)*t,n(i),alpha{i},weight{i});
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
%% Function cf_LogR
function cf = cf_LogR(t,k,alpha,weight)
% Auxiliary function. Calculates the characteristic function of the random
% variable log(R), where R is a ratio of weighted geometric mean and the
% unweighted arithmetic mean of independent gamma random variables with
% their shape parameters given by the vector alpha.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Jun-2017 12:08:24

A = sum(alpha);
cf = (1i*t)*log(k) + GammaLog(A) - GammaLog(A+1i*t);
for j = 1:k
    cf = cf + GammaLog(alpha(j)+1i*t*weight(j)) - GammaLog(alpha(j));
end
cf = exp(cf);

end