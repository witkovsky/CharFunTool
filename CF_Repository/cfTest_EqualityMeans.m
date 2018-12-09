function cf = cfTest_EqualityMeans(t,n,p,q,type)
%% cfTest_EqualityMeans 
%  Characteristic function of the exact null distribution of the
%  multivariate analysis (MVA) test statistic for testing equality of
%  mean vectors of q normal p-dimensional populations (q>1). For more
%  details see Anderson (2003) and/or Marques and Coelho (2011). 
%
%  In particular, here we consider minus log-transformed LRT statistic
%  (Likelihood Ratio Test) or alternatively the minus log-transformed
%  modified LRT statistic. 
%
%  Let X_k ~ N_p(mu_k,Sigma) are p-dimensional random vectors with common
%  covariance matrix Sigma for all k = 1,...,q. We want to test the
%  hypothesis that the mean vectors mu_k are equal for all X_k, k =
%  1,...,q. Then, the null hypothesis is given as
%    H0: mu_1 = ... = mu_q, 
%  i.e. the mean vectors are equal in all q populations. Here, the LRT test
%  statistic is given by  
%    LRT = Lambda^(n/2) = ( det(E) / det(E+H) )^{n/2},
%  and the modified LRT is defined as
%    MLRT = Lambda = det(E) / det(E+H)
%  where Lambda = det(E) / det(E+H) with E = sum_{k=1}^q sum_{j=1}^{n_k}
%  (X_{kj} - bar{X}_k)'*(X_{kj} - bar{X}_k) with E ~ Wishart(n-q,Sigma),
%  and H = sum_{k=1}^q (bar{X}_k - bar{X})'*(bar{X}_k - bar{X}) with H ~
%  Wishart(q-1,Sigma) based on n = n_1 + ... + n_q  samples from the q
%  p-dimensional populations.
%
%  Under null hypothesis, distribution of the LRT statistic is
%    LRT ~  prod_{j=1}^{p} (B_j)^{n/2}, 
%  with B_j ~ Beta((n-q-j+1)/2,(q-1)/2). Here we assume n > p+q-1. 
% 
% SYNTAX
%   cf = cfTest_EqualityMeans(t,n,p,q,type)
%
% INPUTS
%  t    - vector or array of real values, where the CF is evaluated.
%  n    - total number of samples, n = n_1 + ... + n_q from the q
%         p-dimensional populations. If all sampes sizes are equal to n_1,
%         then n = n_1*q .It is assumed that n > p+q-1.  
%  p    - dimension of each of the q populations. It is assumed that each
%         population has the same (equal) dimension p.
%  q    - number of populations, q>1. If empty or missing, the default
%         value of q is q = 2 (p is scalar). 
%  type - specifies the type of the LRT test statistic: type = 'standard'
%         with W = -(n/2)*log(Lambda), or type = 'modified' with W =
%         -log(Lambda). If type is empty or missing, default value is type
%         = 'modified'.
%
% EXAMPLE 1:
% % CF of the log-transformed LRT statistic for equal means
%   t  = linspace(-2,2,201);
%   n  = 12; % say n = 3+4+5
%   p  = 5;
%   q  = 3;
%   type = 'standard';
%   cf = cfTest_EqualityMeans(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing equal means')
%
% EXAMPLE 2:
% % CF of the log-transformed modified LRT statistic for equal means
%   t  = linspace(-10,10,201);
%   n  = 12; % say n = 3+4+5
%   p  = 5;
%   q  = 3;
%   type = 'modified';
%   cf = cfTest_EqualityMeans(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing equal means')
%
% EXAMPLE 3:
% % PDF/CDFF/QF of the log-transformed LRT for equal covariances
%   n  = 12; % say n = 3+4+5
%   p    = 5;
%   q    = 3;
%   type = 'standard';
%   cf   = @(t) cfTest_EqualityMeans(t,n,p,q,type);
%   x    = linspace(0,30,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 4:
% % PDF/CDFF/QF of the log-transformed modified LRT for equal means
%   n  = 12; % say n = 3+4+5
%   p    = 5;
%   q    = 3;
%   type = 'modified';
%   cf   = @(t) cfTest_EqualityMeans(t,n,p,q,type);
%   x    = linspace(0,5,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES
%   [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate
%       Statistical Analysis. New York: Wiley, 3rd Ed., 2003.   
%   [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
%       near-exact distribution theory for the most common likelihood ratio
%       test statistics used in Multivariate Analysis. Test, 2011, 20.1:
%       180-203.
%
% See also: cfTest_Sphericity, cfTest_EqualityCovariances,
%           cfTest_CompoundSymmetry, cfTest_Independence,
%           cfTest_EqualityPopulations

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '09-Dec-2018 09:23:48

%% ALGORITHM
% cf = cfTest_EqualityMeans(t,n,p,q,type)

%% CHECK THE INPUT PARAMETERS
narginchk(3,5);

if nargin < 5, type = []; end
if nargin < 4, q = []; end

if isempty(type)
    type = 'modified';
end

if isempty(q)
    q = 2;
end

%% Characteristic function of the null distribution
N = p+q-1;
if n <= N
error('Sample size n is too small')
end

ind   = (1:p)';
alpha = (n-q-ind+1)/2;
beta  = (q-1)/2;

switch lower(type)
    case {'standard', 's'}
        coef = -n/2;
    case  {'modified', 'm'}
        coef = -1;
    otherwise
        coef = -1;
end

cf   = cf_LogRV_Beta(t,alpha,beta,coef);

end