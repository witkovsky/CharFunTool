function cf = cfTest_EqualityCovariances(t,n,p,q,type)
%% cfTest_EqualityCovariances 
%  Characteristic function of the exact null distribution of the
%  multivariate analysis (MVA) test statistic for testing equality of
%  covariance matrices of q normal p-dimensional populations (q>1). For
%  more details see Anderson (2003) and/or Marques and Coelho (2011). 
%
%  In particular, here we consider minus log-transformed LRT statistic
%  (Likelihood Ratio Test) or alternatively the minus log-transformed
%  modified LRT statistic. 
%
%  Let X_k ~ N_p(mu_k,Sigma_k) are p-dimensional random vectors, for k =
%  1,...,q. We want to test the hypothesis that the covariance matrix Sigma
%  is common for all X_k, k = 1,...,q. Then, the null hypothesis is given
%  as
%    H0: Sigma_1 = ... = Sigma_q, 
%  i.e. the covariance matrices are equal in all q populations. Here, the
%  LRT test statistic is given by
%    LRT = Lambda^(n/2) = (q^{p*q} * prod(det(S_k)) / (det(S))^q)^(n/2),
%  and the modified LRT is defined as
%    MLRT = Lambda = q^{p*q} * prod(det(S_k)) / (det(S))^q,
%  where Lambda = q^{p*q} * prod(det(S_k)) / (det(S))^q where S_k are MLEs
%  of Sigma_k, for k = 1,...,q, and S = S_1 + ... + S_q, based on n samples
%  from each of the the q p-dimensional populations.   
%
%  Under null hypothesis, distribution of the LRT statistic is
%    LRT ~  prod_{k=1}^q prod_{j=1}^{p} (B_{jk})^{n/2}, 
%  with B_{jk} ~ Beta((n-j)/2,(j*(q-1)+2*k-1-q)/2), and we set B_{11} = 1
%  for j=k=1. Here we assume that n > p.
%
%  Hence, the exact characteristic function of the null distribution of
%  minus log-transformed LRT statistic Lambda, say W = -log(LRT) is
%  given by 
%   cf = @(t) cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q)),
%  where k = [1*o,...,q*o] with p-dimensional vector of ones o = [1,...,1]
%  and j = [j_1,...,j_q] with j_k = 1:p. 
%  Similarly, the exact characteristic function of the null distribution of
%  minus log-transformed modified LRT statistic, say W = -log(MLRT) is
%   cf = @(t) cf_LogRV_Beta(-t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q)),
%  where k = [1*o,...,q*o] with p-dimensional vector of ones o = [1,...,1]
%  and j = [j_1,...,j_q] with j_k = 1:p. 
% 
% SYNTAX
%   cf = cfTest_EqualityCovariances(t,n,p,q,type)
%
% INPUTS
%  t    - vector or array of real values, where the CF is evaluated.
%  n    - number of samples from each of the q p-dimensional populations.
%         It is assumed that n > p. 
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
% % CF of the log-transformed LRT statistic for equal covariances
%   t  = linspace(-1,1,201);
%   n  = 10;
%   p  = 5;
%   q  = 3;
%   type = 'standard';
%   cf = cfTest_EqualityCovariances(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing equal covariances')
%
% EXAMPLE 2:
% % CF of the log-transformed modified LRT statistic for equal covariances
%   t  = linspace(-5,5,201);
%   n  = 10;
%   p  = 5;
%   q  = 3;
%   type = 'modified';
%   cf = cfTest_EqualityCovariances(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing equal covariances')
%
% EXAMPLE 3:
% % PDF/CDFF/QF of the log-transformed LRT for equal covariances
%   n    = 10;
%   p    = 5;
%   q    = 3;
%   type = 'standard';
%   cf   = @(t) cfTest_EqualityCovariances(t,n,p,q,type);
%   x    = linspace(0,50,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 4:
% % PDF/CDFF/QF of the log-transformed modified LRT for equal covariances
%   n    = 10;
%   p    = 5;
%   q    = 3;
%   type = 'modified';
%   cf   = @(t) cfTest_EqualityCovariances(t,n,p,q,type);
%   x    = linspace(0,10,201);
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
% See also: cfTest_Independence, cfTest_Sphericity, 
%           cfTest_CompoundSymmetry, cfTest_EqualityMeans,
%           cfTest_EqualityPopulations

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '09-Dec-2018 09:23:48

%% ALGORITHM
% cf = cfTest_EqualityCovariances(t,n,p,q,type)

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
if n <= p 
error('Sample size n is too small')
end

alpha = zeros(p*q,1);
beta  = zeros(p*q,1);
ind   = 0;
for k = 1:q
    for j = 1:p
        ind = ind +1;
        alpha(ind) = (n - j) / 2;
        beta(ind)  = (j*(q-1) + 2*k - 1 - q) / (2*q);       
    end
end

% For k=j=1 the coefficient beta=0, hence log(Beta(alpha,beta))=0. 
alpha = alpha(2:end);
beta  = beta(2:end);

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