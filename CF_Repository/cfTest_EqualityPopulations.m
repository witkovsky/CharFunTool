function cf = cfTest_EqualityPopulations(t,n,p,q,type)
%% cfTest_EqualityPopulations 
%  Characteristic function of the exact null distribution of the
%  multivariate analysis (MVA) test statistic for testing equality of q
%  normal p-dimensional populations. In fact, we simultaneously test
%  equality of means and equality of covariance matrices of q normal
%  p-dimensional populations. For more details see Anderson (2003) and/or
%  Marques and Coelho (2011).
%
%  In particular, here we consider minus log-transformed LRT statistic
%  (Likelihood Ratio Test) or alternatively the minus log-transformed
%  modified LRT statistic. 
%
%  Let X_k ~ N_p(mu_k,Sigma_k), for k = 1,...,q. We want to test the
%  hypothesis that the q normal populations are equally distributed. That
%  is, we want to test that the mean vectors mu_k are equal for all k =
%  1,...,q, as well as the covariance matrices Sigma_k are equal for all k
%  = 1,...,q. Then, the null hypothesis is given as
%    H0: mu_1 = ... = mu_q & Sigma_1 = ... = Sigma_k.
%  Here, the null hypothesis H0 and the LRT statistic can be decomposed:
%    LRT = LRT_Means * LRT_Covariances
%  where (first) LRT_Covariances represents the LRT for testing equality of
%  covariance matrices of given q normal populations, and (second)
%  LRT_Means represents (conditionally) the LRT for testing equality of
%  means of given q normal populations.    
%
%  Under null hypothesis, distributions of LRT_Covariances and LRT_Means
%  are independent, and the distribution of the compound LRT statistic is
%    LRT =  Lambda^(n/2) = LRT_Means * LRT_Covariances, 
%        = (prod_{k=1}^q prod_{j=1}^{p} (B_{jk})^(n/2)) 
%           * (prod_{j=1}^{p} (B_j)^(q*n/2))
%  and the modified LRT is defined as
%    MLRT = Lambda = MLRT_Means * MLRT_Covariances
%         = (prod_{k=1}^q prod_{j=1}^{p} B_{jk}) 
%           * (prod_{j=1}^{p} (B_j)^q)
%  where Lambda = (prod_{k=1}^q prod_{j=1}^p B_{jk})*(prod_{j=1}^{p} B_j^q)
%  where  B_{jk} and B_j are mutually independent beta distributed  
%  random variables. Here we assume that n is equal sample size for each
%  sample, k = 1,...,q, n > p.  
%
%  Hence, the exact characteristic function of the null distribution of
%  minus log-transformed LRT statistic Lambda, say W = -log(LRT) is
%  given by 
%    cf = @(t) cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q)) 
%           .* cf_LogRV_Beta(-(n*q/2)*t, ((n-1)*q-i+1)/2, (q-1)/2),
%  where i = [1, 2, ..., p]', k = [1*o,...,q*o] with p-dimensional vector
%  of ones o = [1,...,1]  and j = [j_1,...,j_q] with j_k = 1:p. 
%  Similarly, the exact characteristic function of the null distribution of
%  minus log-transformed modified LRT statistic, say W = -log(MLRT) is
%    cf = @(t) cf_LogRV_Beta(-t, (n-j)/2, (j*(q-1)+2*k-1-q)/(2*q)) 
%           .* cf_LogRV_Beta(-q*t, ((n-1)*q-i+1)/2, (q-1)/2),
%  where i = [1, 2, ..., p]', k = [1*o,...,q*o] with p-dimensional vector
%  of ones o = [1,...,1]  and j = [j_1,...,j_q] with j_k = 1:p.
% 
% SYNTAX
%   cf = cfTest_EqualityPopulations(t,n,p,q,type)
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
% % CF of the log-transformed LRT statistic for equal populations
%   t  = linspace(-1,1,201);
%   n  = 10;
%   p  = 5;
%   q  = 3;
%   type = 'standard';
%   cf = cfTest_EqualityPopulations(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing equal populations')
%
% EXAMPLE 2:
% % CF of the log-transformed modified LRT statistic for equal populations
%   t  = linspace(-5,5,201);
%   n  = 10;
%   p  = 5;
%   q  = 3;
%   type = 'modified';
%   cf = cfTest_EqualityPopulations(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing equal populations')
%
% EXAMPLE 3:
% % PDF/CDFF/QF of the log-transformed LRT for equal populations
%   n    = 10;
%   p    = 5;
%   q    = 3;
%   type = 'standard';
%   cf   = @(t) cfTest_EqualityPopulations(t,n,p,q,type);
%   x    = linspace(0,60,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 4:
% % PDF/CDFF/QF of the log-transformed modified LRT for equal populations
%   n    = 10;
%   p    = 5;
%   q    = 3;
%   type = 'modified';
%   cf   = @(t) cfTest_EqualityPopulations(t,n,p,q,type);
%   x    = linspace(0,12,201);
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
%           cfTest_EqualityCovariances

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '09-Dec-2018 09:23:48

%% ALGORITHM
% cf = cfTest_EqualityPopulations(t,n,p,q,type)

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

ind   = (1:p)';
alpha = zeros(p*(q+1),1);
beta  = zeros(p*(q+1),1);
coef  = ones(p*(q+1),1); 

alpha(1:p) = ((n-1)*q-ind+1)/2;
beta(1:p)  = (q-1)/2;
coef(1:p)  = q*coef(1:p); 

ind    = p;
for k  = 1:q
    for j = 1:p
        ind = ind +1;
        alpha(ind) = (n - j) / 2;
        beta(ind)  = (j*(q-1) + 2*k - 1 - q) / (2*q);
    end
end

switch lower(type)
    case {'standard', 's'}
        coef = (-n/2)*coef;
    case  {'modified', 'm'}
        coef = -coef;
    otherwise
        coef = -coef;
end

cf   = cf_LogRV_Beta(t,alpha,beta,coef);

end