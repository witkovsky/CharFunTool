function cf = cfTest_Independence(t,n,p,q,type)
%% cfTest_Independence 
%  Characteristic function of the exact null distribution of the
%  multivariate analysis (MVA) test statistic for testing independence q
%  normal populations. For more details see Anderson (2003) and/or Marques
%  and Coelho (2011).
%
%  In particular, here we consider minus log-transformed LRT statistic
%  (Likelihood Ratio Test) or alternatively the minus log-transformed
%  modified LRT statistic. 
%
%  Let X_k ~ N_{p_k}(mu_k,Sigma_k) are p_k dimensional random vectors, for
%  k = 1,...,q. Let us denote X = [X_1,...,X_q] and assume X ~
%  N_p(mu,Sigma). Then, the null hypothesis is given as
%    H0: Sigma = diag(Sigma_1,...,Sigma_q), 
%  i.e. the off-diagonal blocks of Sigma are blocks of zeros. Here, the LRT
%  test statistic is given by 
%    LRT = Lambda^(n/2) = (det(S)/prod(det(S_k)))^(n/2) ,
%  and the modified LRT is defined as
%    MLRT = Lambda = det(S)/prod(det(S_k)),
%  where Lambda = det(S)/prod(det(S_k)) with S is MLE of Sigma, and S_k are
%  MLEs of Sigma_k, for k = 1,...,q, based on n samples from the compound
%  vector X = [X_1,...,X_q]. Here we assume that n > p = sum(p_k).
%
%  Hence, the exact characteristic function of the null distribution of
%  minus log-transformed LRT statistic Lambda, say W = -log(LRT) is
%  given by 
%   cf = @(t) cf_LogRV_Beta(-(n/2)*t, (n-qk-j)/2, qk/2),
%  where qk = [qk_1,...,qk_q] with qk_k = p_{k+1} + ... + p_m, and j =
%  [j_1,...,j_m] with j_k = 1:p_k.
%  Similarly, the exact characteristic function of the null distribution of
%  minus log-transformed modified LRT statistic, say W = -log(MLRT) is
%   cf = @(t) cf_LogRV_Beta(-t, (n-qk-j)/2, qk/2),
%  where qk = [qk_1,...,qk_m] with qk_k = p_{k+1} + ... + p_m, and j =
%  [j_1,...,j_m] with j_k = 1:p_k.
%
% SYNTAX
%   cf = cfTest_Independence(t,n,p,q,type)
%
% INPUTS
%  t    - vector or array of real values, where the CF is evaluated.
%  n    - number of samples from the compound vector X = [X_1,...,X_q]. It
%         is assumed that n > length(X).
%  p    - q-dimensional vector of dimensions of the vectors X_1,...,X_q.
%         If p is scalar value, then the algorithm assumes equal dimension
%         p for all q vectors X_1,...,X_q.
%  q    - number of populations, q > 1. If p is vector with dimension >= 2,
%         then q must be equal to the value q = length(p). If empty or
%         missing, the default value of q is q = length(p) (p is vector)
%         or q = 2 (p is scalar).
%  type - specifies the type of the LRT test statistic: type = 'standard'
%         with W = -(n/2)*log(Lambda), or type = 'modified' with W =
%         -log(Lambda). If type is empty or missing, default value is type
%         = 'modified'.
%
% EXAMPLE 1:
% % CF of the log-transformed LRT statistic for independence
%   t    = linspace(-1,1,201);
%   n    = 20;
%   p    = [4 5 6] ;
%   q    = [];
%   type = 'standard';
%   cf    = cfTest_Independence(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of LRT for testing independence of q normal populations')
%
% EXAMPLE 2:
% % CF of the log-transformed LRT statistic for independence
%   t  = linspace(-1,1,201);
%   n  = 20;
%   p  = 5;
%   q  = 3;
%   type = 'standard';
%   cf = cfTest_Independence(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of LRT for testing independence of q normal populations')
%
% EXAMPLE 3:
% % CF of the log-transformed modified LRT statistic for independence
%   t  = linspace(-5,5,201);
%   n  = 20;
%   p  = 5;
%   q  = 3;
%   type = 'modified';
%   cf = cfTest_Independence(t,n,p,q,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of LRT for testing independence of q normal populations')
%
% EXAMPLE 4:
% % PDF/CDFF/QF of the log-transformed LRT for independence
%   n    = 20;
%   p    = 5;
%   q    = 3;
%   type = 'standard';
%   cf   = @(t) cfTest_Independence(t,n,p,q,type);
%   x    = linspace(30,120,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 5:
% % PDF/CDFF/QF of the log-transformed modified LRT for independence
%   n    = 20;
%   p    = 5;
%   q    = 3;
%   cf   = @(t) cfTest_Independence(t,n,p,q);
%   x    = linspace(2,12,201);
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
%           cfTest_CompoundSymmetry, cfTest_EqualityMeans,
%           cfTest_EqualityPopulations

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '09-Dec-2018 09:23:48

%% ALGORITHM
% cf = cfTest_Independence(t,n,p,q,type)

%% CHECK THE INPUT PARAMETERS
narginchk(3,5);

if nargin < 5, type = []; end
if nargin < 4, q = []; end

if isempty(type)
    type = 'modified';
end

if isempty(q)
    if length(p)>1
        q   = length(p);
    else
        q = 2;
        p = [p p];
    end 
else
    if length(p)>1 && length(p)~=q
        warning('Dimensions mismatch. Dimension of p must be equal to q')
        q = length(p);
    elseif length(p)==1
        p = p*ones(1,q);
    end
end

%% Characteristic function of the null distribution
csp = cumsum(p);
if n <= csp(q)
error('Sample size n is too small')
end

alpha = zeros(csp(q-1),1);
beta  = zeros(csp(q-1),1);
ind   = 0;
for k = 1:(q-1)
    qk = csp(q)-csp(k);
    for j = 1:p(k)
        ind = ind +1;
        alpha(ind) = (n - qk - j) / 2;
        beta(ind)  = qk / 2;
    end
end

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