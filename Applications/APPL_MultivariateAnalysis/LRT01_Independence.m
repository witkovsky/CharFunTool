function [pval,result] = LRT01_Independence(W,n,p,options)
%% LRT01_Independence computes p-value of the log-transformed LRT statistic,
%  W = -log(Lambda), for testing the null hypothesis of independence (under
%  normality assumptions) of m groups of variables (m>1), and/or its null
%  distribution CF/PDF/CDF. In particular, let X_k ~ N_{p_k}(mu_k,Sigma_k)
%  are p_k dimensional random vectors, for k = 1,...,m. Let us denote X =
%  [X_1,...,X_m] and assume X ~ N_p(mu,Sigma). Then, the null hypothesis is
%  given as
%    H0: Sigma = diag(Sigma_1,...,Sigma_m), 
%  i.e. the off-diagonal blocks of Sigma are blocks of zeros. Here, the LRT
%  test statistic is given by 
%    Lambda = det(S) /  prod(det(S_k)),
%  where S is MLE of Sigma, and S_k are MLEs of Sigma_k, for k = 1,...,m, 
%  based on n samples from the compound vector X = [X_1,...,X_m].
%
%  Under null hypothesis, distribution of the test statistic Lambda is
%    Lambda ~  prod_{k=1}^{m-1} prod_{j=1}^{p_k} (B_{jk})^{n/2}, 
%  with B_{jk} ~ Beta((n-q_k-j)/2,q_k/2), where q_k = p_{k+1} + ... + p_m.
%  Here we assume that n > q_k, for all k = 1,...,m-1.
%
%  Hence, the exact characteristic function of the null distribution of
%  minus log-transformed LRT statistic Lambda, say W = -log(Lambda) is
%  given by 
%   cf = @(t) cf_LogRV_Beta(-(n/2)*t, (n-q-j)/2, q/2),
%  where q = [q_1,...,q_m] with q_k = p_{k+1} + ... + p_m, and j =
%  [j_1,...,j_m] with j_k = 1:p_k.
%
% SYNTAX:
%  pval = LRT01_Independence(W,n,p,options)
%  [pval,result] = LRT01_Independence(W,n,p,options)
%  [~,result] = LRT01_Independence([],n,p,options)
%
% INPUTS:
%  W       - observed value of the minus log-transformed LRT statistic
%            W = -log(Lambda). If empty [], the  algorithm evaluates the
%            CF/PDF/CDF and the quantiles of the null distribution of W.  
%  n       - sample size (n > q_k).
%  p       - vector p = [p_1,...,p_k] of dimensions of X_k, k = 1,...m.
%  options - option structure, for more details see cf2DistGP. Moreover,
%            x    - set vector of values where PDF/CDF is evaluated
%            prob - set vector of probabilities for the quantiles.
%            coef - set arbitrary multiplicator of the argument t of
%            the characteristic function. If empty, default value is -n/2
%            (standard value for minus log-transform of LRT). Possible
%            alternative is e.g. coef = -1, leading to W = -(2/n)*log(LRT).
%
% EXAMPLE: (LRT for testing hypothesis about independence)
% % Null distribution of the minus log-transformed LRT statistic
% n = 30;            % sample size
% p = [3 4 5 6 7];   % dimensions of X_k,k = 1,...,m where m = 5
% W = [];            % observed value of W = -log(Lambda)
% clear options;
% % options.coef = -1;
% options.prob = [0.9 0.95 0.99];
% [pval,result] = LRT01_Independence(W,n,p,options)
%
% REFERENCES:
%   [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate
%       Statistical Analysis. New York: Wiley, 3rd Ed., 2003.   
%   [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
%       near-exact distribution theory for the most common likelihood ratio
%       test statistics used in Multivariate Analysis. Test, 2011, 20.1:
%       180-203.
%   [3] WITKOVSKÝ, Viktor. Exact distribution of selected multivariate test
%       criteria by numerical inversion of their characteristic functions.
%       arXiv preprint arXiv:1801.02248, 2018.   

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Jan-2018 12:43:15

%% ALGORITHM
% [pval,result] = LRT01_Independence(W,n,p,options)

%% CHECK THE INPUT PARAMETERS
narginchk(3,4);

if nargin < 4, options = []; end

if ~isfield(options, 'x')
    options.x = [];
end

if ~isfield(options, 'prob')
    options.prob = [];
end

if ~isfield(options, 'coef')
    options.coef = [];
end

if ~isfield(options, 'xMin')
    options.xMin = 0;
end

%% 
m = length(p);
q = cumsum(p);
N = q(m-1);

if n <= N 
error('Sample size n is too small')
end

alpha = zeros(N,1);
beta  = zeros(N,1);
ind   = 0;
for k = 1:(m-1)
    for j = 1:p(k)
        ind = ind +1;
        alpha(ind) = (n - q(k) - j) / 2;
        beta(ind)  = q(k) / 2;       
    end
end


% CHARACTERISTIC FUNCTION CF
coef = options.coef;
if isempty(coef)
    coef = -n/2;  % Set this option for using with W = -log(LRT)
    % coef = -1;  % Set this options for using normalized -log(LRT^(2/n))
end
cf = @(t)cf_LogRV_Beta(t,alpha,beta,coef);

% Evaluate the p-value and PDF/CDF/QF of the log-transformed LRT statistic
if ~isempty(W)
    % P-VALUE
    options.isPlot = false;
    result = cf2DistGP(cf,W,[],options);
    pval   = 1-result.cdf;
else
    % DISTRIBUTION of Lambda PDF/CDF
    pval   = [];
    result = cf2DistGP(cf,options.x,options.prob,options);
end

% Save the parameters of the used beta distributions
result.alpha = alpha;
result.beta  = beta;

end