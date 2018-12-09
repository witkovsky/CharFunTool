function [pval,result] = LRT04_EqualityPopulations(W,n,p,q,options)
%% LRT04_EqualityPopulations computes p-value of the log-transformed LRT 
%  statistic W = -log(Lambda), for testing the null hypothesis of equality
%  of q (q>1) p-dimensional normal populations, and/or its null
%  distribution CF/PDF/CDF. 
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
% SYNTAX:
%  pval = LRT04_EqualityPopulations(W,n,p,q,options)
%  [pval,result] = LRT04_EqualityPopulations(W,n,p,q,options)
%  [~,result] = LRT04_EqualityPopulations(W,n,p,q,options)
%
% INPUTS:
%  W       - observed value of the minus log-transformed LRT statistic
%            W = -log(Lambda). If empty [], the  algorithm evaluates the
%            CF/PDF/CDF and the quantiles of the null distribution of W.  
%  n       - sample size, n > min(p+q-1).
%  p       - common dimension of the vectors X_k, k = 1,...q.
%  q       - number of normal populations, q > 1.
%  options - option structure, for more details see cf2DistGP. Moreover,
%            x    - set vector of values where PDF/CDF is evaluated
%            prob - set vector of probabilities for the quantiles.
%            type - specifies the type of the LRT test statistic: type =
%            'standard' with W = -(n/2)*log(Lambda), or type = 'modified'
%            with W = -log(Lambda). If type is empty or missing, default
%            value is type  = 'modified'.
%
% EXAMPLE: (LRT for testing hypothesis on equality of populations)
% % Null distribution of the minus log-transformed LRT statistic
% n = 30;           % total sample size
% p = 8;            % dimension of X_k, k = 1,...,q where q = 5
% q = 5;            % number of populations 
% W = [];           % observed value of W = -log(Lambda)
% clear options;
% options.type = 'standard';
% % options.type = 'modified';
% options.prob = [0.9 0.95 0.99];
% [pval,result] = LRT04_EqualityPopulations(W,n,p,q,options)
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

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 09-Dec-2018 15:34:58

%% ALGORITHM
% [pval,result] = LRT04_EqualityPopulations(W,n,p,q,options)

%% CHECK THE INPUT PARAMETERS
narginchk(4,5);

if nargin < 5, options = []; end

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

if ~isfield(options, 'type')
    options.type = [];
end

type = options.type;

%% 
if isempty(type)
    type = 'standard';
end

cf = @(t) cfTest_EqualityPopulations(t,n,p,q,type);

% Evaluate p-value or PDF/CDF/QF
if ~isempty(W)
    % P-VALUE
    options.isPlot = false;
    result = cf2DistGP(cf,W,[],options);
    pval   = 1-result.cdf;
else
    % PDF/CDF/QF otf the test statistic
    pval   = [];
    result = cf2DistGP(cf,options.x,options.prob,options);
end

end