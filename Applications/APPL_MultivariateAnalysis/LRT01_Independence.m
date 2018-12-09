function [pval,result] = LRT01_Independence(W,n,p,options)
%% LRT01_Independence computes p-value of the log-transformed LRT statistic,
%  W = -log(Lambda), for testing the null hypothesis of independence (under
%  normality assumptions) of m groups of variables (m>1), and/or its null
%  distribution CF/PDF/CDF. 
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
%            type - specifies the type of the LRT test statistic: type =
%            'standard' with W = -(n/2)*log(Lambda), or type = 'modified'
%            with W = -log(Lambda). If type is empty or missing, default
%            value is type  = 'modified'.
%
% EXAMPLE: (LRT for testing hypothesis about independence)
% % Null distribution of the minus log-transformed LRT statistic
% n = 30;            % sample size
% p = [3 4 5 6 7];   % dimensions of X_k,k = 1,...,m where m = 5
% W = [];            % observed value of W = -log(Lambda)
% clear options;
% options.type = 'standard';
% % options.type = 'modified';
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

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 09-Dec-2018 15:34:58

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

if ~isfield(options, 'type')
    options.type = [];
end

type = options.type;

%% 
q = length(p);

if isempty(type)
    type = 'standard';
end

cf = @(t) cfTest_Independence(t,n,p,q,type);

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