function [pval,result] = LRT05_Sphericity(W,n,p,options)
%% LRT05_Sphericity computes p-value of the log-transformed LRT statistic,
%  for testing the null hypothesis on the shape of the unknown covariance
%  matrix Sigma of p-dimensional normal distribution, and/or its null
%  distribution CF/PDF/CDF.  
%
%  Let X ~ N_p(mu,Sigma) with p >= 2. Then the null hypothesis is given as  
%    H0: Sigma = sigma^2 I_p (sigma^2 > 0 is unspecified).
%  Here, the LRT test statistic is given by  
%    LRT = Lambda^(n/2) = (det(S) / trace((1/p)*S)^p)^(n/2),
%  and the modified LRT is defined as
%    MLRT = Lambda = det(S) / trace((1/p)*S)^p,
%  where Lambda = det(S) / trace((1/p)*S)^p, where S is MLE of Sigma based
%  on sample size n from X ~ N_p(mu,Sigma). 
%
%  Under null hypothesis, distribution of the test statistic LRT is
%    LRT ~  prod_{j=2}^{p} (B_j)^{n/2}, 
%  with B_j ~ Beta((n-j)/2,(j-1)/p + (j-1)/2), where j = [2,...,p]. Here
%  we assume that n > p. 
%
%  Hence, the exact characteristic function of the null distribution of
%  minus log-transformed LRT statistic Lambda, say W = -log(LRT) is
%  given by 
%   cf = @(t) cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j-1)/p + (j-1)/2),
%  where j = [2,..., p]'.
%  Similarly, the exact characteristic function of the null distribution of
%  minus log-transformed modified LRT statistic, say W = -log(MLRT) is
%   cf = @(t) cf_LogRV_Beta(-t, (n-j)/2, (j-1)/p + (j-1)/2),
%  where j = [2,..., p]'.
%
% SYNTAX:
%  pval = LRT05_Sphericity(W,n,p,options)
%  [pval,result] = LRT05_Sphericity(W,n,p,options)
%  [~,result] = LRT05_Sphericity(W,n,p,options)
%
% INPUTS:
%  W       - observed value of the minus log-transformed LRT statistic
%            W = -log(Lambda). If empty [], the  algorithm evaluates the
%            CF/PDF/CDF and the quantiles of the null distribution of W.  
%  n       - sample size, n > p.
%  p       - dimension of the vector X ~ N_p(mu,Sigma).
%  options - option structure, for more details see cf2DistGP. Moreover,
%            x    - set vector of values where PDF/CDF is evaluated
%            prob - set vector of probabilities for the quantiles.
%            type - specifies the type of the LRT test statistic: type =
%            'standard' with W = -(n/2)*log(Lambda), or type = 'modified'
%            with W = -log(Lambda). If type is empty or missing, default
%            value is type  = 'modified'.
%
% EXAMPLE: (LRT for testing hypothesis on sphericity of covariance matrix)
% % Null distribution of the minus log-transformed LRT statistic
% n = 30;           % total sample size
% p = 8;            % dimension of X ~ N_p(mu,Sigma)
% W = [];           % observed value of W = -log(Lambda)
% clear options;
% % options.coef = -1;
% options.prob = [0.9 0.95 0.99];
% [pval,result] = LRT05_Sphericity(W,n,p,options)
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
% [pval,result] = LRT05_Sphericity(W,n,p,options)

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
if isempty(type)
    type = 'standard';
end

cf = @(t) cfTest_Sphericity(t,n,p,type);

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