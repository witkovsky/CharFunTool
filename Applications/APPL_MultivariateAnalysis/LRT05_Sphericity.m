function [pval,result] = LRT05_Sphericity(W,n,p,options)
%% LRT05_Sphericity computes p-value of the log-transformed LRT statistic,
%  W = -log(Lambda), for testing the null hypothesis on the shape of the
%  unknown covariance matrix Sigma of p-dimensional normal distribution,
%  and/or its null distribution CF/PDF/CDF. In particular, let X ~
%  N_p(mu,Sigma), then, the null hypothesis is given as  
%    H0: Sigma = sigma^2 I_p (sigma^2 unspecified).
%  Here, the LRT test statistic is given by  
%    Lambda = ( det(S) / trace((1/p)*S)^p )^{n/2},
%  where S is MLE of Sigma based on sample size n from X ~ N_p(mu,Sigma).
%
%  Under null hypothesis, distribution of the test statistic Lambda is
%    Lambda ~  prod_{j=2}^{p} (B_j)^{n/2}, 
%  with B_j ~ Beta((n-j)/2,(j-1)/p + (j-1)/2), where j = [2,...,p]. Here
%  we assume that n > p.  
%
%  Hence, the exact characteristic function of the null distribution of
%  minus log-transformed LRT statistic Lambda, say W = -log(Lambda) is
%  given by 
%   cf = @(t) cf_LogRV_Beta(-(n/2)*t, (n-j)/2, (j-1)/p + (j-1)/2),
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
%            coef - set arbitrary multiplicator of the argument t of
%            the characteristic function. If empty, default value is -n/2
%            (standard value for minus log-transform of LRT). Possible
%            alternative is e.g. coef = -1, leading to W = -(2/n)*log(LRT).
%
% EXAMPLE: (LRT for testing hypothesis on sphericity of covariance matrix)
% % Null distribution of the minus log-transformed LRT statistic
% n = 30;           % total sample size
% p = 8;            % dimension of X ~ N_p(mu,Sigma)
% W = [];           % observed value of W = -log(Lambda)
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

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Jan-2018 12:43:15

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

%% 
if n <= p 
error('Sample size n is too small')
end


% CHARACTERISTIC FUNCTION CF
coef = options.coef;
if isempty(coef)
    coef = -n/2;  % Set this option for using with W = -log(LRT)
    % coef = -1;  % Set this options for using normalized -log(LRT^(2/n))
end
ind   = (2:p)';
alpha = (n-ind)/2;
beta  = (ind-1)/p + (ind-1)/2;
cf    = @(t)cf_LogRV_Beta(t,alpha,beta,coef);

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