function cf = cfTest_Sphericity(t,n,p,type)
%% cfTest_Sphericity 
%  Characteristic function of the exact null distribution of the
%  multivariate analysis (MVA) test statistic for testing the shape of the
%  unknown covariance matrix Sigma of p-dimensional normal distribution.
%  For more details see Anderson (2003) and/or Marques and Coelho (2011). 
%
%  In particular, here we consider minus log-transformed LRT statistic
%  (Likelihood Ratio Test) or alternatively the minus log-transformed
%  modified LRT statistic. 
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
% SYNTAX
%   cf = cfTest_Sphericity(t,n,p,type)
%
% INPUTS
%  t    - vector or array of real values, where the CF is evaluated.
%  n    - the sample size from the p-dimensional population. It is assumed
%         that n > p.   
%  p    - dimension of the population (p>=2).
%  type - specifies the type of the LRT test statistic: type = 'standard'
%         with W = -(n/2)*log(Lambda), or type = 'modified' with W =
%         -log(Lambda). If type is empty or missing, default value is type
%         = 'modified'.
%
% EXAMPLE 1:
% % CF of the log-transformed LRT statistic for testing sphericity
%   t  = linspace(-2,2,201);
%   n  = 10; 
%   p  = 5;
%   type = 'standard';
%   cf = cfTest_Sphericity(t,n,p,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing sphericity')
%
% EXAMPLE 2:
% % CF of the log-transformed modified LRT statistic for testing sphericity
%   t  = linspace(-10,10,201);
%   n  = 10; 
%   p  = 5;
%   type = 'modified';
%   cf = cfTest_Sphericity(t,n,p,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing sphericity')
%
% EXAMPLE 3:
% % PDF/CDFF/QF of the log-transformed LRT for testing sphericity
%   n  = 10; 
%   p  = 5;
%   type = 'standard';
%   cf   = @(t) cfTest_Sphericity(t,n,p,type);
%   x    = linspace(0,30,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 4:
% % PDF/CDFF/QF of the log-transformed modified LRT for testing sphericity
%   n  = 10; 
%   p  = 5;
%   type = 'modified';
%   cf   = @(t) cfTest_Sphericity(t,n,p,type);
%   x    = linspace(0,6,201);
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
% cf = cfTest_Sphericity(t,n,p,type)

%% CHECK THE INPUT PARAMETERS
narginchk(3,4);

if nargin < 4, type = []; end

if isempty(type)
    type = 'modified';
end

%% Characteristic function of the null distribution
if n <= p 
error('Sample size n is too small')
end

ind   = (2:p)';
alpha = (n-ind)/2;
beta  = (ind-1)/p + (ind-1)/2;

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