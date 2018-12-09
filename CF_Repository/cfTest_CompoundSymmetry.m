function cf = cfTest_CompoundSymmetry(t,n,p,type)
%% cfTest_CompoundSymmetry 
%  Characteristic function of the exact null distribution of the
%  multivariate analysis (MVA) test statistic for testing the shape of the
%  unknown covariance matrix Sigma of p-dimensional normal distribution.
%  Here it is assumed that under null hypothesis the covariance matric has
%  the compound symmetry structure, i.e. Sigma = sigma^2*((1-rho)*I+rho*J),
%  where I is the p-dimensional identity matrix and J is pxp matrix of
%  ones. For more details see Anderson (2003), Marques and Coelho (2011),
%  Nagar and Castaneda(2006), Marques and Coelho (2013 and 2018).
%
%  In particular, here we consider minus log-transformed LRT statistic
%  (Likelihood Ratio Test) or alternatively the minus log-transformed
%  modified LRT statistic. 
%
%  Let X ~ N_p(mu,Sigma) with p >= 2. Then the null hypothesis is given as  
%    H0: Sigma = sigma^2*((1-rho)*I_p + rho*J_p) (sigma^2>0 unspecified).
%  Here, the LRT test statistic is given by  
%    LRT = Lambda^(n/2) = ((p-1)^(p-1) * p^p * det(S) / ...
%                         (trace(J*S) * trace((p*I-J)*S)^(p-1))^(n/2),
%  and the modified LRT is defined as
%    MLRT = Lambda = (p-1)^(p-1) * p^p * det(S) / ...
%                         (trace(J*S) * trace((p*I-J)*S)^(p-1),
%  where S is MLE of Sigma based on sample size n from X ~ N_p(mu,Sigma). 
%
%  Under null hypothesis, distribution of the test statistic LRT is
%    LRT ~  prod_{j=2}^{p} (B_j)^{n/2}, 
%  with B_j ~ Beta((n-j+1)/2,(j-2)/(p-1) + (j-1)/2), where j = [2,...,p].
%  Here we assume that n > p. 
%
%  Hence, the exact characteristic function of the null distribution of
%  minus log-transformed LRT statistic Lambda, say W = -log(LRT) is
%  given by 
%   cf = @(t) cf_LogRV_Beta(-(n/2)*t, (n-j+1)/2, (j-2)/(p-1) + (j-1)/2),
%  where j = [2,..., p]'.
%  Similarly, the exact characteristic function of the null distribution of
%  minus log-transformed modified LRT statistic, say W = -log(MLRT) is
%   cf = @(t) cf_LogRV_Beta(-t, (n-j+1)/2, (j-2)/(p-1) + (j-1)/2),
%  where j = [2,..., p]'.
% 
% SYNTAX
%   cf = cfTest_CompoundSymmetry(t,n,p,type)
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
% % CF of the log-transformed LRT statistic for testing compound symmetry
%   t  = linspace(-2,2,201);
%   n  = 10; 
%   p  = 5;
%   type = 'standard';
%   cf = cfTest_CompoundSymmetry(t,n,p,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing compound symmetry')
%
% EXAMPLE 2:
% % CF of the log-transformed modified LRT statistic for compound symmetry
%   t  = linspace(-10,10,201);
%   n  = 10; 
%   p  = 5;
%   type = 'modified';
%   cf = cfTest_CompoundSymmetry(t,n,p,type);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of test statistic for testing compound symmetry')
%
% EXAMPLE 3:
% % PDF/CDFF/QF of the log-transformed LRT for testing compound symmetry
%   n  = 10; 
%   p  = 5;
%   type = 'standard';
%   cf   = @(t) cfTest_CompoundSymmetry(t,n,p,type);
%   x    = linspace(0,25,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 4:
% % PDF/CDFF/QF of the log-transformed modified LRT for compound symmetry
%   n  = 10; 
%   p  = 5;
%   type = 'modified';
%   cf   = @(t) cfTest_CompoundSymmetry(t,n,p,type);
%   x    = linspace(0,5,201);
%   prob = [0.9 0.95 0.99];
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES
%  [1] ANDERSON, Theodore Wilbur. An Introduction to Multivariate
%      Statistical Analysis. New York: Wiley, 3rd Ed., 2003.   
%  [2] MARQUES, Filipe J.; COELHO, Carlos A.; ARNOLD, Barry C. A general
%      near-exact distribution theory for the most common likelihood ratio
%      test statistics used in Multivariate Analysis. Test, 2011, 20.1:
%      180-203. 
%  [3] NAGAR, D.K., CASTANEDA, E.M. Distribution and percentage
%      points of LRC for testing multisample compound symmetry in the
%      bivariate and the trivariate cases. Metron 64, no. 2 (2006): 217-238.       
%  [4] MARQUES, Filipe J.; COELHO, Carlos A. The multisample block-diagonal
%      equicorrelation and equivariance test, AIP Conference Proceedings
%      1558 (2013) pp. 793-796.  
%  [5] MARQUES, Filipe J.; COELHO, Carlos A. Testing simultaneously
%      different covariance block diagonal structures - the multisample
%      case. 2018 Working paper. 
%
% See also: cfTest_CompoundSymmetry, cfTest_EqualityCovariances,
%           cfTest_CompoundSymmetry, cfTest_Independence,
%           cfTest_EqualityPopulations

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '09-Dec-2018 09:23:48

%% ALGORITHM
% cf = cfTest_CompoundSymmetry(t,n,p,type)

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
alpha = (n-ind+1)/2;
beta  = (ind-2)/(p-1) + (ind-1)/2;

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