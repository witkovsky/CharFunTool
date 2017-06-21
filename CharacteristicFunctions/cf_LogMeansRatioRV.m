function cf = cf_LogMeansRatioRV(t,n,alpha,coef,niid)
%%cf_LogMeansRatioRV Characteristic function of a linear combination (resp.
%  convolution) of independent log-transformed random variables (RVs)
%  W_i = log(R_i), where each R_i = G_i/A_i is a ratio of G_i (the
%  geometric mean) and A_i (the arithmetic mean) of a random sample
%  X_i1,...,X_in from a gamma distribution. 
%
%  Here, we define G_i = (X_i1 * ... * X_in)^(1/n_i) and A_i = (X_i1 +...+
%  X_in)/n_i, with X_ij ~ Gamma(alpha_i,beta_i) for all j = 1,...,n_i,
%  where alpha_i is the shape parameter and beta_i is the rate parameter of
%  the gamma distribution in the i-th population. Note that the R_i random
%  variables are scale invariant, so the distribution does not depend on
%  the rate (scale) parameters beta_i.
%  
%  That is, cf_LogMeansRatioRV evaluates the characteristic function of a
%  random variable Y  = coef_1*W_1 + ... + coef_N*W_N, such that
%  cf_Y(t) =  cf_W_1(coef_1*t) * ... * cf_W_N(coef_N*t), where cf_W_i(t)
%  is CF of W_i = log(R_i), and  R_i ~ MeansRatio(alpha_i,n_i), where
%  MeansRatio(alpha_i,n_i) denotes the distribution of R_i.
%
%  The MeansRatio distribution of R_i ~ MeansRatio(alpha_i,n_i), with
%  R_i in (0,1), is defined by R_i ~ (Prod_{j=1}^{n_i-1} B_ij)^(1/n_i),
%  where B_ij ~ Beta{alpha_i, j/n_i)} for j = 1,...,n_i-1, i.e. B_ij follow
%  independent Beta distributions for all i = 1,...,N and j = 1,...,n_i-1.
%  Alternatively, log(R_i) ~ (log(B_i1) + ... + log(B_{i,n_i-1}))/n_i.
%
% SYNTAX
%  cf = cf_LogMeansRatioRV(t,n,alpha,coef,niid)
%
% INPUTS
%  t     - vector or array of real values, where the CF is evaluated.
%  n     - vector of sample size parameters n = (n_1,...,n_N). If empty,
%          default value is n = (1,...,1). 
%  alpha - vector of the 'shape' parameters alpha = (alpha_1,...,alpha_N).
%          If empty, default value is alpha = (1,...,1). 
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically
%          distributed random variable. If empty, default value is n = 1.  
%
% EXAMPLE 1:
% % CF of log MeansRatio RV with alpha = 7/2 and n = 5 
%   n     = 5;
%   alpha = 7/2;
%   t     = linspace(-100,100,201);
%   cf    =  cf_LogMeansRatioRV(t,n,alpha);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of log MeansRatio RV with alpha = 7/2 and n = 5')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of minus log MeansRatio RVs 
%   n     = [5 7 10];
%   alpha = [7 10 15]/2;
%   coef  = -[5 7 10]/22;
%   t     = linspace(-100,100,201);
%   cf    =  cf_LogMeansRatioRV(t,n,alpha,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a weighted linear combination of minus log MeansRatio RVs')
%
% EXAMPLE 3:
% % PDF/CDF of minus log MeansRatio RV, alpha = 7/2 and n = 5, from its CF
%   n     = 5;
%   alpha = 7/2;
%   coef  = -1;
%   cf    = @(t) cf_LogMeansRatioRV(t,n,alpha,coef);
%   x = linspace(0,0.6)';
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 4:
% % Compare the the exact distribution with the Bartlett's approximation
%   df    = 3;
%   n     = 25;
%   DF    = n*df; 
%   alpha = df/2;
%   const = (1 + 1/(3*(n-1))*(n/df - 1/DF))/DF;
%   coef  = -1/const;
%   cf    = @(t) cf_LogMeansRatioRV(t,n,alpha,coef);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%   x = result.x;
%   figure;plot(x,result.cdf,x,chi2cdf(x,n-1));grid
%   title('Exact CDF vs. the Bartlett approximation')
%   xlabel('corrected test statistic')
%   ylabel('CDF')
%   disp(prob)
%   disp(result.qf)
%   disp(chi2inv(prob,n-1))
%
% EXAMPLE 5:
% % Exact Critical Values for Bartlett's Test for Homogeneity of Variances
% % See and compare the selected results in Glaser (1976b, Table 1)
%   n     = 3;
%   df    = [4 5 6 7 8 9 10 11 14 19 24 29 49 99];
%   alpha = df/2;
%   coef  = -1;
%   prob  = [0.9 0.95 0.99 0.999 0.9999];
%   clear options
%   options.N = 2^12;
%   options.SixSigmaRule = 15;
%   options.xMin = 0;
%   options.isPlot = false;
%   critW = zeros(length(df),length(prob));
%   for i = 1:length(df)
%      cf    = @(t) cf_LogMeansRatioRV(t,n,alpha(i),coef);
%      result = cf2DistGP(cf,[],prob,options);
%      critW(i,:) = result.qf;
%   end;
%   critR = exp(-critW);
%   disp(['n = ',num2str(n)])
%   disp('alpha = ')
%   disp(1-prob)
%   disp(critR)
%
%  REMARKS
%  The MeansRatio distribution is the distribution of the Bartlett's test
%  statistitic for testing homogeneity of variances of n populations, based
%  on equal sample sizes of size m, i.e. with df = m-1 for all j =
%  1,...,n. Let DF = n*df. The Bartlett (1937) test statistic is defined by
%  R = Prod((df*S^2_j/sigma^2)^(df/DF)) / (df/DF)*Sum(df*S^2_j/sigma^2)
%  = Prod((S^2_j)^(1/n)) / Sum(S^2_j)/n, where S^2_j = (1/(m-1))*Sum(X_jk -
%  mean(X_jk))^2. Notice that df*S^2_j/sigma^2 ~ Chi^2_df = Gamma(df/2,1/2)
%  for all j = 1,...,n. The exact critical values can be calculated from
%  the distribution of W = -log(R) by inverting the characteristic function
%  cf_LogMeansRatioRV.  
%  
%  The corrected version of the test statistic W with better convergence to
%  the asymptotic chi-square distribution with n-1 degrees of freedom is
%  Chi2 = (log(S^2_pool) - Sum(log(S^2_j))/n) / const, where the correction
%  constant is const = (1 + 1/(3*(n-1))*(sum(1/df) - 1/DF))/DF.
%
%  WIKIPEDIA
%  https://en.wikipedia.org/wiki/Bartlett%27s_test
%
%  REFERENCES
%  Glaser, R. E. (1976a). The ratio of the geometric mean to the arithmetic
%  mean for a random sample from a gamma distribution. Journal of the
%  American Statistical Association, 71(354), 480-487.   
%
%  Glaser, R. E. (1976b). Exact critical values for Bartlett's test for
%  homogeneity of variances. Journal of the American Statistical
%  Association, 71(354), 488-490.  
%
%  Chao, M. T., & Glaser, R. E. (1978). The exact distribution of
%  Bartlett's test statistic for homogeneity of variances with unequal
%  sample sizes. Journal of the American Statistical Association, 73(362),
%  422-426.   

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Jun-2017 17:18:39

%% ALGORITHM
% cf = cf_LogMeansRatioRV(t,n,alpha,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, n = []; end
if nargin < 2, alpha = []; end

%%
if isempty(alpha), alpha = 1; end
if isempty(n), n = 1; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Check size of the parameters
[errorcode,alpha,n,coef] = distchck(3,alpha(:)',n(:)',coef(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:length(coef)
    if n(i)>1
        beta  = (1:(n(i)-1))/n(i);
        cf = cf .* cf_LogBetaRV(coef(i)*t,alpha(i),beta,1/n(i));
    end
end
cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end