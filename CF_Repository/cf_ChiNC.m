function cf = cf_ChiNC(t,df,delta,coef,niid,tol)
%% cf_ChiNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central chi-distributed random variables, with
%  distributions ChiNC(df_i,delta_i), with df_i > 0 degrees of freedom and
%  the noncentrality parameters delta_i >= 0, for all i = 1,...,N.
%
%  cf_ChiNC evaluates the characteristic function cf(t) of Y = coef_i*X_1 +
%  ... + coef_N*X_N, where X_i ~ ChiNC(df_i,delta_i) are inedependent chi
%  distributed RVs, with df_i >0 degrees of freedom and the non-centrality
%  parameters delta_i > 0, for i = 1,...,N.
%
%  The characteristic function of X ~ ChiNC(df,delta) is Poisson mixture
%  (with the intensity parameter lambda = delta^2) of CFs of the central
%  chi-distributed RVs of the form  
%   cf(t) = cf_ChiNC(t,df,delta) = 
%         = exp(-delta^2/2)*sum_{j=1}^Inf (delta^2/2)^j/j! * cf_Chi(df+2*j) 
%  where cf_Chi(df+2*j) are the CFs of central chi-distributed RVs with
%  df+2*j degrees of freedom.  
%
%  Hence, the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is cf_Y(t) = cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters df_i and delta_i.
%
% SYNTAX
%  cf = cf_ChiNC(t,df,delta,coef,niid,tol)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the  degrees of freedom df > 0. If empty, default
%          value is df = 1.  
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0. Notice that the noncentrality
%          parameter delta can be interpreted as a square root of the sum
%          of squared means, delta = sqrt(sum_{i=1}^df mu_i^2).
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
%  tol    - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%   https://en.wikipedia.org/wiki/Rice_distribution
%
% REMARKS:
%  The noncentral chi distribution is a generalization of the chi
%  distribution. If Z_i ~ N(mu_i,1) are k independent, normally distributed
%  random  variables with means mu_i and variances sigma_i^2 = 1, then the
%  statistic X = sqrt(sum_{i=1}^{k} Z_i^2) is distributed according to the
%  noncentral chi distribution with df = k degrees of freedom and the
%  non-centrality parameter delta = sqrt(sum_{i=1}^{k} mu_i^2). Notice that
%  here, the non-centrality parameter delta is stated as a square root of
%  the non-centrality parameter of the chi-squared distribution.
%  The probability density function (pdf) of the noncentral
%  chi-distribution is 
%  pdf(x;df,delta) = exp{-(x^2+delta^2)/2} * x^df/2 * delta^(1-df/2) * ...
%                    besseli(df/2-1,delta*x)
% where besseli is a modified Bessel function of the first kind.
%
% NOTE:
%  A noncentral chi distribution with 2 degrees of freedom is equivalent to
%  a Rice distribution with sigma =1. For more details see Wikipedia 
%  https://en.wikipedia.org/wiki/Rice_distribution   
%
% EXAMPLE 1:
% % CF of the non-central chi-distributed RV with df=3 and delta=1
%   df    = 3;
%   delta = 1;
%   t = linspace(-10,10,201);
%   cf = cf_ChiNC(t,df,delta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the non-central chi-distributed RV with df=3 and delta=1')
%
% EXAMPLE 2:
% % CDF/PDF of the non-central chi-distributed RV with df=3 and delta=1
%   df    = 3;
%   delta = 1;
%   cf = @(t) cf_ChiNC(t,df,delta);
%   clear options;
%   options.xMin = 0;
%   x = linspace(0,5);
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3:
% % CDF/PDF of the non-central chi-distributed RV
% % non-central ChiSquare RVs 
%   df    = [3 4 5];
%   delta = [0 1 2];
%   coef  = [1 1 1];
%   cf = @(t) cf_ChiNC(t,df,delta,coef);
%   clear options;
%   options.xMin = 0;
%   x = linspace(0,12);
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,x,12,options)
%
% SEE ALSO:  cf_Chi, cf_LogRV_ChiNC, cf_LogRV_Rice, cf_LogRV_RayleighNC,
% cf_LogRV_MaxwellBoltzmannNC 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-Oct-2018 11:57:58
% Rev.: 28-Apr-2020 13:47:42

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);

if nargin < 6, tol = []; end
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, delta = []; end
if nargin < 2, df = []; end

if isempty(df), df = 1; end
if isempty(delta), delta = 0; end
if isempty(coef), coef = 1; end
if isempty(tol), tol = 1e-12; end

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,df,delta] = distchck(3,coef(:),df(:),delta(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent nc RVs
szt = size(t);
t   = t(:);
cf  = ones(length(t),1);
for i = 1:length(coef)
    cf = cf .* cf_ncChi(coef(i)*t,df(i),delta(i),tol);
end
cf = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end
end

%% Function funCF
function f = cf_ncChi(t,df,delta,tol)
% cf_ncChi Characteristic function of the distribution of the
% non-central chi distributed RV with df degrees of freedom
% and the non-centrality parameter delta > 0. 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 05-Oct-2018 10:30:02

%%
f = 0;
%delta = delta/2;  % This is alternative definition standard noncentrality
                   % parameter defined as for the chi-square distribution
delta  = delta^2/2;
if delta == 0   % Deal with the central distribution
    f = cf_Chi(t,df);
elseif delta > 0
    % Sum the Poisson series of CFs of central iid Chi distributed RVs,
    % poisspdf(j,delta) .* cf_Chi(t,df+2*j)
    j0 = floor(delta/2);
    p0 = exp(-delta + j0 .* log(delta) - gammaln(j0 + 1));
    f  = f + p0 * cf_Chi(t,df+2*j0);
    p  = p0;
    j  = j0-1;
    while j >= 0 && p > tol
        p = p * (j+1) / delta;
        f = f + p * cf_Chi(t,df+2*j);
        j = j - 1;
    end
    p  = p0;
    j  = j0+1;
    i  = 0;
    while p > tol && i <= 5000
        p = p * delta / j;
        f = f + p * cf_Chi(t,df+2*j);
        j = j + 1;
        i = i + 1;
    end
    if (i == 5000)
        warning(message('NoConvergence'));
    end
else
    error('delta should be nonnegative');
end
end