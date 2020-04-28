function cf = cf_BetaNC(t,alpha,beta,delta,coef,niid,tol,type)
%% cf_BetaNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central BETA random variables (Type I and Type II), with
%  distributions Beta(alpha_i,beta_i,delta_i), defined on the interval
%  (0,1), and specified by the parameters alpha_i, beta_i, and the
%  noncentrality parameters delta_i.
%
%  The noncentral beta distribution has two types. The Type I is the
%  distribution  of the random variable B1 = X1/(X1+X2), X1 ~
%  Gamma(alpha,gamma,delta) and X2 ~ Gamma(beta,gamma). The Type II
%  noncentral beta distribution is the distribution of the ratio random
%  variable B2 = X1/(X1+X2), where X1 ~ Gamma(alpha,gamma) and  X2 ~
%  Gamma(beta,gamma,delta).
%
%  That is, cf_BetaNC evaluates the characteristic function cf(t) of  Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Beta(alpha_i,beta_i,delta_i) are
%  independent RVs, with the shape parameters alpha_i > 0 and beta_i >0,
%  and the noncentrality parameters delta_i > 0, for i = 1,...,N.
%
%  For Type I noncentral beta distribution, the characteristic function of
%  X ~ Beta(alpha,beta,delta) is Poisson mixture of the central Beta CFs of
%  the form
%   cf(t) = cf_BetaNC(t,alpha,beta,delta) =
%      = exp(-delta/2) sum_{j=1}^Inf (delta/2)^j/j! * cf_beta(alpha+j,beta)
%  where cf_beta(alpha+j,beta) are the CFs of central Beta RVs with
%  parameters alpha+j and beta.
%
%  For Type I noncentral beta distribution, the characteristic function of
%  X ~ Beta(alpha,beta,delta) is Poisson mixture of the central Beta CFs of
%  the form
%   cf(t) = cf_BetaNC(t,alpha,beta,delta) =
%      = exp(-delta/2) sum_{j=1}^Inf (delta/2)^j/j! * cf_beta(alpha,beta+j)
%  where cf_beta(alpha,beta+j) are the CFs of central Beta RVs with
%  parameters alpha+j and beta.
%
%  Hence,the characteristic function of Y  =
%  coef(1)*X_1 + ... + coef(N)*X_N is
%   cf(t) =  cf_X_1(coef(1)*t) * ... * cf_X_N(coef(N)*t),
%  where X_i ~ BetaNC(alpha(i),beta(i),delta(i)) with cf_X_i(t).
%
% SYNTAX
%  cf = cf_BetaNC(t,alpha,beta,delta,coef,niid,tol,type)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.
%  beta  - vector of the 'shape' parameters beta > 0. If empty, default
%          value is beta = 1.
%  delta - vector of the non-centrality parameters delta > 0. If empty,
%          default value is delta = 0.
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%  type  - indicator of the type of the noncentral distribution (Type I = 1
%          or Type II = 2). If empty, default value is type = 1. 
%
% WIKIPEDIA:
%   https://en.wikipedia.org/wiki/Noncentral_beta_distribution
%
% EXAMPLE 1:
% % CF of non-central Beta RV with delta = 1
%   alpha = 1;
%   beta  = 3;
%   delta = 1;
%   t = linspace(-50,50,201);
%   type = 1;
%   cf1 = cf_BetaNC(t,alpha,beta,delta,coef,[],[],type);
%   figure; plot(t,real(cf1),t,imag(cf1)),grid
%   title('CF of Type I and II Beta RV')
%   hold on
%   type = 2;
%   cf2 = cf_BetaNC(t,alpha,beta,delta,coef,[],[],type);
%   plot(t,real(cf2),t,imag(cf2)); hold off    
%
% EXAMPLE 2:
% % CDF/PDF of non-central Beta RV with delta = 1
%   alpha = 1;
%   beta  = 3;
%   delta = 1;
%   cf = @(t) cf_BetaNC(t,alpha,beta,delta);
%   clear options
%   options.xMin = 0;
%   options.xMax = 1;
%   result = cf2DistGP(cf,[],[],options)
%
% EXAMPLE 3:
% % CDF/PDF of the linear combination of non-central Beta RVs
%   alpha = [5 4 3];
%   beta  = [3 4 5];
%   delta = [0 1 2];
%   coef  = 1/3;
%   cf    = @(t) cf_BetaNC(t,alpha,beta,delta,coef);
%   clear options;
%   options.xMin = 0;
%   options.xMax = 1;
%   result = cf2DistGP(cf,[],[],options)

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 06-Jul-2018 15:10:17
% Rev.: 28-Apr-2020 13:47:42

%% CHECK THE INPUT PARAMETERS
narginchk(1, 8);

if nargin < 8, type = []; end
if nargin < 7, tol = []; end
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, delta = []; end
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end

if isempty(alpha), alpha = 1; end
if isempty(beta), beta = 1; end
if isempty(delta), delta = 0; end
if isempty(coef), coef = 1; end
if isempty(tol), tol = 1e-12; end
if isempty(type), type = 1; end

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,alpha,beta,delta] = ...
    distchck(4,coef(:),alpha(:),beta(:),delta(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent nc F RVs
szt = size(t);
t   = t(:);
cf  = ones(length(t),1);
for i = 1:length(coef)
    cf = cf .* cf_ncBeta(coef(i)*t,alpha(i),beta(i),delta(i),tol,type);
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
function f = cf_ncBeta(t,alpha,beta,delta,tol,type)
% cf_ncBeta Characteristic function of the non-central Beta
% distribution with the parameters alpha and beta, and the non-centrality
% parameter delta > 0.

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Feb-2018 13:53:40

%%
f = 0;
delta  = delta/2;
if delta == 0   % Deal with the central distribution
    f = cf_Beta(t,alpha,beta);
elseif delta > 0
    if type == 1 % Type I noncentral beta distribution
        % Sum the Poisson series of CFs of central iid Beta RVs,
        % poisspdf(j,delta) .* cf_Beta(t,alpha+j,beta)
        j0 = floor(delta/2);
        p0 = exp(-delta + j0 .* log(delta) - gammaln(j0 + 1));
        f  = f + p0 * cf_Beta(t,alpha+j0,beta);
        p  = p0;
        j  = j0-1;
        while j >= 0 && p > tol
            p = p * (j+1) / delta;
            f = f + p * cf_Beta(t,alpha+j,beta);
            j = j - 1;
        end
        p  = p0;
        j  = j0+1;
        i  = 0;
        while p > tol && i <= 5000
            p = p * delta / j;
            f = f + p * cf_Beta(t,alpha+j,beta);
            j = j + 1;
            i = i + 1;
        end
        if (i == 5000)
            warning(message('NoConvergence'));
        end
    else % Type II noncentral beta distribution
        % Sum the Poisson series of CFs of central iid Beta RVs,
        % poisspdf(j,delta) .* cf_Beta(t,alpha,beta+j)
        j0 = floor(delta/2);
        p0 = exp(-delta + j0 .* log(delta) - gammaln(j0 + 1));
        f  = f + p0 * cf_Beta(t,alpha,beta+j0);
        p  = p0;
        j  = j0-1;
        while j >= 0 && p > tol
            p = p * (j+1) / delta;
            f = f + p * cf_Beta(t,alpha,beta+j);
            j = j - 1;
        end
        p  = p0;
        j  = j0+1;
        i  = 0;
        while p > tol && i <= 5000
            p = p * delta / j;
            f = f + p * cf_Beta(t,alpha,beta+j);
            j = j + 1;
            i = i + 1;
        end
        if (i == 5000)
            warning(message('NoConvergence'));
        end
    end
else
    error('delta should be nonnegative');
end
end