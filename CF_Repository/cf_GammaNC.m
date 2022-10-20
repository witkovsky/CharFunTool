function cf = cf_GammaNC(t,alpha,beta,delta,coef,niid,tol)
%% cf_GammaNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Gamma random variables with distributions
%  Gamma(alpha_i,beta_i,delta_i), where the shape parameters alpha_i > 0
%  and the rate parameters beta_i > 0, and the non-centrality parameters
%  delta_i >= 0 for i = 1,...,N.
%
%  That is, cf_GammaNC evaluates the characteristic function cf(t)
%  of Y = coef_i*X_1 +...+ coef_N*X_N, where X_i ~
%  Gamma(alpha_i,beta_i,delta_i) are inedependent RVs, with parameters
%  alpha_i > 0, beta_i > 0, and delta_i >= 0, for i = 1,...,N.   
%
%  The characteristic function of Y = X with X ~ Gamma(alpha,beta,delta) is
%  Poisson mixture of CFs of the central Gamma RVs of the form
%   cf(t) = cf_GammaNC(t,alpha,beta,delta) = 
%         = exp(-delta/2) sum_{j=1}^Inf (delta/2)^j/j! *
%           * cf_Gamma(alpha+j,beta) 
%  where cf_Gamma(alpha+j,beta) are the CFs of central Gamma RVs with
%  parameters alpha+j and beta. 
%
% SYNTAX
%  cf = cf_GammaNC(t,alpha,beta,delta,coef,niid,tol)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'rate' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0.
%  coef  - vector of the coefficients of the linear combination of the
%          gamma distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Gamma_distribution
%
% EXAMPLE 1:
% % CF of the non-central Gamma RV with delta = 1
%   alpha = 3/2;
%   beta  = 1/2;
%   delta = 1;
%   t = linspace(-10,10,201);
%   cf = cf_GammaNC(t,alpha,beta,delta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the non-central Gamma RV with delta = 1')
%
% EXAMPLE 2:
% % CDF/PDF of the non-central Gamma RV
%   alpha = 3/2;
%   beta  = 1/2;
%   delta = 1;
%   cf = @(t) cf_GammaNC(t,alpha,beta,delta);
%   clear options;
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options)
%
% EXAMPLE 3:
% % CDF/PDF of the linear combination of the non-central Gamma RVs 
%   alpha = [3 4 5]/2;
%   beta  = [1 1 2]/2;
%   delta = [0 1 2];
%   coef  = [1 1 1]/3;
%   cf = @(t) cf_GammaNC(t,alpha,beta,delta,coef);
%   clear options;
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Oct-2022 14:31:48

%% CHECK THE INPUT PARAMETERS
narginchk(1, 7);

if nargin < 7, tol   = []; end
if nargin < 6, niid  = []; end
if nargin < 5, coef  = []; end
if nargin < 4, delta = []; end
if nargin < 3, beta  = []; end
if nargin < 2, alpha = []; end

if isempty(alpha), alpha = 1; end
if isempty(beta), beta = 1; end
if isempty(delta), delta = 0; end
if isempty(coef), coef = 1; end
if isempty(tol), tol = 1e-12; end

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,alpha,beta,delta] = ...
    distchck(4,coef(:),alpha(:),beta(:),delta(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent nc RVs
szt = size(t);
t   = t(:);
cf  = ones(length(t),1);
for i = 1:length(coef)
    cf = cf .* cf_ncGamma(coef(i)*t,alpha(i),beta(i),delta(i),tol);
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
function f = cf_ncGamma(t,alpha,beta,delta,tol)
% cf_ncGamma Characteristic function of the distribution of the
% non-central Gamma RV with parameters alpha > 0, beta > 0, 
% and the non-centrality parameter delta >= 0. 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Oct-2022 14:31:48

%%
f = 0;
delta  = delta/2;
if delta == 0   % Deal with the central distribution
    f = cf_Gamma(t,alpha,beta);
elseif delta > 0
    % Sum the Poisson series of CFs of central iid Gamma RVs,
    % poisspdf(j,delta) .* cf_Gamma(t,alpha+j,beta)
    j0 = floor(delta/2);
    p0 = exp(-delta + j0 .* log(delta) - gammaln(j0 + 1));
    f  = f + p0 * cf_Gamma(t,alpha+j0,beta);
    p  = p0;
    j  = j0-1;
    while j >= 0 && p > tol
        p = p * (j+1) / delta;
        f = f + p * cf_Gamma(t,alpha+j,beta);
        j = j - 1;
    end
    p  = p0;
    j  = j0+1;
    i  = 0;
    while p > tol && i <= 5000
        p = p * delta / j;
        f = f + p * cf_Gamma(t,alpha+j,beta);
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