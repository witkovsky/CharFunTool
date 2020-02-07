function cf = cf_SkewNormal(t,mu,sigma,alpha,coef,niid)
%cf_SkewNormal
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Skew-Normal random variables.     
%   
%  The Skew-Normal distribution is a continuous probability distribution
%  that generalises the normal distribution to allow for non-zero skewness.
%  This distribution was first introduced by O'Hagan and Leonard (1976),
%  see also Azzalini (1985). 
% 
%  In particular, if X ~ SN(mu,sigma,alpha), with real the location
%  parameter mu, the positive real scale parameter sigma, and the real
%  shape parameter alpha. 
% 
%  The probability density function (PDF) with the parameters mu, sigma and
%  alpha becomes 
%    pdf(x)= 2/sigma * phi((x-mu)/sigma) * Phi(alpha*(x-mu)/sigma), 
%  where phi denote the standard Normal (Gaussian) PDF and Phi is its CDF.
%
%  Note that the skewness of the distribution is limited to the interval
%  (-1,1). When alpha=0, the skewness vanishes, and we obtain the standard
%  Normal density, as alpha increases (in absolute value), the skewness of
%  the distribution increases, when alpha -> infty, the density converges
%  to half-normal (or folded normal) density function.
%
%  cf_SkewNormal evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ SN(mu_i,sigma_i,alpha_i) are
%  inedependent Skew-Normal RVs with the parameters mu_i, sigma_i, and
%  alpha_i, for i  = 1,...,N. 
%
%  The characteristic function of X ~ SN(mu,sigma,alpha) is defined by
%   cf(t) = cf_SkewNormal(t|mu,sigma,alpha)
%    = exp(1i*t*mu - t^2*sigma^2/2) * (1 + 1i*erfi(t*sigma*delta/sqrt(2))),
%  where delta = alpha/sqrt(1+alpha^2) and erfi(z) is the imaginary
%  error function (which is related with the Fadeeva function w(z). In
%  particular, erfi(z) = -1i*(1 - exp(z^2)*w(-z)).   
%
%  Based on that, we use the following representation 
%   cf(t) = 2*exp(1i*t*mu - t^2*sigma^2/2) - ...
%           exp( 1i*t*mu - t^2*sigma^2*(1-delta.^2)/2 ) * ...
%           Fadeeva(-t*sigma*delta/sqrt(2));
%  Note that this expression is more suitable for practical purposesis as
%  it is numerically more stable for large t. See also erfiZX.m and
%  Fadeeva.m.   
% 
%  Hence, the characteristic function of Y is  
%   cf(t) = Prod ( cf_SkewNormal(t | mu_i,sigma_i,alpha_i) ).
%
% SYNTAX:
%  cf = cf_SkewNormal(t,mu,sigma,alpha,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the 'location' parameters mu in R. If empty, default
%          value is mu = 0. 
%  sigma - vector of the scale parameters of the Half-Normal random
%          variables. If sigma is scalar, it is assumed that all scale
%          parameters are equal. If empty, default value is sigma = 1. 
%  alpha - vector of the 'shape' parameters alpha in R. If empty, default
%          value is alpha = 0. 
%  coef  - vector of the coefficients of the linear combination of the
%          Half-Normal random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Skew_normal_distribution
%
% EXAMPLE 1:
% % CF of the distribution of Skew-Normal RV with mu=0, sigma=1, alpha=1
%   mu    = 0;
%   sigma = 1;
%   alpha = 1;
%   t     = linspace(-5,5,501);
%   cf    = cf_SkewNormal(t,mu,sigma,alpha);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the Skew-Normal RV with mu=0, sigma=1, alpha=1')
%
% EXAMPLE 2: 
% % CF of a linear combination of the Skew-Normal RVs
%   mu    = [0 0 0];
%   sigma = [1 0.5 0.5];
%   alpha = [0 0.5 1];
%   t     = linspace(-5,5,501);
%   cf    = cf_SkewNormal(t,mu,sigma,alpha);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of the Skew-Normal RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent  Skew-Normal RVs
%   mu    = [0 0 0];
%   sigma = [1 0.5 0.5];
%   alpha = [0 0.5 1];
%   cf    = @(t) cf_SkewNormal(t,mu,sigma,alpha);
%   clear options
%   options.N = 2^10;
%   x = linspace(5,40,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,[],prob,options);
%
% See also: cf_Normal, erfiZX

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Feb-2020 13:08:06
%% ALGORITHM
%  cf = cf_SkewNormal(t,mu,sigma,alpha,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid  = []; end
if nargin < 5, coef  = []; end
if nargin < 4, alpha = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

if isempty(alpha)
    alpha = 0;
end

if isempty(sigma)
    sigma = 1;
end

if isempty(mu)
    mu = 0;
end

if isempty(coef) 
    coef = 1;
end

%%
[errorcode,coef,mu,sigma,alpha] = ...
    distchck(3,coef(:)',mu(:)',sigma(:)',alpha(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end
szt = size(t);
t   = t(:);

% CF of the linear combination of the Skew-Normal RVs
delta = alpha ./ sqrt(1+alpha.^2);
% cf =  prod(exp(1i*t*(coef.*mu) - t.^2*(coef.*sigma).^2/2) .* ...
%     (1 + 1i*erfiZX(t*(sigma.*delta)/sqrt(2))),2);

% CF based on using Fadeeva function (numerically more stable)
aux1 = 1i*t*(coef.*mu);
aux2 = (coef.*sigma).^2/2;
cf =  prod(2 * exp(aux1 - t.^2*aux2) - ...
    exp(aux1 - t.^2*(aux2.*(1-delta.^2))) .* ...
    Fadeeva(-t*(coef.*sigma.*delta)/sqrt(2)),2);

cf(t==0) = 1;
cf   = reshape(cf,szt);

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end
end