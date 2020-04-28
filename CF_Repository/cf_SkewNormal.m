function cf = cf_SkewNormal(t,mu,sigma,lambda,coef,niid)
%cf_SkewNormal
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Skew-Normal random variables.     
%   
%  The Skew-Normal distribution is a continuous probability distribution
%  that generalises the normal distribution to allow for non-zero skewness.
%
%  The Skew-Normal distribution was first introduced by O'Hagan and Leonard
%  in [1], for more details see works of Azzalini, e.g. [2]. Pewsey gave
%  the CF of the Skew-Normal distribution and discussed the wrapped
%  Skew-Normal distribution on the circle. The CF of the Skew-Normal
%  distribution was formally proved by Kim and Genton in [4]. For more
%  details see also [5].
% 
%  In particular, if X ~ SN(mu,sigma,lambda), with real the location
%  parameter mu, the positive real scale parameter sigma, and the real
%  shape parameter lambda. 
% 
%  The probability density function (PDF) with the parameters mu, sigma and
%  lambda becomes 
%    pdf(x)= 2/sigma * phi((x-mu)/sigma) * Phi(lambda*(x-mu)/sigma), 
%  where phi denote the standard Normal (Gaussian) PDF and Phi is its CDF.
%
%  Note that the skewness of the distribution is limited to the interval
%  (-1,1). When lambda=0, the skewness vanishes, and we obtain the standard
%  Normal density, as lambda increases (in absolute value), the skewness of
%  the distribution increases, when lambda -> infty, the density converges
%  to half-normal (or folded normal) density function.
%
%  cf_SkewNormal evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ SN(mu_i,sigma_i,lambda_i) are
%  inedependent Skew-Normal RVs with the parameters mu_i, sigma_i, and
%  lambda_i, for i  = 1,...,N. 
%
%  The characteristic function of X ~ SN(mu,sigma,lambda) was defined in
%  [3] by 
%   cf(t) = cf_SkewNormal(t|mu,sigma,lambda)
%    = exp(1i*t*mu - t^2*sigma^2/2) * (1 + 1i*erfi(t*sigma*delta/sqrt(2))),
%  where delta = lambda/sqrt(1+lambda^2) and erfi(z) is the imaginary
%  error function (which is related with the Faddeeva function w(z). In
%  particular, erfi(z) = -1i*(1 - exp(z^2)*w(-z)).   
%
%  Based on that, we use the following representation 
%   cf(t) = 2*exp(1i*t*mu - t^2*sigma^2/2) - ...
%           exp( 1i*t*mu - t^2*sigma^2*(1-delta.^2)/2 ) * ...
%           Faddeeva(-t*sigma*delta/sqrt(2));
%  Note that this expression is more suitable for practical purposesis as
%  it is numerically more stable for large t. See also erfiZX.m and
%  Faddeeva.m.   
% 
%  Hence, the characteristic function of Y is  
%   cf(t) = Prod ( cf_SkewNormal(t | mu_i,sigma_i,lambda_i) ).
%
% SYNTAX:
%  cf = cf_SkewNormal(t,mu,sigma,lambda,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the 'location' parameters mu in R. If empty, default
%          value is mu = 0. 
%  sigma - vector of the scale parameters of the Half-Normal random
%          variables. If sigma is scalar, it is assumed that all scale
%          parameters are equal. If empty, default value is sigma = 1. 
%  lambda - vector of the 'shape' parameters lambda in R. If empty, default
%          value is lambda = 0. 
%  coef  - vector of the coefficients of the linear combination of the
%          Half-Normal random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Skew_normal_distribution
%
% EXAMPLE 1:
% % CF of the distribution of Skew-Normal RV with mu=0, sigma=1, lambda=1
%   mu     = 0;
%   sigma  = 1;
%   lambda = 1;
%   t      = linspace(-5,5,501);
%   cf     = cf_SkewNormal(t,mu,sigma,lambda);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the Skew-Normal RV with mu=0, sigma=1, lambda=1')
%
% EXAMPLE 2: 
% % CF of a linear combination of the Skew-Normal RVs
%   mu     = [0 0 0];
%   sigma  = [1 0.5 0.5];
%   lambda = [0 0.5 1];
%   t      = linspace(-5,5,501);
%   cf     = cf_SkewNormal(t,mu,sigma,lambda);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of the Skew-Normal RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Skew-Normal RVs
%   mu     = [0 0 0];
%   sigma  = [1 0.5 0.5];
%   lambda = [0 0.5 1];
%   cf     = @(t) cf_SkewNormal(t,mu,sigma,lambda);
%   clear options
%   options.N = 2^10;
%   x = linspace(5,40,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,[],prob,options);
%
% EXAMPLE 4:
% % PDF/CDF combination of wrapped Skew-Normal distributions on the circle
%   mu     = [0 0 0];
%   sigma  = [1 0.5 0.5];
%   lambda = [0 0.5 1];
%   cf     = @(t) cf_SkewNormal(t,mu,sigma,lambda);
%   clear options
%   options.isCircular = true;
%   options.correctedCDF = true;
%   result = cf2DistGP(cf,[],[],options);
%   disp(result);
%   angle  = result.x;
%   radius = result.pdf;
%   figure; polarplot(angle,radius);
%   ax = gca; ax.ThetaAxisUnits = 'radians';
%   title('PDF of a Wrapped Skew-Normal Distribution on the Circle')
%
% REFERENCES:
% [1] O'Hagan, A. and Leonard, T., 1976. Bayes estimation subject to
%     uncertainty about parameter constraints. Biometrika, 63(1), 201-203.  
% [2] Azzalini, A., 1985. A class of distributions which includes the
%     normal ones. Scandinavian Journal of Statistics, 171-178.
% [3] Pewsey, A., 2000. The wrapped skew-normal distribution on the circle.
%     Communications in Statistics-Theory and Methods, 29(11), 2459-2472.
% [4] Kim, H.M. and Genton, M.G., 2011. Characteristic functions of scale
%     mixtures of multivariate skew-normal distributions. Journal of
%     Multivariate Analysis, 102(7), 1105-1117.  
% [5] Gupta, A.K., Nguyen, T.T. and Sanqui, J.A.T., 2004. Characterization
%     of the skew-normal distribution. Annals of the Institute of
%     Statistical Mathematics, 56(2), 351-360.  
%
% SEE ALSO: cf_Normal, erfiZX

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Feb-2020 13:08:06
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_SkewNormal(t,mu,sigma,lambda,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid  = []; end
if nargin < 5, coef  = []; end
if nargin < 4, lambda = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

if isempty(lambda)
    lambda = 0;
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

%% CHECK THE INPUT PARAMETERS
[errorcode,coef,mu,sigma,lambda] = ...
    distchck(3,coef(:)',mu(:)',sigma(:)',lambda(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end
szt = size(t);
t   = t(:);

%% CF of the linear combination of the Skew-Normal RVs
delta = lambda ./ sqrt(1+lambda.^2);
% cf =  prod(exp(1i*t*(coef.*mu) - t.^2*(coef.*sigma).^2/2) .* ...
%     (1 + 1i*erfiZX(t*(sigma.*delta)/sqrt(2))),2);

% CF based on using Faddeeva function (numerically more stable)
aux1 = 1i*t*(coef.*mu);
aux2 = (coef.*sigma).^2/2;
cf =  prod(2 * exp(aux1 - t.^2*aux2) - ...
    exp(aux1 - t.^2*(aux2.*(1-delta.^2))) .* ...
    Faddeeva(-t*(coef.*sigma.*delta)/sqrt(2)),2);

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