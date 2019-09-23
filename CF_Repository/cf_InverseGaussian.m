function cf = cf_InverseGaussian(t,mu,lambda,coef,niid)
%% cf_InverseGaussian 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent INVERSE-GAUSSIAN random variables.
%
%  The Inverse-Gaussian distribution (also known as the Wald distribution)
%  is a two-parameter family of continuous probability distributions with
%  support on (0,Infinity), with the parameter mu > 0 (the mean parameter)
%  and the parameter lambda > 0 (the shape parameter).
%
%  As lambda tends to infinity, the Inverse-Gaussian distribution becomes
%  more like a normal (Gaussian) distribution. The inverse Gaussian
%  distribution has several properties analogous to a Gaussian
%  distribution, however the name can be misleading.
%
%  cf_InverseGaussian evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ InverseGaussian(mu_i,lambda_i) are
%  inedependent (symmetric) Inverse-Gaussian distributed RVs, with mu_i > 0
%  and lambda_i > 0, for i = 1,...,N.
%
%  The characteristic function of X ~ InverseGaussian(mu,lambda) is given
%  by  
%   cf(t) = exp(lambda/mu *(1 - sqrt(1-2i*mu^2*t/lambda))).
%  Hence, the characteristic function of Y is
%   cf(t) = Prod ( cf_InverseGaussian(coef_i*t,mu_i,lambda_i) ).
%
% SYNTAX:
%  cf = cf_InverseGaussian(t,mu,lambda,coef,niid)
% 
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  mu     - vector of mean parameters, mu > 0. If empty, default
%           value is mu = 1. 
%  lambda - vector of shape parameters, lambda_i > 0. If empty, default
%           value is lambda = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef_i * X_i is independently and identically distributed random
%          variable. If empty, default value is niid = 1. 
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
%
% EXAMPLE 1:
%  % CF of Inverse-Gaussian RVs with mu = 1 and lambda = 3
%  mu     = 1;
%  lambda = 3;
%  t      = linspace(-20,20,201);
%  cf     = cf_InverseGaussian(t,mu,lambda);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('CF of Inverse-Gaussian RVs with mu = 1 and lambda = 3')
%
% EXAMPLE 2:
%  % CDF/PDF of Inverse-Gaussian RVs with mu = 1 and lambda = 3
%  mu     = 1;
%  lambda = 3;
%  cf     = @(t) cf_InverseGaussian(t,mu,lambda);
%  x      = linspace(0,4);
%  prob   = [0.9 0.95 0.975 0.99];
%  clear options;
%  options.N = 2^12;
%  options.xMin = 0;
%  result = cf2DistGP(cf,x,prob,options);
%  disp(result)
%
% EXAMPLE 3:
%  % CF of a linear combination of independent Inverse-Gaussian RVs
%  mu     = [1 1 2 2 3];
%  lambda = [1 2 3 4 5];
%  t      = linspace(-2,2,201);
%  cf     = cf_InverseGaussian(t,mu,lambda);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('CF of a linear combination of Inverse-Gaussian RVs')
%
% EXAMPLE 4:
%  % CDF/PDF of a linear combination of independent Inverse-Gaussian RVs
%  mu     = [1 1 2 2 3];
%  lambda = [1 2 3 4 5];
%  cf   = @(t) cf_InverseGaussian(t,mu,lambda);
%  x    = linspace(0,30);
%  prob = [0.9 0.95 0.975 0.99];
%  clear options;
%  options.N = 2^12;
%  options.xMin = 0;
%  result = cf2DistGP(cf,x,prob,options);
%  disp(result)

% (c) 2019 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Sep-2019 17:58:27

%% ALGORITHM
% cf = cf_InverseGaussian(t,df,mu,sigma,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5,   niid = []; end
if nargin < 4,   coef = []; end
if nargin < 3, lambda = []; end
if nargin < 2,     mu = []; end

if isempty(mu), mu = 1; end
if isempty(lambda), lambda = 1; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Equal size of the parameters
[errorcode,coef,mu,lambda] = ...
    distchck(3,coef(:)',mu(:)',lambda(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function 
szt   = size(t);
t     = t(:);
cf    = exp(sum(bsxfun(@times,(1-sqrt(1-2i*t*(coef.*mu.^2./lambda))),...
        lambda./mu),2));
cf    = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end
