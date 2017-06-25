function cf = cf_TSPSymmetric(t,theta,mu,sigma,coef,niid)
%% cf_TSPSymmetric
%  Characteristic function of a linear combination (resp. convolution) of
%  independent symmetric (location and scale shifted) TWO-SIDED-POWER (TSP)
%  random variables. 
%
%  That is, cfS_TSPSymmetric evaluates the characteristic function cf(t) of
%  Y = sum_{i=1}^N coef_i * (mu_i + sigma_i * X_i), where X_i ~
%  TSP(theta_i) are inedependent RVs, with symmetric TSP distributions
%  defined on the interval (-1,1) with zero mean and variance Var(X_i) =
%  2*theta_i*gamma(theta_i)/gamma(3+theta_i), where theta_i > 0 are shape
%  parameters for all i = 1,...,N. 
%
%  The characteristic function of the random variable mu + sigma*X, where
%  X ~ TSP(theta) is given by 
%   cf(t) = cfS_TSPSymmetric(t,theta,mu,sigma) 
%         = 1/2 * exp(1i*t*mu) * (hypergeom1F1(1,1+theta,1i*t*sigma) + ...
%           hypergeom1F1(1,1+theta,-1i*t*sigma)). 
%
%  Hence, the characteristic function of Y  = coef_1*(mu_1+sigma_1*X_1)
%  +...+ coef_N*(mu_N+sigma_N*X_N) is cf_Y(t) = exp(1i*mu*t) *
%  (cf_1(coef_1*sigma_1*t) *...* cf_N(coef_N*sigma_N*t)), where cf_i(t) is
%  the characteristic function of X_i ~ TSP(theta_i).
%
% SYNTAX:
%  cf = cf_TSPSymmetric(t,theta,mu,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  theta - vector of the shape parameters, theta > 0. If theta is scalar,
%          it is assumed that all parameters theta are equal. If empty,
%          default value is theta = 1.  
%  mu    - vector of location parameters, mu in Real. If empty, default
%          value is mu = 0. 
%  sigma - vector of scale parameters, sigma_i > 0. If empty, default value
%          is sigma = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef_i * X_i is independently and identically distributed random
%          variable. If empty, default value is niid = 1. 
%
% SPECIAL CASES:
%  - theta = 1/2; Arcsine distribution on (-1,1):     cf(t) = besselj(0,t).
%  - theta = 1;   Rectangular distribution on (-1,1): cf(t) = sin(t)/t;
%
% EXAMPLE 1: 
% % CF of the symmetric TSP distribution with theta = 3/2 on (-1,1)
%   theta = 3/2;
%   t = linspace(-50,50,501);
%   cf = cf_TSPSymmetric(t,theta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the the symmetric TSP distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the the symmetric TSP distribution on (-1,1)
%   theta = 3/2;
%   cf = @(t) cf_TSPSymmetric(t,theta);
%   x = linspace(-1,1,101);
%   xRange = 2;
%   clear options
%   options.N = 2^8;
%   options.dt = 2*pi/xRange;
%   result = cf2DistGP(cf,x,[],options)
%
% EXAMPLE 3: 
% % CF of the weighted linear combinantion of TSP RVs
%   theta = [1 2 3 4 5]/2;
%   mu    = [1 2 0 0 0];
%   sigma = [1 2 3 4 5]/5;
%   coef  = 1/5;
%   t = linspace(-50,50,501);
%   cf = cf_TSPSymmetric(t,theta,mu,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the weighted linear combinantion of TSP RVs')
%
% EXAMPLE 4: 
% % CDF/PDF of the weighted linear combinantion of TSP RVs
%   theta = [1 2 3 4 5]/2;
%   mu    = 0;
%   sigma = [5 4 3 2 1];
%   coef  = 1/5;
%   t = linspace(-50,50,501);
%   cf = @(t) cf_TSPSymmetric(t,theta,mu,sigma,coef);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
% [2] VAN DORP, R.J., KOTZ, S. (2003). Generalizations of two-sided power
%     distributions and their convolution. Communications in
%     Statistics-Theory and Methods, 32(9), 1703-1723. 

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jun-2017 18:25:56

%% ALGORITHM
%cf = cf_TSPSymmetric(t,theta,mu,sigma,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6,  niid = []; end
if nargin < 5,  coef = []; end
if nargin < 4, sigma = []; end
if nargin < 3,    mu = []; end
if nargin < 2, theta = []; end

%%
if isempty(theta), theta = 1; end
if isempty(coef),   coef = 1; end
if isempty(mu),       mu = 0; end
if isempty(sigma), sigma = 1; end
if isempty(niid),   niid = 1; end

%% Equal size of the parameters
[errorcode,coef,theta,mu,sigma] = ...
    distchck(4,coef(:)',theta(:)',mu(:)',sigma(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt = size(t);
t   = abs(t(:));

cf  = 1;
for i = 1:length(coef)
    cf = cf .* exp(1i*t*mu(i)) .* ...
        ( Hypergeom1F1(1,1+theta(i),1i*t*coef(i)*sigma(i)) + ...
        Hypergeom1F1(1,1+theta(i),-1i*t*coef(i)*sigma(i)) )/2;
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