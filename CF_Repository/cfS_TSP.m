function cf = cfS_TSP(t,theta,mu,sigma,coef,niid)
%% cfS_TSP
%  Characteristic function of the TWO-SIDED-POWER (TSP) distribution with
%  shape parameter theta > 0.
%
%  cfS_TSP is an ALIAS of the more general function cf_TSPSymmetric, used
%  to evaluate the characteristic function of a linear combination of
%  independent (location-scale) TSP distributed random variables.
%
%  The characteristic function of the random variable X ~ TSP(theta) is 
%   cf(t) = 1/2 * (hypergeom1F1(1,1+theta,1i*t) + ...
%                  hypergeom1F1(1,1+theta,-1i*t)).
%
% SYNTAX:
%  cf = cfS_TSP(t,theta)
%  cf = cfS_TSP(t,theta,mu,sigma,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  theta - the shape parameter, theta > 0. If empty, default value is theta
%          = 1.  
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
%
% EXAMPLE 1:
% % CF of the symmetric TSP distribution with theta = 3/2 on (-1,1)
%   theta = 3/2;
%   t = linspace(-50,50,501);
%   cf = cfS_TSP(t,theta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the the symmetric TSP distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the the symmetric TSP distribution on (-1,1)
%   theta = 3/2;
%   cf = @(t) cfS_TSP(t,theta);
%   x = linspace(-1,1,101);
%   xRange = 2;
%   clear options
%   options.N = 2^8;
%   options.dt = 2*pi/xRange;
%   result = cf2DistGP(cf,x,[],options)
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
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%% ALGORITHM
narginchk(1, 6);
if nargin < 6,  niid = []; end
if nargin < 5,  coef = []; end
if nargin < 4, sigma = []; end
if nargin < 3,    mu = []; end
if nargin < 2, theta = []; end

cf = cf_TSPSymmetric(t,theta,mu,sigma,coef,niid);

end