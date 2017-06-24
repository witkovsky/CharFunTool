function cf = cfS_Beta(t,theta)
%% cfS_Beta
%  Characteristic function of the zero-mean symmetric BETA distribution
%  defined on the interval (-1,1). 
%
%  cfS_Beta is an ALIAS NAME of the more general function cf_BetaSymmetric,
%  used to evaluate the characteristic function of a linear combination of
%  independent BETA distributed random variables.
%
%  The characteristic function of X ~ BetaSymmetric(theta) is defined by
%   cf(t) = cf_BetaSymmetric(t,theta) 
%          = gamma(1/2+theta) * (t/2)^(1/2-theta) * besselj(theta-1/2,t).
%
% SYNTAX:
%  cf = cfS_Beta(t,theta);
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  theta - vector of the 'shape' parameters theta > 0. If empty, default
%          value is theta = 1. 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Beta_distribution
%
% EXAMPLE 1:
% % CF of the symmetric Beta distribution with theta = 3/2 on (-1,1)
%   theta = 3/2;
%   t = linspace(-50,50,501);
%   cf = cfS_Beta(t,theta);
%   figure; plot(t,cf),grid
%   title('CF of the the symmetric Beta distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the the symmetric Beta distribution on (-1,1)
%   theta = 3/2;
%   cf = @(t) cfS_Beta(t,theta);
%   x = linspace(-1,1,101);
%   xRange = 2;
%   clear options
%   options.N = 2^8;
%   options.dt = 2*pi/xRange;
%   result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
%   WITKOVSKY V. (2016). Numerical inversion of a characteristic
%   function: An alternative tool to form the probability distribution of
%   output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24

%% ALGORITHM
narginchk(1, 2);
if nargin < 2, theta  = []; end
if isempty(theta), theta = 1; end

cf = cf_BetaSymmetric(t,theta);

end