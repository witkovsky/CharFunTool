function cf = cfS_Beta(t,theta)
%cfS_Beta(t,theta) evaluates the characteristic function cf(t) of
% the symmetric zero-mean Beta distribution with shape parameter theta >0,
% defined on the interval (-1,1), i.e. symmetric beta distribution with
% zero mean and variance VAR = 1/(1+2*theta). The standard deviation is
% given by STD = sqrt(1/(1+2*theta))), i.e. 
%   cf(t) = cfS_Beta(t,theta) 
%         = gamma(1/2+theta) * (t/2)^(1/2-theta) * besselj(theta-1/2,t).
% Special cases (for specific values of the shape parameter theta):
% 1) theta = 1/2; Arcsine distribution on (-1,1):     cf(t) = besselj(0,t).
% 2) theta = 1;   Rectangular distribution on (-1,1): cf(t) = sin(t)/t;
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Beta_distribution
%
% SYNTAX
%  cf = cfS_Beta(t,theta)
%
% EXAMPLE1 (CF of the symmetric Beta distribution with theta = 3/2 on (-1,1))
%  theta = 3/2;
%  t = linspace(-50,50,501);
%  cf = cfS_Beta(t,theta);
%  figure; plot(t,cf),grid
%  title('CF of the the symmetric Beta distribution on (-1,1)')
%
% EXAMPLE2 (PDF/CDF of the the symmetric Beta distribution on (-1,1))
%  theta = 3/2;
%  cf = @(t) cfS_Beta(t,theta);
%  x = linspace(-1,1,101);
%  xRange = 2;
%  clear options
%  options.N = 2^8;
%  options.dt = 2*pi/xRange;
%  result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Feb-2017 14:58:48

%% ALGORITHM
%cf = cfS_Beta(t,theta);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, theta = []; end
if isempty(theta), theta = 1; end


%% Characteristic function of the symmetric Beta distribution
szt = size(t);
t   = t(:);

cf = min(1,gamma(0.5+theta)*(0.5*t).^(0.5-theta).*besselj(theta-0.5,t));
cf = reshape(cf,szt);
cf(t==0) = 1;

end