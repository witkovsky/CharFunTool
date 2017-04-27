function cf = cfS_TSP(t,theta)
%cfS_TSP(t,theta) evaluates the characteristic function cf(t) of
% the symmetric zero-mean Two-Sided-Power (TSP) distribution with shape
% parameter theta > 0, defined on the interval (-1,1), i.e. symmetric TSP
% distribution with MEAN = 0, and variance VAR =
% 2*theta*gamma(theta)/gamma(3+theta). The standard deviation is given by
% STD = sqrt(VAR), i.e.  
%   cf(t) = cfS_TSP(t,theta) 
%         = (theta/2)*t^(-2*theta)* (exp(-1i*t))*(1i*t)^theta *
%         (gamma(theta)-gammainc(theta,-1i*t,'upper')
%         +exp(1i*t))*(-1i*t)^theta *
%         (gamma(theta)-gammainc(theta,1i*t,'upper'))),   
% where gammainc(theta,1i*t,'upper') = exp(-1i*t)*(1i*t)^theta
% *hypergeomU(1,1+theta,1i*t).
%
% Special cases (for specific values of the shape parameter theta):
% 1) theta = 1/2; Arcsine distribution on (-1,1):     cf(t) = besselj(0,t).
% 2) theta = 1;   Rectangular distribution on (-1,1): cf(t) = sin(t)/t;
% For more details see:
% van Dorp, R.J., Kotz, S. (2003). Generalizations of two-sided power
% distributions and their convolution. Communications in Statistics-Theory
% and Methods, 32(9), 1703-1723. 
%
% SYNTAX
%  cf = cfS_TSP(t,theta)
%
% EXAMPLE1 (CF of the symmetric TSP distribution with theta = 3/2 on (-1,1))
%  theta = 3/2;
%  t = linspace(-50,50,501);
%  cf = cfS_TSP(t,theta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the the symmetric TSP distribution on (-1,1)')
%
% EXAMPLE2 (PDF/CDF of the the symmetric TSP distribution on (-1,1))
%  theta = 3/2;
%  cf = @(t) cfS_TSP(t,theta);
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
% [2] van Dorp, R.J., Kotz, S. (2003). Generalizations of two-sided power
%     distributions and their convolution. Communications in
%     Statistics-Theory and Methods, 32(9), 1703-1723. 

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Feb-2017 14:58:48

%% ALGORITHM
%cf = cfS_TSP(t,theta);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, theta = []; end
if isempty(theta), theta = 1; end


%% Characteristic function of the symmetric TSP distribution
szt = size(t);
t   = abs(t(:));

% cf = (theta/2) * ((t.^(-2*theta)).* ...
%      (exp(-1i*t).*(1i*t).^theta + exp(1i*t).*(-1i*t).^theta) *gamma(theta) ...
%       - (KummerU(1,1+theta,1i*t) + KummerU(1,1+theta,-1i*t)));
  cf = (1/2) * ((hypergeom1F1(1,1+theta,1i*t) + hypergeom1F1(1,1+theta,-1i*t)));
%cf = min(1,cf);
cf = reshape(cf,szt);
cf(t==0) = 1;

end