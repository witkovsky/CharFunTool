%function cf = cfS_Arcsine(t)
%cfS_Arcsine(t) evaluates the characteristic function cf(t) of
% the symmetric zero-mean Arcsine distribution on the interval
% (-1,1)(U-shaped distribution with zero mean and variance 1/2), i.e.   
%   cf(t) = cfS_Arcsine(t) = besselj(0,t);
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Arcsine_distribution
%
% SYNTAX
%  cf = cfS_Arcsine(t)
%
% EXAMPLE1 (CF of the symmetric Arcsine distribution on (-1,1))
%  t = linspace(-50,50,501);
%  cf = cfS_Arcsine(t);
%  figure; plot(t,cf),grid
%  title('CF of the the Arcsine distribution on (-1,1)')
%
% EXAMPLE2 (PDF/CDF of the symmetric Arcsine distribution on (-1,1))
%  cf = @(t) cfS_Arcsine(t);
%  x = linspace(-1,1,501);
%  xRange = 2;
%  clear options
%  options.N = 2^12;
%  options.dt = 2*pi/xRange;
%  result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfS_Arcsine(t);