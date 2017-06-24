function cf = cfS_Rectangular(t)
%% cfS_Rectangular
%  Characteristic function of the zero-mean symmetric RECTANGULAR
%  distribution defined on the interval (-1,1). 
%
%  cfS_Rectangular is a version resp. ALIAS of the more the more general
%  function cf_ArcsineSymmetric, used to evaluate the characteristic
%  function of a linear combination of independent RECTANGULAR distributed
%  random variables.
%
%  The characteristic function of X ~ RectangularSymmetric is defined by
%   cf(t) = sinc(t) = sin(t)/t;
%
% SYNTAX:
%  cf = cfS_Rectangular(t);
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
%
% EXAMPLE 1:
% % CF of the Rectangular distribution on (-1,1)
%   t = linspace(-50,50,501);
%   cf = cfS_Rectangular(t);
%   figure; plot(t,cf),grid
%   title('CF of the Rectangular distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the Rectangular distribution on (-1,1)
%   cf = @(t) cfS_Rectangular(t);
%   x = linspace(-2,2,101);
%   xRange = 2;
%   clear options
%   options.N = 2^5;
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
cf = cf_RectangularSymmetric(t);

end