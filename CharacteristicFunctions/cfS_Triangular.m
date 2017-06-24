function cf = cfS_Triangular(t)
%cfS_Triangular(t)

%  Characteristic function of the zero-mean symmetric TRIANGULAR
%  distribution defined on the interval (-1,1). 
%
%  cfS_Triangular is a version resp. ALIAS of the more the more general
%  function cf_TriangularSymmetric, used to evaluate the characteristic
%  function of a linear combination of independent TRIANGULAR distributed
%  random variables.
%
%  The characteristic function of X ~ TriangularSymmetric is defined by
%   cf(t) = (2-2*cos(t))/t^2;;
%
% SYNTAX:
%  cf = cfS_Triangular(t);
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Triangular_distribution
%
% EXAMPLE 1:
% % CF of the symmetric Triangular distribution on (-1,1)
%   t = linspace(-50,50,501);
%   cf = cfS_Triangular(t);
%   figure; plot(t,cf),grid
%   title('CF of the the symmetric Triangular distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the the symmetric Triangular distribution on (-1,1)
%   cf = @(t) cfS_Triangular(t);
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
cf = cf_TriangularSymmetric(t);

end