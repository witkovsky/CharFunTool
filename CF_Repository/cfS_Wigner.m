function cf = cfS_Wigner(t,coef,niid)
%% cfS_Wigner
%  Characteristic function of a linear combination (resp. convolution) of
%  the zero-mean symmetric WIGNER distribution defined on the interval
%  (-1,1).  
%
%  cfS_Wigner is an ALIAS of the more general function cf_WignerSemicircle,
%  used to evaluate the characteristic function of a linear combination of
%  independent symmetric WIGNER SEMICIRCLE distributed random variables. 
%
%  The characteristic function of the symmetric WIGNER distribution on
%  (-1,1) is 
%   cf(t) = 2*besselj(1,t)/t.
%
% SYNTAX:
%  cf = cfS_Wigner(t);
%  cf = cfS_Wigner(t,coef,niid);
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1. 
% 
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
%
% EXAMPLE 1:
% % CF of the symmetric WIGNER distribution on (-1,1)
%   t = linspace(-50,50,501);
%   cf = cfS_Wigner(t);
%   figure; plot(t,cf),grid
%   title('CF of the the WIGNER distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of the symmetric WIGNER distribution on (-1,1)
%   cf   = @(t) cfS_Wigner(t);
%   x    = linspace(-1,1,501);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = -1;
%   options.xMax = 1;
%   result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
%  WITKOVSKY V. (2016). Numerical inversion of a characteristic
%  function: An alternative tool to form the probability distribution of
%  output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Sep-2018 00:45:54
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 3);
if nargin < 3, niid = []; end
if nargin < 2, coef = []; end

cf = cf_WignerSemicircle(t,[],[],coef,niid);

end