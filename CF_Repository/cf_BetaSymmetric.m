function cf = cf_BetaSymmetric(t,theta,coef,niid)
%% cf_BetaSymmetric
%  Characteristic function of a linear combination (resp. convolution) of
%  independent zero-mean symmetric BETA random variables defined on the
%  interval (-1,1).
%
%  That is, cf_BetaSymmetric evaluates the characteristic function
%  cf(t) of  Y = sum_{i=1}^N coef_i * X_i, where X_i~BetaSymmetric(theta_i)
%  are independent RVs defined on (-1,1), for all i = 1,...,N.
%
%  The characteristic function of X ~ BetaSymmetric(theta) is defined by
%   cf(t) = cf_BetaSymmetric(t,theta) 
%          = gamma(1/2+theta) * (t/2)^(1/2-theta) * besselj(theta-1/2,t).
%
% SYNTAX
%  cf = cf_BetaSymmetric(t,theta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  theta - vector of the 'shape' parameters theta > 0. If empty, default
%          value is theta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
% 
% SPECIAL CASES:
% 1) theta = 1/2; Arcsine distribution on (-1,1):     cf(t) = besselj(0,t).
% 2) theta = 1;   Rectangular distribution on (-1,1): cf(t) = sin(t)/t;
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Beta_distribution
%
% EXAMPLE 1:
% % CF of the symmetric Beta distribution with theta = 3/2 on (-1,1))
%   theta = 3/2;
%   t = linspace(-50,50,201);
%   cf = cf_BetaSymmetric(t,theta);
%   figure; plot(t,cf),grid
%   title('CF of the symmetric Beta distribution on (-1,1)')
%
% EXAMPLE 2:
% % CF of a linear combination of independent Beta RVs
%   t = linspace(-20,20,201);
%   theta = [3 3 4 4 5]/2;
%   coef = [1 2 3 4 5]/15;
%   cf = cf_BetaSymmetric(t,theta,coef);
%   figure; plot(t,cf),grid
%   title('CF of a linear combination of independent Beta RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a weighted linear combination of independent Beta RVs
%   theta = [3 3 4 4 5]/2;
%   coef = [1 2 3 4 5]/15;
%   cf   = @(t) cf_BetaSymmetric(t,theta,coef);
%   x    = linspace(-1,1,201);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.N = 2^12;
%   options.xMin = -1;
%   options.xMax = 1;
%   result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
%   WITKOVSKY V. (2016). Numerical inversion of a characteristic
%   function: An alternative tool to form the probability distribution of
%   output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%cf = cf_BetaSymmetric(t,theta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid  = []; end
if nargin < 3, coef  = []; end
if nargin < 2, theta = []; end

if isempty(theta), theta = 1; end
if isempty(coef),  coef  = 1; end
if isempty(niid),  niid  = 1; end


[errorcode,coef,theta] = distchck(2,coef(:)',theta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt = size(t);
t   = t(:);

if length(coef)==1
    cf   = real(exp(GammaLog(0.5+theta) + (0.5-theta).*log(0.5*coef*t)) ...
        .* besselj(theta-0.5,coef*t));
else
    aux1 = t*coef;
    aux2 = ones(size(t))*theta;
    cf   = real(prod(exp(GammaLog(0.5+aux2) + (0.5-aux2).*log(0.5*aux1))...
        .* besselj(aux2-0.5,aux1),2));
end

   
cf   = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end