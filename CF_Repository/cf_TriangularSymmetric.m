function cf = cf_TriangularSymmetric(t,coef,niid)
%% cf_TriangularSymmetric
%  Characteristic function of a linear combination (resp. convolution) of
%  independent zero-mean symmetric TRIANGULAR random variables defined on
%  the interval (-1,1).
%
%  That is, cf_TriangularSymmetric evaluates the characteristic function
%  cf(t) of  Y = sum_{i=1}^N coef_i * X_i, where X_i ~ TriangularSymmetric
%  are independent uniformly distributed RVs defined on (-1,1), for all i =
%  1,...,N. 
%
%  The characteristic function of X ~ TriangularSymmetric is defined by
%   cf(t) = cf_TriangularSymmetric(t) = (2-2*cos(t))/t^2;
%
% SYNTAX
%  cf = cf_TriangularSymmetric(t,coef,niid)
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
%  https://en.wikipedia.org/wiki/Triangular_distribution
%
% EXAMPLE 1:
% % CF of the Triangular distribution on (-1,1)
%   t = linspace(-50,50,201);
%   cf = cf_TriangularSymmetric(t);
%   figure; plot(t,cf),grid
%   title('CF of the Triangular distribution on (-1,1)')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of independent Triangular RVs
%   t = linspace(-20,20,201);
%   coef = [1 2 3 4 5]/15;
%   cf = cf_TriangularSymmetric(t,coef);
%   figure; plot(t,cf),grid
%   title('CF of a weighted linear combination of Triangular RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a weighted linear combination of independent Triangular RVs
%   coef = [1 2 3 4 5]/15;
%   cf   = @(t) cf_TriangularSymmetric(t,coef);
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
%cf = cf_TriangularSymmetric(t,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, niid = []; end
if nargin < 2, coef = []; end

if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Characteristic function
szt  = size(t);
t    = t(:);
coef = coef(:)';
cf   = prod((2-2*cos(t*coef))./(t*coef).^2,2);

cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end