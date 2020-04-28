function cf = cf_TrapezoidalSymmetric(t,lambda,coef,niid)
%% cf_TrapezoidalSymmetric
%  Characteristic function of a linear combination (resp. convolution) of
%  independent zero-mean symmetric TRAPEZOIDAL random variables defined on
%  the interval (-1,1).
%
%  That is, cf_TrapezoidalSymmetric evaluates the characteristic function
%  cf(t) of  Y = sum_{i=1}^N coef_i * X_i, where X_i ~
%  TrapezoidalSymmetric(lambda_i) are independent RVs defined on (-1,1),
%  for all i = 1,...,N.
%
%  The characteristic function of X ~ TrapezoidalSymmetric(lambda), with
%  E(X) = 0 and Var(X) = (1+lambda^2)/6, is  
%   cf(t) = cf_TrapezoidalSymmetric(t) 
%         = cf_RectangularSymmetric(w*t))*cf_RectangularSymmetric((1-w)*t);
%         = (sin(w*t)./(w*t)).*(sin((1-w)*t)./((1-w)*t)),
%  where w  = (1+lambda)/2.
%
% SYNTAX
%  cf = cf_TrapezoidalSymmetric(t,lambda,coef,niid)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  lambda - parameter of the offset, 0 <= lambda <=1. If empty, default
%           value is lambda = 0.  
%  coef   - vector of the coefficients of the linear combination of the
%           Beta distributed random variables. If coef is scalar, it is
%           assumed that all coefficients are equal. If empty, default value
%           is coef = 1.
%  niid   - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%           sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%           coef(i) * X_i is independently and identically distributed
%           random variable. If empty, default value is niid = 1.   
%
% EXAMPLE 1:
% % CF of the symmetric Trapezoidal distribution with lambda = 1/2 
%   lambda = 1/2;
%   t = linspace(-50,50,201);
%   cf = cf_TrapezoidalSymmetric(t,lambda);
%   figure; plot(t,cf),grid
%   title('CF of the symmetric Trapezoidal distribution on (-1,1)')
%
% EXAMPLE 2:
% % CF of a linear combination of independent Trapezoidal RVs
%   t = linspace(-20,20,201);
%   lambda = [3 3 4 4 5]/7;
%   coef = [1 2 3 4 5]/15;
%   cf = cf_TrapezoidalSymmetric(t,lambda,coef);
%   figure; plot(t,cf),grid
%   title('CF of a linear combination of independent Trapezoidal RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a weighted linear combination of independent Trapezoidal RVs
%   lambda = [3 3 4 4 5]/7;
%   coef = [1 2 3 4 5]/15;
%   cf   = @(t) cf_TrapezoidalSymmetric(t,lambda,coef);
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
%cf = cf_TrapezoidalSymmetric(t,theta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid  = []; end
if nargin < 3, coef  = []; end
if nargin < 2, lambda = []; end

if isempty(lambda), lambda = 1; end
if isempty(coef), coef  = 1; end
if isempty(niid), niid  = 1; end

[errorcode,coef,lambda] = distchck(2,coef(:)',lambda(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt = size(t);
t   = t(:);
w1  = coef.*(1/2+lambda/2);
w2  = coef.*(1/2-lambda/2);

cf  = prod((sin(t*w1)./(t*w1)) .* (sin(t*w2)./(t*w2)),2);
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