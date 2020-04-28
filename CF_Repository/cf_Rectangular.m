function cf = cf_Rectangular(t,a,b,coef,niid)
%% cf_Rectangular
%  Characteristic function of a linear combination (resp. convolution) of
%  independent RECTANGULAR (continuous UNIFORM) random variables defined on
%  the interval (a,b).  
%
%  That is, cf_Rectangular evaluates the characteristic function
%  cf(t) of  Y = sum_{i=1}^N coef_i * X_i, where X_i ~ R(a_i,b_i) are
%  independent uniformly distributed RVs defined on (a_i,b_i), for all i =
%  1,...,N.  
%
%  The characteristic function of RECTANGULAR RV X ~ R(a,b) is defined by
%   cf(t) = cf_Rectangular(t,a,b) = (exp(1i*t*b) -exp(1i*t*a))/(1i*t*(b-a))
%
% SYNTAX
%   cf = cf_Rectangular(t,a,b,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  a     - vector of the parameters a (lower bounds of the interval). If
%          empty, default value is a = 0.  
%  b     - vector of the parameters b (upper bounds of the interval). Here,
%          it is assumed that a < b. If empty, default value is b = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          RECTANGULAR distributed random variables. If coef is scalar, it
%          is assumed that all coefficients are equal. If empty, default
%          value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
%
% WOLFRAM MATH WORLD: 
%  http://mathworld.wolfram.com/UniformDistribution.html
%
% EXAMPLE 1:
% % CF of the Rectangular distribution on (0,1)
%   t  = linspace(-50,50,201);
%   cf = cf_Rectangular(t);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the Rectangular distribution on (0,1)')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of independent Rectangular RVs
%   t    = linspace(-10,10,201);
%   a    = [0 0 0 0 0];
%   b    = [1 2 3 4 5];
%   coef = [1 2 3 4 5]/15;
%   cf   = cf_Rectangular(t,a,b,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of a weighted linear combination of Rectangular RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a weighted linear combination of independent Rectangular RVs
%   a    = [0 0 0 0 0];
%   b    = [1 2 3 4 5];
%   coef = [1 2 3 4 5]/15;
%   cf   = @(t) cf_Rectangular(t,a,b,coef);
%   x    = linspace(0,5,201);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.N = 2^12;
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
%   WITKOVSKY V. (2016). Numerical inversion of a characteristic
%   function: An alternative tool to form the probability distribution of
%   output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Oct-2018 11:26:23
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%cf = cf_Rectangular(t,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, b    = []; end
if nargin < 2, a    = []; end

if isempty(a), a = 0; end
if isempty(b), b = 1; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Check size of the parameters
[errorcode,coef,a,b] = distchck(3,coef(:)',a(:)',b(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt = size(t);
t   = t(:);
cf  = prod((exp(1i*t*(b.*coef))-exp(1i*t*(a.*coef))) ./ ... 
     (1i*t*(coef.*(b-a))),2);
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