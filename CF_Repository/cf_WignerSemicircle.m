function cf = cf_WignerSemicircle(t,mu,R,coef,niid)
%% cf_WignerSemicircle
%  Characteristic function of a linear combination (resp. convolution) of
%  independent WIGNER SEMICIRCLE random variables  defined on the interval
%  (mu-R,mu+R). 
%
%  That is, cf_WignerSemicircle evaluates the characteristic function
%  cf(t) of  Y = sum_{i=1}^N coef_i * X_i, where X_i ~ WignerSemicircle
%  are independent RVs defined on (mu_i-R_i,mu_i+R_i), for all i = 1,...,N.
%
%  The characteristic function of X ~ WignerSemicircle(mu,R) is defined by
%   cf(t) = 2*exp(1i*t*mu).*besselj(1,R*t)./(R*t);
%
% SYNTAX
%  cf = cf_WignerSemicircle(t,mu,R,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the 'location' parameters mu in R. If empty, default
%          value is mu = 0.  
%  R     - vector of the 'radius' parameters R > 0. If empty, default
%          value is sigma = 1. 
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
% % CF of the symmetric Wigner Semicircle distribution on (-1,1)
%   t = linspace(-30,30,501);
%   cf = cf_WignerSemicircle(t);
%   figure; plot(t,cf),grid
%   title('CF of the the Wigner Semicircle distribution on (-1,1)')
%
% EXAMPLE 2:
% % PDF/CDF of a Wigner Semicircle RVs
%   cf   = @(t) cf_WignerSemicircle(t);
%   x    = linspace(-1,1,201);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = -1;
%   options.xMax = 1;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3:
% % CF of a linear combination of independent Wigner Semicircle RVs
%   t = linspace(-1,1,501);
%   mu   = [0 1 2 1 0];
%   R    = [1 1 2 2 3];
%   coef = [1 2 3 4 5];
%   cf = cf_WignerSemicircle(t,mu,R,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of a linear combination of independent Wigner Semicircle RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent Wigner Semicircle RVs
%   mu   = [0 1 2 1 0];
%   R    = [1 1 2 2 3];
%   coef = [1 2 3 4 5];
%   cf   = @(t) cf_WignerSemicircle(t,mu,R,coef);
%   x    = linspace(-20,40,201);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
%   WITKOVSKY V. (2016). Numerical inversion of a characteristic
%   function: An alternative tool to form the probability distribution of
%   output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Sep-2018 00:32:58
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%cf = cf_WignerSemicircle(t,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, R    = []; end
if nargin < 4, mu    = []; end

if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end
if isempty(mu),     mu = 0; end
if isempty(R),       R = 1; end

%% Characteristic function
[errorcode,coef,mu,R] = distchck(3,coef(:)',mu(:)',R(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

szt  = size(t);
t    = t(:);
cf   = prod(2*exp(1i*t*(mu.*coef)) .* ...
       besselj(1,t*(R.*coef)) ./ (t*(R.*coef)),2);
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