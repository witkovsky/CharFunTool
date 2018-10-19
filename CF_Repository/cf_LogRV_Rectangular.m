function cf = cf_LogRV_Rectangular(t,a,b,coef,niid)
%% cf_LogRV_Rectangular
%  Characteristic function of a linear combination (resp. convolution) of
%  independent  LOG-TRANSFORMED RECTANGULAR (continuous UNIFORM) random
%  variables defined on the interval (a,b), with b > a >= 0.  
%
%  That is, cf_LogRV_Rectangular evaluates the characteristic function
%  cf(t) of  Y = sum_{i=1}^N coef_i * log(X_i), where X_i ~ R(a_i,b_i) are
%  independent uniformly distributed RVs defined on (a_i,b_i), for all i =
%  1,...,N.  
%
%  The characteristic function of Y = log(X), with X ~ R(a,b) is
%  defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)). 
%  That is, the characteristic function can be derived from expression for
%  the r-th moment of X, E(X^r) by using (1i*t) instead of r. Here, the
%  r-th raw moments of the RECTANGULAR distribution on the interval (a,b),
%  with 0 <= a < b, are given analytically by (for more details see, e.g.,
%  http://mathworld.wolfram.com/UniformDistribution.html): 
%    mu_r = (b^(r+1) -  a^(r+1))/((r+1)*(b-a)).
%
%  Hence, the characteristic function of Y = log(X), with X ~
%  R(a,b) is  defined by   
%    cf_Y(t) = (b^(1i*t+1) -  a^(1i*t+1)) / ((1i*t+1)*(b-a)).
%  Then, the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters a(i) and b(i).
%
% SYNTAX
%   cf = cf_LogRV_Rectangular(t,a,b,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  a     - vector of the parameters a (lower bounds of the interval). Here,
%          it is assumed that 0 <= a. If empty, default value is a = 0. 
%  b     - vector of the parameters b (upper bounds of the interval). Here,
%          it is assumed that 0 <= a < b. If empty, default value is b = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          RECTANGULAR distributed random variables. If coef is scalar, it
%          is assumed that all coefficients are equal. If empty, default
%          value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
%
% WOLFRAM MATH WORLD: 
%  http://mathworld.wolfram.com/UniformDistribution.html
%
% EXAMPLE 1:
% % CF of the LOG-TRANSFORMED RECTANGULAR distribution on (0,1)
%   t  = linspace(-20,20,201);
%   cf = cf_LogRV_Rectangular(t);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the Rectangular distribution on (0,1)')
%
% EXAMPLE 2:
% % CF of a combination of LOG-TRANSFORMED RECTANGULAR RVs
%   t    = linspace(-5,5,201);
%   a    = [0 0 0 0 0];
%   b    = [1 2 3 4 5];
%   coef = [1 2 3 4 5]/15;
%   cf   = cf_LogRV_Rectangular(t,a,b,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of a combination of LOG-TRANSFORMED RECTANGULAR RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a combination of LOG-TRANSFORMED RECTANGULAR RVs
%   a    = [0 0 0 0 0];
%   b    = [1 2 3 4 5];
%   coef = [1 2 3 4 5]/15;
%   cf   = @(t) cf_LogRV_Rectangular(t,a,b,coef);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options)
%
% REFERENCES:
%   WITKOVSKY V. (2016). Numerical inversion of a characteristic
%   function: An alternative tool to form the probability distribution of
%   output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Oct-2018 11:26:23

%% ALGORITHM
%cf = cf_LogRV_Rectangular(t,a,b,coef,niid)

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

if any(a<0)
    error(message('Negative lower bounds'));
end
    
%% Characteristic function
szt = size(t);
t   = t(:);
cf  = prod((b.^(1i*t*coef +1) - a.^(1i*t*coef +1)) ./ ((1i*t+1)*(b-a)),2);
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