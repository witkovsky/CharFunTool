function cf = cf_Rayleigh(t,sigma,coef,niid)
%cf_Rayleigh 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Rayleigh random variables.     
%  
%  The Rayleigh distribution is a continuous probability distribution for
%  positive-valued random variables. It is a chi distribution in two
%  degrees of freedom.
%
%  cf_Rayleigh evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Rayleigh(sigma_i) are inedependent
%  Rayleigh RVs with the scale parameters sigma_i > 0, for i  = 1,...,N.
%
%  The characteristic function of X ~ Rayleigh(sigma) is defined by
%   cf_Rayleigh(t) = cf_Chi(sigma*t,df=2), 
%  where cf_Chi(t,df) denotes the characteristic function of the Chi
%  distribution with df degrees of freedom. Hence, the characteristic
%  function of Y is 
%   cf(t) = Prod ( cf_Rayleigh(coef_i*t,sigma_i) )
%
% SYNTAX:
%  cf = cf_Rayleigh(t,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Rayleigh random
%          variables.  If sigma is scalar, it is assumed that all scale
%          parameters are equal. If empty, default value is sigma = 1. 
%  coef  - vector of the coefficients of the linear combination of the
%          Rayleigh random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Rayleigh_distribution
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% NOTES (from Wikipedia)
%  A Rayleigh distribution is often observed when the overall magnitude of
%  a vector is related to its directional components. One example where the
%  Rayleigh distribution naturally arises is when wind velocity is analyzed
%  in two dimensions. Assuming that each component is uncorrelated,
%  normally distributed with equal variance, and zero mean, then the
%  overall wind speed (vector magnitude) will be characterized by a
%  Rayleigh distribution. A second example of the distribution arises in
%  the case of random complex numbers whose real and imaginary components
%  are independently and identically distributed Gaussian with equal
%  variance and zero mean. In that case, the absolute value of the complex
%  number is Rayleigh-distributed.
%
% EXAMPLE 1:
% % CF of the distribution of Rayleigh RV with sigma = 3
%   sigma = 3;
%   t     = linspace(-5,5,501);
%   cf    = cf_Rayleigh(t,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the Rayleigh RV with sigma = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent Rayleigh RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   t     = linspace(-1,1,501);
%   cf    = cf_Rayleigh(t,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Rayleigh RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Rayleigh RVs
%   sigma = [1 2 3 4 5];
%   coef  = [1 1 1 1 1];
%   cf    = @(t) cf_Rayleigh(t,sigma,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(5,40,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_Chi, cf_MaxwellBoltzmann

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Oct-2018 13:47:29
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_Rayleigh(t,sigma,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid  = []; end
if nargin < 3, coef  = []; end
if nargin < 2, sigma = []; end

if isempty(sigma)
    sigma = 1;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and sigma 
[errorcode,coef,sigma] = distchck(2,coef(:),sigma(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

% CF of the linear combination of the Rayleigh RVs (expressed by using Chi)
df = 2;
cf = cf_Chi(t,df,sigma.*coef,niid);

end