function cf = cf_MaxwellBoltzmann(t,sigma,coef,niid)
%cf_MaxwellBoltzmann 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Maxwell-Boltzmann distributed random variables.     
%  
%  The Maxwell-Boltzmann distribution is a continuous probability
%  distribution for positive-valued random variables. It is a scaledd chi
%  distribution with three degrees of freedom.
%
%  cf_MaxwellBoltzmann evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ MaxwellBoltzmann(sigma_i) are
%  inedependent Maxwell-Boltzmann distributed RVs with the scale parameters
%  sigma_i > 0, for i  = 1,...,N.
%
%  The characteristic function of X ~ MaxwellBoltzmann(sigma) is defined by
%   cf_MaxwellBoltzmann(t) = cf_Chi(sigma*t,df=3), 
%  where cf_Chi(t,df) denotes the characteristic function of the Chi
%  distribution with df degrees of freedom. Hence, the characteristic
%  function of Y is 
%   cf(t) = Prod ( cf_MaxwellBoltzmann(coef_i*t,sigma_i) )
%
% SYNTAX:
%  cf = cf_MaxwellBoltzmann(t,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  sigma - vector of the scale parameters of the Maxwell-Boltzmann
%          distributed random variables.  If sigma is scalar, it is assumed
%          that all scale parameters are equal. If empty, default value is
%          sigma = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          Maxwell-Boltzmann distributed random variables. If coef is
%          scalar, it is assumed that all coefficients are equal. If empty,
%          default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% NOTES (from Wikipedia)
%  A Maxwell-Boltzmann distribution was first defined and used for describing
%  particle speeds in idealized gases, where the particles move freely
%  inside a stationary container without interacting with one another,
%  except for very brief collisions in which they exchange energy and
%  momentum with each other or with their thermal environment.
%  Mathematically, the Maxwell–Boltzmann distribution is the chi
%  distribution with three degrees of freedom (the components of the
%  velocity vector in Euclidean space), with a sigma parameter measuring
%  speeds in units proportional to the square root of T/m (the ratio of
%  temperature and particle mass).       
%
% EXAMPLE 1:
% % CF of the distribution of Maxwell-Boltzmann RV with sigma = 3
%   sigma = 3;
%   t     = linspace(-3,3,501);
%   cf    = cf_MaxwellBoltzmann(t,sigma);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the Maxwell-Boltzmann RV with sigma = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent Maxwell-Boltzmann RVs
%   sigma = [1 2 3];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_MaxwellBoltzmann(t,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Maxwell-Boltzmann RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Maxwell-Boltzmann RVs
%   sigma = [1 2 3];
%   coef  = [1 1 1];
%   cf    = @(t) cf_MaxwellBoltzmann(t,sigma,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,20,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_Chi, cf_Rayleigh 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Oct-2018 13:47:29
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_MaxwellBoltzmann(t,sigma,coef,niid)

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

% CF of the linear combination of the Maxwell-Boltzmann RVs (by using Chi)
df = 3;
cf = cf_Chi(t,df,sigma.*coef,niid);

end