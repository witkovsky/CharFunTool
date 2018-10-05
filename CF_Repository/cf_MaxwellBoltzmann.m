function cf = cf_MaxwellBoltzmann(t,scale,coef,niid)
%cf_MaxwellBoltzmann 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Maxwell-Boltzmann distributed random variables.     
%  
%  The Maxwell-Boltzmann distribution is a continuous probability
%  distribution for positive-valued random variables. It is a chi
%  distribution in three degrees of freedom.
%
%  cf_MaxwellBoltzmann evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ MaxwellBoltzmann(scale_i) are
%  inedependent Maxwell-Boltzmann distributed RVs with the scale parameters
%  scale_i > 0, for i  = 1,...,N.
%
%  The characteristic function of X ~ MaxwellBoltzmann(scale) is defined by
%   cf_MaxwellBoltzmann(t) = cf_Chi(scale*t,df=3), 
%  where cf_Chi(t,df) denotes the characteristic function of the Chi
%  distribution with df degrees of freedom. Hence, the characteristic
%  function of Y is 
%   cf(t) = Prod ( cf_MaxwellBoltzmann(t,scale_i) )
%
% SYNTAX:
%  cf = cf_MaxwellBoltzmann(t,scale,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  scale - vector of the scale parameters of the Maxwell-Boltzmann
%          distributed random variables.  If scale is scalar, it is assumed
%          that all scale parameters are equal. If empty, default value is
%          scale = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          Maxwell-Boltzmann distributed random variables. If coef is
%          scalar, it is assumed that all coefficients are equal. If empty,
%          default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
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
%  velocity vector in Euclidean space), with a scale parameter measuring
%  speeds in units proportional to the square root of T/m (the ratio of
%  temperature and particle mass).       
%
% EXAMPLE 1:
% % CF of the distribution of Maxwell-Boltzmann RV with scale = 3
%   scale = 3;
%   t     = linspace(-3,3,501);
%   cf    = cf_MaxwellBoltzmann(t,scale);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the Maxwell-Boltzmann RV with scale = 3')
%
% EXAMPLE 2: 
% % CF of a linear combination of independent Maxwell-Boltzmann RVs
%   scale = [1 2 3];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_MaxwellBoltzmann(t,scale,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Maxwell-Boltzmann RVs')
%
% EXAMPLE 3:
% % PDF/CDF of a linear combination of independent Maxwell-Boltzmann RVs
%   scale = [1 2 3];
%   coef  = [1 1 1];
%   cf    = @(t) cf_MaxwellBoltzmann(t,scale,coef);
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

%% ALGORITHM
%  cf = cf_MaxwellBoltzmann(t,scale,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid  = []; end
if nargin < 3, coef  = []; end
if nargin < 2, scale = []; end

if isempty(scale)
    scale = 1;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and scale 
[errorcode,coef,scale] = distchck(2,coef(:),scale(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

% CF of the linear combination of the Maxwell-Boltzmann RVs (by using Chi)
df = 3;
cf = cf_Chi(t,df,scale.*coef,niid);

end