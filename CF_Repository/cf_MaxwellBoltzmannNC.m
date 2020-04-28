function cf = cf_MaxwellBoltzmannNC(t,scale,delta,coef,niid,tol)
%cf_MaxwellBoltzmannNC 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Maxwell-Boltzmann distributed random variables,
%  with the scale parameters scale_i > 0, and the non-centrality parameters
%  delta_i >= 0, for i  = 1,...,N.      
%  
%  The non-central Maxwell-Boltzmann distribution is a continuous
%  probability distribution for positive-valued random variables. It is a
%  special case of the non-central chi distribution in three degrees of
%  freedom.  
%
%  cf_MaxwellBoltzmannNC evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i~ MaxwellBoltzmannNC(scale_i,delta_i)
%  are inedependent non-central Maxwell-Boltzmann distributed RVs with the
%  scale parameters scale_i > 0,  and the non-centrality parameters
%  delta_i >= 0, for i  = 1,...,N.      
%
%  The characteristic function of X ~ MaxwellBoltzmannNC(scale,delta) is
%  defined by 
%   cf_MaxwellBoltzmannNC(t) = cf_ChiNC(scale*t,df=3,delta/scale), 
%  where by cf_ChiNC(t,df,delta) we denote the characteristic function of
%  the noncentral chi distribution with df degrees of freedom and the
%  non-centrality parameter delta. Hence, the characteristic function of Y
%  is  
%   cf(t) = Prod ( cf_MaxwellBoltzmannNC(t,scale_i,delta_i) )
%
% SYNTAX:
%  cf = cf_MaxwellBoltzmannNC(t,scale,delta,coef,niid,tol)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  scale - vector of the scale parameters of the Maxwell-Boltzmann
%          distributed random variables.  If scale is scalar, it is assumed
%          that all scale parameters are equal. If empty, default value is
%          scale = 1.
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0. Notice that the noncentrality
%          parameter delta can be interpreted as a square root of the sum
%          of standardized squared means, delta = sqrt(sum_{i=1}^3
%          mu_i^2/sigma^2_i). 
%  coef  - vector of the coefficients of the linear combination of the
%          Maxwell-Boltzmann distributed random variables. If coef is
%          scalar, it is assumed that all coefficients are equal. If empty,
%          default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%  https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
%
% NOTES (from Wikipedia)
%  Mathematically, the Maxwell–Boltzmann distribution with scale = 1 is the
%  chi distribution with three degrees of freedom (the components of the
%  velocity vector in Euclidean space).     
%
% EXAMPLE 1:
% % CF of the distribution of Maxwell-Boltzmann RV with scale = 3
%   scale = 1;
%   delta = 5;
%   t     = linspace(-5,5,501);
%   cf    = cf_MaxwellBoltzmannNC(t,scale,delta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the non-central Maxwell-Boltzmann RV')
%
% EXAMPLE 2:
% % PDF/CDF of the distribution of Maxwell-Boltzmann RV with scale = 3
%   scale = 1;
%   delta = 5;
%   cf    = @(t) cf_MaxwellBoltzmannNC(t,scale,delta);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a linear combination of independent Maxwell-Boltzmann RVs
%   scale = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_MaxwellBoltzmannNC(t,scale,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Maxwell-Boltzmann RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent Maxwell-Boltzmann RVs
%   scale = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   cf    = @(t) cf_MaxwellBoltzmannNC(t,scale,delta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,25,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% SEE ALSO: cf_ChiNC, cf_MaxwellBoltzmann 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-Oct-2018 23:48:44
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
%  cf = cf_MaxwellBoltzmannNC(t,scale,delta,coef,niid,tol)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, tol   = []; end
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 3, delta = []; end
if nargin < 2, scale = []; end

if isempty(scale)
    scale = 1;
end

if isempty(delta)
    delta = 0;
end

if isempty(coef) 
    coef = 1;
end

% Check/set equal dimensions for the vectors coef and scale 
[errorcode,coef,scale,delta] = distchck(3,coef(:),scale(:),delta(:));
if errorcode > 0
        error(message('InputSizeMismatch'));
end

%% CF OF THE NON-CENTRAL Maxwell-Boltzmann RVs (expressed by cf_ChiNC) 
%  Here, we assume delta = sqrt(sum(mu_i^2/sigma_i^2))
%  Alternatively, if delta = sqrt(sum(mu_i^2)), set delta = delta./sigma;

df = 3;
cf = cf_ChiNC(t,df,delta,scale.*coef,niid,tol);

end