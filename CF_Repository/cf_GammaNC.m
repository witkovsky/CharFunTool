function cf = cf_GammaNC(t,alpha,beta,delta,coef,niid,tol)
%% cf_GammaNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Gamma random variables with distributions
%  Gamma(alpha_i,beta_i,delta_i), where the shape parameters alpha_i > 0
%  and the rate parameters beta_i > 0, and the non-centrality parameters
%  delta_i >= 0 for i = 1,...,N.
%
%  That is, cf_GammaNC evaluates the characteristic function cf(t)
%  of Y = coef_i*X_1 +...+ coef_N*X_N, where X_i ~
%  Gamma(alpha_i,beta_i,delta_i) are inedependent RVs, with parameters
%  alpha_i > 0, beta_i > 0, and delta_i >= 0, for i = 1,...,N.   
%
%  The characteristic function of X ~ Gamma(alpha,beta,delta) is
%  Poisson mixture of CFs of the central Gamma RVs of the form
%   cf(t) = cf_GammaNC(t,alpha,beta,delta) = 
%         = exp(-delta/2) sum_{j=1}^Inf (delta/2)^j/j! *
%           * cf_Gamma(alpha+j,beta) 
%  where cf_Gamma(alpha+j,beta) are the CFs of central Gamma RVs with
%  parameters alpha+j and beta. Equivalenly, we get
%   cf(t) = cf_GammaNC(t,alpha,beta,delta) = 
%         = exp( i*t*delta/(2*beta*(1-i*t/beta)) ) * (1-i*t/beta)^(-alpha).
%
% SYNTAX
%  cf = cf_GammaNC(t,alpha,beta,delta,coef,niid,tol)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'rate' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0.
%  coef  - vector of the coefficients of the linear combination of the
%          gamma distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Gamma_distribution
%
% EXAMPLE 1:
% % CF of the non-central Gamma RV with delta = 1
%   alpha = 3/2;
%   beta  = 1/2;
%   delta = 1;
%   t = linspace(-10,10,201);
%   cf = cf_GammaNC(t,alpha,beta,delta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the non-central Gamma RV with delta = 1')
%
% EXAMPLE 2:
% % CDF/PDF of the non-central Gamma RV
%   alpha = 3/2;
%   beta  = 1/2;
%   delta = 1;
%   cf = @(t) cf_GammaNC(t,alpha,beta,delta);
%   clear options;
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options)
%
% EXAMPLE 3:
% % CDF/PDF of the linear combination of the non-central Gamma RVs 
%   alpha = [3 4 5]/2;
%   beta  = [1 1 2]/2;
%   delta = [0 1 2];
%   coef  = [1 1 1]/3;
%   cf = @(t) cf_GammaNC(t,alpha,beta,delta,coef);
%   clear options;
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Oct-2022 23:50:16

%% CHECK THE INPUT PARAMETERS
narginchk(1, 7);

if nargin < 7, tol   = []; end
if nargin < 6, niid  = []; end
if nargin < 5, coef  = []; end
if nargin < 4, delta = []; end
if nargin < 3, beta  = []; end
if nargin < 2, alpha = []; end

if isempty(alpha), alpha = 1; end
if isempty(beta), beta = 1; end
if isempty(delta), delta = 0; end
if isempty(coef), coef = 1; end
if isempty(tol), tol = 1e-12; end

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,alpha,beta,delta] = ...
    distchck(4,coef(:)',alpha(:)',beta(:)',delta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent nc RVs
szt   = size(t);
t     = t(:);
cf    = 1;
aux   = 1 - bsxfun(@times,1i*t,coef./beta);
if any(delta~=0)
    aux = bsxfun(@power,aux,-alpha) .* ...
        exp(bsxfun(@times,1i*t,0.5*coef.*delta./beta)./aux);
else
    aux = bsxfun(@power,aux,-alpha);
end
cf   = cf .* prod(aux,2);
cf   = reshape(cf,szt);

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end