function cf = cf_Gumbel(t,mu,beta,coef,niid)
%% cf_Gumbel 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent GUMBEL (Generalized Extreme Value distribution Type-I)
%  random variables with location mu (real) and scale beta > 0.
%
%  That is, cf_Gumbel evaluates the characteristic function cf(t)  of  Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Gumbel(mu_i,beta_i) are
%  inedependent RVs, with real location parameters mu_i and the scale
%  parameters beta_i > 0, for i = 1,...,N. 
%
%  The characteristic function of Y is defined by
%   cf(t) = Prod( gamma(1 - 1i*t*coef(i)*beta(i)) * exp(1i*t*coef(i)*mu(i)) )
%
% SYNTAX:
%  cf = cf_Gumbel(t,mu,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of real location parameter. If empty, default value is 
%          mu = 0.   
%  beta  - vector of the scale parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          GUMBEL random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% SPECIAL CASES:
%   1) If X is an exponential RV with mean 1, then -log(X) has a standard
%      Gumbel distribution (with mu = 0 and beta = 1). 
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Gumbel_distribution.
%
% EXAMPLE 1:
% % CF of the Gumbel RV
%   mu   = 1;
%   beta = 2;
%   t    = linspace(-10,10,201);
%   cf   = cf_Gumbel(t,mu,beta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the Gumbel RVs')
%
% EXAMPLE 2:
% % PDF/CDF of the linear combination of Gumbel RV from its CF by cf2DistGP
%   mu   = [0 2 -1];
%   beta = [1 2 1];
%   coef = [1 2 3];
%   cf   = @(t) cf_Gumbel(t,mu,beta,coef);
%   clear options
%   options.N = 2^10;
%   result = cf2DistGP(cf,[],[],options);

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Aug-2018 16:57:55

%% ALGORITHM
% cf = cf_Gumbel(t,mu,beta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, mu = []; end

%%
if isempty(mu), mu = 0; end
if isempty(beta), beta = 1; end
if isempty(coef), coef = 1; end


%% Equal size of the parameters   

[errorcode,coef,mu,beta] = distchck(3,coef(:)',mu(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt = size(t);
t   = t(:);
cf  = prod(GammaZX(1-bsxfun(@times,1i*t,coef.*beta)) .* ...
    exp(bsxfun(@times,1i*t,coef.*mu)),2);    
cf = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end