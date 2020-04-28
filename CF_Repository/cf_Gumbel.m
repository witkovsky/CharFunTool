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
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% SPECIAL CASES:
%   1) If X is an exponential RV with mean 1, then -log(X) has a standard
%      Gumbel distribution (with mu = 0 and beta = 1). 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Gumbel_distribution.
%  For more details see also Marques, Coelho, De Carvalho (2015).
%
% EXAMPLE 1:
% % CF of the the linear combination of Gumbel RVs
%   mu   = [2 3];
%   beta = [5 6];
%   coef = [1 1];
%   t    = linspace(-1,1,201);
%   cf   = cf_Gumbel(t,mu,beta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the Gumbel RVs')
%
% EXAMPLE 2 (Scenario I in Marques etal (2015)):
% % PDF/CDF of the linear combination of Gumbel RVs 
%   mu   = [2 3];
%   beta = [5 6];
%   coef = [1 1];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Gumbel(t,mu,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%
% EXAMPLE 3 (Scenario II in Marques etal (2015)):
% % PDF/CDF of the linear combination of Gumbel RVs 
%   mu   = [-4 -1 2 3];
%   beta = [0.1 0.2 0.3 0.4];
%   coef = [1 2 3 4];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Gumbel(t,mu,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%
% EXAMPLE 4 (Scenario III in Marques etal (2015)):
% % PDF/CDF of the linear combination of Gumbel RVs 
%   mu   = [-10 10 20 30 40];
%   beta = [1 2 3 4 5];
%   coef = [1/2 1 3/4 5 1];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Gumbel(t,mu,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%
% REFERENCES:
%  [1] Marques, F.J., Coelho, C.A., De Carvalho, M. (2015). On the
%      distribution of linear combinations of independent Gumbel random
%      variables. Statistics and Computing, 25(3), pp.683-701.   

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Aug-2018 16:57:55
% Rev.: 28-Apr-2020 13:47:42

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