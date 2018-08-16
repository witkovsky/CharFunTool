function cf = cf_Laplace(t,mu,beta,coef,niid)
%% cf_Laplace  
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LAPLACE random variables with location parameter mu (real)
%  and scale parameter beta > 0. 
%
%  That is, cf_Laplace evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Laplace (mu_i,beta_i) are
%  inedependent RVs, with real location parameters mu_i and the scale
%  parameters beta_i > 0, for i = 1,...,N. 
%
%  The characteristic function of Y is defined by
%   cf(t) = Prod( exp(1i*t*coef(i)*mu(i)) ./ (1 + (t*coef(i)*beta(i)).^2 )
%
% SYNTAX:
%  cf = cf_Laplace(t,mu,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of real location parameters. If empty, default value is 
%          mu = 0.   
%  beta  - vector of the scale parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          LAPLACE random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Laplace_distribution.
%
% EXAMPLE 1:
% % CF of the Laplace RV
%   mu   = 0;
%   beta = 1;
%   t    = linspace(-10,10,201);
%   cf   = cf_Laplace(t,mu,beta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the Laplace RVs')
%
% EXAMPLE 2:
% % PDF/CDF of the Laplace RV
%   mu   = 0;
%   beta = 1;
%   x    = linspace(-5,5,101);
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Laplace(t,mu,beta);
%   result = cf2DistGP(cf,x,prob);
%
% EXAMPLE 3:
% % PDF/CDF of the linear combination of Laplace RVs 
%   mu   = [-4 -1 2 3];
%   beta = [0.1 0.2 0.3 0.4];
%   coef = [1 2 3 4];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Laplace(t,mu,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%
% EXAMPLE 4:
% % PDF/CDF of the linear combination of Laplace  RVs 
%   mu   = [-10 10 20 30 40];
%   beta = [1 2 3 4 5];
%   coef = [1/2 1 3/4 5 1];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Laplace(t,mu,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Aug-2018 16:00:43

%% ALGORITHM
% cf = cf_Laplace(t,mu,beta,coef,niid)

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
cf  = prod(exp(bsxfun(@times,1i*t,coef.*mu)) ./ ...
    (1+bsxfun(@times,t.^2,(coef.*beta).^2) ),2);    
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