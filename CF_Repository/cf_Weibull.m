function cf = cf_Weibull(t,alpha,beta,coef,niid,tol)
%% cf_Weibull 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent Weibull random variables with the parameters alpha (scale
%  parameters, alpha > 0) and beta (shape parameters, beta > 0), for real
%  (vector) argument t. 
%
%  That is, cf_Weibull evaluates the characteristic function cf(t) of
%  Y = sum_{i=1}^N coef_i * X_i, where X_i ~ Weibull(alpha_i,beta_i) are
%  inedependent Weibull RVs, with the parameters alpha_i > 0 and beta_i,
%  for i = 1,...,N. 
%
%  The characteristic function of Y is defined by
%   cf(t) = Prod(cf_i(coef(i)*t,alpha_i,beta_i)),
%  where cf_i(coef(i)*t,alpha_i,beta_i) is CF of the random variable
%  Y_i = coef(i)*X_i.
%
% SYNTAX:
%  cf = cf_Weibull(t,alpha,beta,,coef,niid,tol)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  alpha  - vector of the 'scale' parameters alpha > 0. If empty, default
%           value is alpha = 1.  
%  beta   - vector of the 'shape' parameters beta > 0. If empty, default
%           value is beta = 1.  
%  coef   - vector of the coefficients of the linear combination of the
%           WEIBULL random variables. If coef is scalar, it is assumed
%           that all coefficients are equal. If empty, default value is
%           coef = 1.
%  niid   - scalar convolution coeficient niid, such that Z = Y +...+ Y is
%           sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%           coef(i) * X_i is independently and identically distributed
%           random variable. If empty, default value is niid = 1.   
%  tol    - parameter of the relative tolerance RelTol used for
%           integration. If empty, default value is tol = 1e-6.  
% 
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Weibull_distribution.
%
% EXAMPLE 1:
%  alpha = 1;
%  beta  = 0.5;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cf_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('CF of the Weibull distribution')
%
% EXAMPLE 2: 
%  alpha = [1 1 2];
%  beta  = [0.5 1 3];
%  coef  = [1 1 1]/3;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cf_Weibull(t,alpha,beta,coef);
%  plot(t,real(cf),t,imag(cf));grid
%  title('CF of a linear combination of Weibull RVs')
%
% EXAMPLE 3: 
%  alpha = [1 1 2];
%  beta  = [0.5 1 3];
%  coef  = [1 1 1]/3;
%  x     = linspace(0,10,201);
%  prob  = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^10;
%  cf = @(t) cf_Weibull(t,alpha,beta,coef);
%  result = cf2DistGP(cf,x,prob,options)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Oct-2017 12:47:41
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf_Weibull(t,alpha,beta,,coef,niid,tol)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, tol   = []; end
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 2, beta  = []; end
if nargin < 2, alpha = []; end

%%
if isempty(alpha)
    alpha = 1;
end

if isempty(beta)
    beta = 1;
end

if isempty(coef)
    coef = 1;
end

if isempty(niid)
    niid = 1;
end

if isempty(tol)
    tol = 1e-6;
end

%% Check/set equal dimensions for the vectors coef, df, and ncp

[errorcode,coef,alpha,beta] = distchck(3,coef(:),alpha(:),beta(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function

szt = size(t);
n   = length(coef);
t   = t(:);
cf  = 1;
for i = 1:n
    cf = cf .* cfX_Weibull(coef(i)*t,alpha(i),beta(i),tol);
end
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