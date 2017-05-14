function cf = cf_Normal(t,mu,sigma,coef,n)
%%cf_Normal Characteristic function of a linear combination (resp.
%  convolution) of independent Normal random variables, X ~
%  Normal(mu,sigma).  
%
%  In particular, cf_Normal(t,alpha,beta,coef) evaluates the characteristic
%  function cf(t) of Y = coef(1) * X_1 + ... + coef(N) * X_N, where X_i ~
%  Normal(mu(i),sigma(i)), and mu(i) is real location parameter and
%  sigma(i)>0 represent the 'scale' parameter of the Normal distribution.
%
%  The characteristic function of Y is defined by
%   cf(t) = exp(1i*(coef'*mu)*t - (1/2)*(coef^2'*sigma^2)*t^2 )
%
% SYNTAX:
%  cf = cf_Normal(t,alpha,beta,coef,n)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the 'location' parameters mu in R. If empty, default
%          value is mu = 0.  
%  sigma - vector of the 'scale' parameters sigma > 0. If empty, default
%          value is sigma = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          Normal random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  n     - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of n iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i, with X_i ~ Normal(mu(i),sigma(i))) independently
%          and identically distributed random variables. If empty, default
%          value is n = 1.    
%
% REMARK:
%   The characteristic function of a lienar combination of independent
%   Normal random variables has well known exact form.  
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Normal_distribution.
%
% EXAMPLE 1:
% % CF of a linear combination of K=100 independent Norma RVs)
%   coef = 1./(((1:50) - 0.5)*pi).^2;
%   figure; plot(coef,'.-');grid on
%   title('Coefficients of the linear combination of Normal RVs')
%   mu = linspace(-3,3,50);
%   sigma = linspace(0.1,1.5,50);
%   t = linspace(-100,100,2001);
%   cf = cf_Normal(t,mu,sigma,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the linear combination of Normal RVs')
%
% EXAMPLE 2:
% % PDF/CDF from the CF by cf2DistGP
%   mu    = linspace(-3,3,50);
%   sigma = linspace(0.1,1.5,50);
%   coef  = 1./(((1:50) - 0.5)*pi).^2;
%   cf    = @(t) cf_Normal(t,mu,sigma,coef);
%   clear options
%   options.N = 2^10;
%   options.SixSigmaRule = 8;
%   prob = [0.01 0.05 0.1 0.5 0.9 0.950 .99];
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-May-2017 18:11:50

%% ALGORITHM
% cf = cf_Normal(t,mu,sigma,coef,n)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, n = []; end
if nargin < 4, coef = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

%%
if isempty(sigma) && ~isempty(mu)
    sigma = 1;
elseif isempty(sigma) && ~isempty(coef)
    sigma = 1;
elseif ~any(sigma)
    sigma = 1;
end

if isempty(mu) && ~isempty(coef)
    mu = 0;
elseif isempty(mu) && ~isempty(sigma)
    mu = 0;
end

if isempty(coef) && ~isempty(sigma)
    coef = 1;
elseif isempty(coef) && ~isempty(mu)
    coef = 1;
end

if isempty(n)
    n = 1;
end

%% Equal size of the parameters   
[errorcode,coef,mu,sigma] = distchck(3,coef(:)',mu(:)',sigma(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end
szt      = size(t);
t        = t(:);

mean = coef*mu';
var  = sum((coef.*sigma).^2);
cf   = exp(1i*mean*t - (var/2)*t.^2);
cf   = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(n)
    if isscalar(n)
        cf = cf .^ n;
    else
        error('n should be a scalar (positive integer) value');
    end
end

end