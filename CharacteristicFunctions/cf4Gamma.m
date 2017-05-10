function cf = cf4Gamma(t,k,theta,coef,n)
%%CF4GAMMA Characteristic function of the distribution of a convolution (or
%  a linear combination) of independent GAMMA random variables X_i ~ GAMMA.
%
%  In particular, cf4Gamma(t,k,theta,coef) evaluates the characteristic
%  function cf(t) of Y = coef(1) * X_1 + ... + coef(N) * X_N, where X_i ~
%  GAMMA(k(i),theta(i)), and k(i) and  theta(i) represent the 'shape' and
%  the 'scale' parameters of the GAMMA distribution.  
%
%  The characteristic function of Y is defined by
%   cf(t) = Prod( (1 - i*t*theta(i)*coef(i))^(-k(i)) )
%
%  Alternative parametrization of the GAMMA distribution is using the
%  'rate' parameter instead of the 'shape' parameter (rate is the
%  reciprocal value of the shape parameter). For more details see, e.g,
%  WIKIPEDIA: https://en.wikipedia.org/wiki/Gamma_distribution.
%
% SYNTAX
%  cf = cf4Gamma(t,k,theta,coef,n)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  k     - vector or array of the shape parameters k>0. If empty, default
%          value is k = 1. 
%  theta - vector of the non-centrality parameters of the the chi-squared
%          random variables. If ncp is scalar, it is assumed that all
%          non-centrality parameters are equal. If empty, default value is
%          ncp = 0. 
%  coef  - vector of the coefficients of the lienar combination of the
%          chi-squared random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  n     - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of n iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * GAMMA(k(i),theta(i))) is independently and identically
%          distributed random variable. If empty, default value is n = 1.   
%
% EXAMPLE 1:
% % CF of a linear combination of K=100 independent Gamma RVs)
%   K = 50;
%   t = linspace(-10,10,201);
%   idx = 1:K;
%   coef = 1./((idx - 0.5)*pi).^2;
%   figure; plot(idx,coef,'.-')
%   title('Coefficients of the linear combination of GAMMA RVs')
%   k = 5/2;
%   theta = 2;
%   cf = cf4Gamma(t,k,theta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the linear combination of GAMMA RVs')
%
% EXAMPLE 2:
% % PDF/CDF from the CF by cf2DistGP
%   k = 5/2;
%   theta = 2;
%   coef = 1./(((1:50) - 0.5)*pi).^2;
%   cf = @(t) cf4Gamma(t,k,theta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options);

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-May-2017 18:11:50

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, n = []; end
if nargin < 4, coef = []; end
if nargin < 3, theta = []; end
if nargin < 2, k = []; end

%%
if isempty(theta) && ~isempty(k)
    theta = 1;
elseif isempty(theta) && ~isempty(coef)
    theta = 1;
elseif ~any(theta)
    theta = 1;
end

if isempty(k) && ~isempty(coef)
    k = 1;
elseif isempty(k) && ~isempty(theta)
    k = 1;
end

if isempty(coef) && ~isempty(theta)
    coef = 1;
elseif isempty(coef) && ~isempty(k)
    coef = 1;
end

if isempty(n)
    n = 1;
end

%% Equal size of the parameters   
if ~isempty(coef) && isscalar(k) && isscalar(theta) && isempty(n)
    coef = sort(coef);
    m    = length(coef);
    [coef,idx] = unique(coef);
    k = k * diff([idx;m+1]);
end

[errorcode,coef,k,theta] = distchck(3,coef(:)',k(:)',theta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

% Special treatment for linear combinations with large number of RVs
szcoefs  = size(coef);
szcoefs  = szcoefs(1)*szcoefs(2);
szt      = size(t);
sz       = szt(1)*szt(2);
szcLimit = ceil(1e3 / (sz/2^16));
idc = 1:fix(szcoefs/szcLimit)+1;

%% ALGORITHM
t        = t(:);
idx0     = 1;
cf       = 1;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux  = bsxfun(@times,t,theta(idx).*coef(idx));
    aux  = bsxfun(@power,(1 - 1i * aux),-k(idx));
    cf   = cf .* prod(aux,2);
end
cf = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(n)
    if isscalar(n)
        cf = cf .^ n;
    else
        error('n should be a scalar (positive integer) value');
    end
end

end