function cf = cf_Student(t,df,mu,sigma,coef,niid)
%% cf_Student 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent (location and scale shifted) STUDENT's t random variables.
%
%  That is, cf_Student evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * (mu_i + sigma_i * X_i), where X_i ~ t(df_i) are
%  inedependent (symmetric) t-distributed RVs, with df_i > 0 degrees of
%  freedom, for i = 1,...,N.
%
%  The characteristic function of the random variable mu + sigma*X, where
%  X ~ t(df) is given by 
%   cf(t) = exp(1i*t*mu) * besselk(df/2,abs(sigma*t)*sqrt(df),1) * ...
%           exp(-abs(sigma*t)*sqrt(df)) * (sqrt(df)*abs(aigma*t))^(df/2)...
%           / 2^(df/2-1)/gamma(df/2).
%
%  Hence, the characteristic function of Y  = coef_1*(mu_1+sigma_1*X_1)
%  +...+ coef_N*(mu_N+sigma_N*X_N) is cf_Y(t) = exp(1i*mu*t) *
%  (cf_1(coef_1*sigma_1*t) *...* cf_N(coef_N*sigma_N*t)), where cf_i(t) is
%  the characteristic function of X_i ~ t(df_i).
%
% SYNTAX:
%  cf = cf_Student(t,df,mu,sigma,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the degrees of freedom of the t-distributed random
%          variables.  If df is scalar, it is assumed that all degrees of
%          freedom are equal. If empty, default value is df = 1.
%  mu    - vector of location parameters, mu in Real. If empty, default
%          value is mu = 0. 
%  sigma - vector of scale parameters, sigma_i > 0. If empty, default value
%          is sigma = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef_i * X_i is independently and identically distributed random
%          variable. If empty, default value is niid = 1. 
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Student%27s_t-distribution.
%
% EXAMPLE 1:
%  % CF of a linear combination of independent Student's t RVs
%  coef = 1./(1:50);
%  df   = 50:-1:1;
%  t    = linspace(-1,1,201);
%  cf   = cf_Student(t,df,coef);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the linear combination of t RVs')
%
% EXAMPLE 2:
%  % CDF/PDF of a linear combination of independent Student's t RVs
%  coef = 1./(1:50);
%  df   = 50:-1:1;
%  cf   = @(t) cf_Student(t,df,coef);
%  x    = linspace(-50,50);
%  prob = [0.9 0.95 0.975 0.99];
%  clear options;
%  options.N = 2^12;%   
%  result = cf2DistGP(cf,x,prob,options);
%  disp(result)
%
% REFERENCES:
%   WITKOVSKY, V.: On the exact computation of the density and of the
%   quantiles of linear combinations of t and F random variables. Journal
%   of Statistical Planning and Inference 94 (2001), 1–13.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf_Student(t,df,mu,sigma,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6,  niid = []; end
if nargin < 5,  coef = []; end
if nargin < 4, sigma = []; end
if nargin < 3,    mu = []; end
if nargin < 2,    df = []; end

%%
if isempty(df) && ~isempty(coef)
    df = 1;
end

if isempty(coef) && ~isempty(df)
    coef = 1;
end

if isempty(mu), mu = 0; end
if isempty(sigma), sigma = 1; end
if isempty(niid), niid = 1; end

%% Equal size of the parameters
[errorcode,coef,df,mu,sigma] = ...
    distchck(4,coef(:)',df(:)',mu(:)',sigma(:)');
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

%% Characteristic function 
df2   = df/2;
t     = t(:);
o     = ones(length(t),1);
idx0  = 1;

for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux  = bsxfun(@times,abs(t),abs(sqrt(df(idx)).*coef(idx).*sigma(idx)));
    aux  = - aux + bsxfun(@times,df2(idx),log(aux)) + ...
        log(besselk(o*df2(idx),aux,1));
    aux = bsxfun(@plus,aux,(-log(2)*(df2(idx)-1))-gammaln(df2(idx)));
    aux = 1i*t*(coef(idx).*mu(idx)) + aux;
    cf  = exp(sum(aux,2));
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
