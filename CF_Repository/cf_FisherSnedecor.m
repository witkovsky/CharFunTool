function cf = cf_FisherSnedecor(t,df1,df2,coef,niid,tol)
%% cf_FisherSnedecor 
%  Characteristic function of the distribution of a linear combination of
%  independent random variables with the central FISHER-SNEDECOR
%  F-distribution.   
%
%  That is, cf_FisherSnedecor evaluates the characteristic function cf(t)
%  of  Y = sum_{i=1}^N coef_i * X_i, where X_i ~ F(df1_i,df2_i) are
%  inedependent RVs, with df1_i and df2_i degrees of freedom, for i =
%  1,...,N.
%
%  The characteristic function of X ~ F(df1,df2), the F-distribution with
%  df1 and df2 degrees of freedom, evaluated at the real t from
%  (-inf,+inf), is defined by
%   CF(t) = (gamma(df1/2+df2/2)/gamma(df2/2)) * ...
%           HypergeometricU(df1/2, 1-df2/2, -1i*(df2/df1)*t),
%  where HypergeometricU(a,b,z) denotes the confluent hypergeometric
%  function of the second kind U(a,b,z). 
%
%  Hence, the characteristic function of Y  = coef_1*X_1 +...+ coef_N*X_N
%  is cf_Y(t) =  cf_1(coef_1*t) *...* cf_N(coef_N*t), where cf_i(t) 
%  is the characteristic function of X_i ~ F(df1_i,df2_i).
%
% SYNTAX
%  cf = cf_FisherSnedecor(t,df1,df2,coef,niid,tol)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df1   - vector of the  degrees of freedom df1 > 0. If empty, default
%          value is df1 = 1.  
%  df2   - vector of the  degrees of freedom df2 > 0. If empty, default
%          value is df2 = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.
%  tol   - relative tolerance used for integration.  If empty, default
%          value is tol = 1e-6.  
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/F-distribution.
%
% EXAMPLE 1: 
%   % CF of a linear combination of independent F RVs
%   df1   = 3:12;
%   df2   = 14:-1:5;
%   coef  = 1/10;
%   t     = linspace(-5,5,501);
%   cf    = cf_FisherSnedecor(t,df1,df2,coef);
%   figure; plot(t,real(cf),t,imag(cf))
%   title('Characteristic function of the linear combination of F RVs')
%
% EXAMPLE 2: 
%   % PDF/CDF  of a linear combination of independent F RVs
%   df1   = 3:12;
%   df2   = 14:-1:5;
%   coef  = 1/10;
%   cf    = @(t) cf_FisherSnedecor(t,df1,df2,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%
% REFERENCES:
% [1] PHILLIPS, P.C.B. The true characteristic function of the F
%     distribution. Biometrika (1982), 261-264. 
% [2] WITKOVSKY, V.: On the exact computation of the density and of the
%     quantiles of linear combinations of t and F random variables. Journal
%     of Statistical Planning and Inference 94 (2001), 1–13.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Sep-2019 21:57:11
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf_FisherSnedecor(t,df1,df2,coef,niid,tol)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, tol = []; end
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, df2 = []; end
if nargin < 2, df1 = []; end

if isempty(df2) && ~isempty(df1)
    df2 = 1;
elseif isempty(df2) && ~isempty(coef)
    df2 = 1;
elseif ~any(df2)
    df2 = 1;
end

if isempty(df1) && ~isempty(coef)
    df1 = 1;
elseif isempty(df1) && ~isempty(df2)
    df1 = 1;
end

if isempty(coef) && ~isempty(df2)
    coef = 1;
elseif isempty(coef) && ~isempty(df1)
    coef = 1;
end

if isempty(tol), tol = 1e-6; end
reltol = tol;

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,df1,df2] = distchck(3,coef(:)',df1(:)',df2(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent F RVs
szt = size(t);
t   = t(:);
cf  = ones(length(t),1);
for i = 1:length(coef)
    cf = cf .* integral(@(x) bsxfun(@times,funCF(df1(i),df2(i), ...
        coef(i)*t,(x/(1-x))^2),2*x/(1-x)^3),0,1,'ArrayValued',true, ...
        'RelTol',reltol);
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

%% Function funCF
function f = funCF(df1,df2,t,x)
%FUNCF Integrand function of the integral representation of the
%  characteristic function CF of the F-distribution with df1 and df2
%  degrees of freedom at the real argument t.
%
%  The characteristic function of F-distribution with df1 and df2 degrees 
%  of freedom, evaluated at the real t from (-inf,+inf), is defined by
%   CF(t) = (gamma(df1/2+df2/2)/gamma(df2/2)) * ...
%           HypergeometricU(df1/2, 1-df2/2, -1i*(df2/df1)*t),
%  where HypergeometricU(a,b,z) denotes the confluent hypergeometric
%  function of the second kind: U(a,b,z). 
%  
%  Here we use an integral representation of the hypergeometric function 
%  U(a,b,z), defined for purly complex argument z as (VW2016):
%   U(a,b,z) = gamma(1-b)/gamma(a-b+1) * Integral_0^inf cfFun(a,b,x) dx. 
%
% SYNTAX:
%   f = funCF(df1,df2,t,x)
%
% EXAMPLE
%  df1 = 5;
%  df2 = 4;
%  t   = 1:5;
%  x   = linspace(0,1);
%  f   = funCF(df1,df2,t,x);
%  plot(x,real(f),x,imag(f))

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Apr-2016 13:40:09

%% ALGORITHM FUNCF
t  = t(:);
nt = length(t);
o1 = ones(nt,1);
x  = x(:)';
nx = length(x);
o2 = ones(1,nx);
a  = df1/2;
b  = 1 - df2/2;
c  = (gammaln(a-b+1) - gammaln(1-b) - gammaln(a));
z  = -(df2/df1)*t;
f  = zeros(nt,nx);

% z == 0
id = (z==0);
if any(id)
    f(id,:) = exp(c + (a-1).*log(x) + (b-a-1).*log(1+x));
end

% abs(z) >= 1
id = (abs(z)>=1);
if any(id)
    zi = -1i./z(id);
    f(id,:) = exp(c + log(zi*o2) + (a-1)*log(zi*x) + ...
        (b-a-1)*log(1+zi*x) - o1(id)*x);
end

% z > 0 & z < 1
id = (z>0 & z<1);
if any(id)
    f(id,:) = exp(c + log(-1i) + (a-1)*log(-1i*o1(id)*x) + ...
        (b-a-1)*log(1-1i*o1(id)*x) - (z(id))*x);
end

% z < 0 & z > -1
id = (z<0 & z>-1);
if any(id)
    f(id,:) = exp(c + log(1i) + (a-1)*log(1i*o1(id)*x) + ...
        (b-a-1)*log(1+1i*o1(id)*x) + (z(id))*x);
end

end