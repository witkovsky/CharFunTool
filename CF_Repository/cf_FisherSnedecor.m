function cf = cf_FisherSnedecor(t,df1,df2,coef,niid)
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
%           KummerU(df1/2, 1-df2/2, -1i*(df2/df1)*t),
%  where KummerU(a,b,z) denotes the confluent hypergeometric function of
%  the second kind U(a,b,z).
%
%  Hence, the characteristic function of Y  = coef_1*X_1 +...+ coef_N*X_N
%  is cf_Y(t) =  cf_1(coef_1*t) *...* cf_N(coef_N*t), where cf_i(t) 
%  is the characteristic function of X_i ~ F(df1_i,df2_i).
%
% SYNTAX
%  cf = cf_FisherSnedecor(t,df1,df2,coef,niid)
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
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/F-distribution.
%
% EXAMPLE 1: 
%   % CF of a linear combination of independent F RVs
%   df1   = 3:12;
%   df2   = 14:-1:5;
%   coef  = 1/10;
%   t     = linspace(-15,15,501);
%   cf    = cf_FisherSnedecor(t,df1,df2,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on
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

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Feb-2023 11:25:10

%% ALGORITHM
% cf = cf_FisherSnedecor(t,df1,df2,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
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

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,df1,df2] = distchck(3,coef(:)',df1(:)',df2(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent F RVs
szt = size(t);
t   = t(:);
o   = ones(size(t));

%  a    = o * (df1/2);
%  b    = o * (1-df2/2);
%  z    = -1i * t * (coef .* df2 ./df1);
%  cf = prod(gamma((df1+df2)/2) ./ gamma(df2/2)) * prod(KummerU(a,b,z),2);

cf  = prod(gamma(df1/2+df2/2)./gamma(df2/2));
cf  = cf .* prod(KummerU(o*(df1/2), o*(1-df2/2), -1i*t*(coef.*df2./df1)),2);

cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end