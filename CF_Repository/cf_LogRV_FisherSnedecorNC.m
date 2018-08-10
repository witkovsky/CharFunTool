function cf = cf_LogRV_FisherSnedecorNC(t,df1,df2,delta,coef,niid,tol)
%% cf_LogRV_FisherSnedecorNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED non-central Fisher-Snedecor random
%  variables, with distributions F(df1_i,df2_i,delta_i).
%
%  That is, cf_LogRV_FisherSnedecorNC evaluates the characteristic function
%  cf(t) of  Y = coef_i*log(X_1) +...+ coef_N*log(X_N), where X_i ~
%  F(df1_i,df2_i,delta_i) are inedependent RVs, with df1_i and df2_i
%  degrees of freedom, and the noncentrality parameters delta_i >0, for i =
%  1,...,N.  
%
%  The characteristic function of Y = log(X) with X ~ F(df1,df2,delta) is
%  Poisson mixture of the CFs of the shifted log-transformed central F RVs
%  of the form  
%   cf(t) = cf_LogRV_FisherSnedecorNC(t,df1,df2,delta) = 
%         = exp(-delta/2) sum_{j=1}^Inf (delta/2)^j/j! .*
%           exp(1i*t*(df1+2*j)/df1) .* 
%           cf_LogRV_FisherSnedecor(t,df1+2*j,df2) 
%  where cf_LogRV_FisherSnedecor(t,df1,df2) denotes CF of log-transformed
%  centrally distributed F RVs with parameters df1 and df2. For more
%  details on  the non-central F distribution see cf_FisherSnedecorNC. 
%  Alternatively, 
%   cf(t) = (df2/df1)^(1i*t) .* gamma(df1/2 + 1i*t) / gamma(df1/2) .* ...
%           gamma(df2/2 - 1i*t) / gamma(df2/2) .* ...
%           1F1(-1i*t;df1/2;-delta/2),
%  where 1F1(a;b;z) is the confluent hypergeometric function, also known as
%  the Kummer function M(a,b,z). 
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN 
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters df1_i, df2_i, and delta_i.
%
% SYNTAX
%  cf = cf_LogRV_FisherSnedecorNC(t,df1,df2,delta,coef,niid,tol)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df1   - vector of the  degrees of freedom df1 > 0. If empty, default
%          value is df1 = 1.  
%  df2   - vector of the  degrees of freedom df2 > 0. If empty, default
%          value is df2 = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Noncentral_F-distribution
%
% EXAMPLE 1:
% % CF of the log-transformed non-central F RV with delta = 1 and coef = -1
%   df1   = 3;
%   df2   = 5;
%   delta = 1;
%   coef  = -1;
%   t = linspace(-10,10,201);
%   cf = cf_LogRV_FisherSnedecorNC(t,df1,df2,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of minus log-transformed F RV')
%
% EXAMPLE 2:
% % CDF/PDF of the minus log-transformed non-central F RV with delta = 1
%   df1   = 3;
%   df2   = 5;
%   delta = 1;
%   coef  = -1;
%   cf = @(t) cf_LogRV_FisherSnedecorNC(t,df1,df2,delta,coef);
%   clear options;
%   options.N  = 2^12;
%   result = cf2DistGP(cf,[],[],options)
%
% EXAMPLE 3:
% % CDF/PDF of the linear combination of log-transformed non-central F RVs
%   df1   = [5 4 3];
%   df2   = [3 4 5];
%   delta = [0 1 2];
%   coef  = -1/3;
%   cf = @(t) cf_LogRV_FisherSnedecorNC(t,df1,df2,delta,coef);
%   clear options;
%   options.N  = 2^12;
%   result = cf2DistGP(cf,[],[],options)
%
% SEE ALSO: cf_FisherSnedecorNC

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-Aug-2018 15:46:49

%% CHECK THE INPUT PARAMETERS
narginchk(1, 7);

if nargin < 7, tol = []; end
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, delta = []; end
if nargin < 3, df2 = []; end
if nargin < 2, df1 = []; end

if isempty(df1), df1 = 1; end
if isempty(df2), df2 = 1; end
if isempty(delta), delta = 0; end
if isempty(coef), coef = 1; end
if isempty(tol), tol = 1e-12; end

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,df1,df2,delta] = ...
    distchck(4,coef(:),df1(:),df2(:),delta(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent nc F RVs
szt = size(t);
t   = t(:);
cf  = ones(length(t),1);
for i = 1:length(coef)
    cf = cf .* cf_ncLogRVFisherSnedecor(coef(i)*t,df1(i),df2(i),delta(i),tol);
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
function f = cf_ncLogRVFisherSnedecor(t,df1,df2,delta,tol)
% cf_ncLogRVFisherSnedecor Characteristic function of the distribution of
% the of the distribution of the log-transformed non-central
% Fisher-Snedecor RV with df1 and df2 degrees of freedom and the
% non-centrality parameter delta > 0. 

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Feb-2018 13:53:40

%%
f = 0;
delta  = delta/2;
if delta == 0   % Deal with the central distribution
    f = cf_LogRV_FisherSnedecor(t,df1,df2);
elseif delta > 0
    % Sum the Poisson series of CFs of independent log-transformed F RVs,
    % poisspdf(j,delta).*cf_LogRV_FisherSnedecor(t*(df1+2*j)/df1,df1+2*j,df2)
    j0 = floor(delta/2);
    p0 = exp(-delta + j0 .* log(delta) - gammaln(j0 + 1));
    f  = f + p0 * exp(1i*t*log((df1+2*j0)/df1)) .* ... 
         cf_LogRV_FisherSnedecor(t,df1+2*j0,df2);
    p  = p0;
    j  = j0-1;
    while j >= 0 && p > tol
        p = p * (j+1) / delta;
        f = f + p * exp(1i*t*log((df1+2*j)/df1)) .* ... 
            cf_LogRV_FisherSnedecor(t,df1+2*j,df2);
        j = j - 1;
    end
    p  = p0;
    j  = j0+1;
    i  = 0;
    while p > tol && i <= 5000
        p = p * delta / j;
        f = f + p * exp(1i*t*log((df1+2*j)/df1)) .* ... 
            cf_LogRV_FisherSnedecor(t,df1+2*j,df2);
        j = j + 1;
        i = i + 1;
    end
    if (i == 5000)
        warning(message('NoConvergence'));
    end
else
    error('delta should be nonnegative');
end
end