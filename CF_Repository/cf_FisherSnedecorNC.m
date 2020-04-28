function cf = cf_FisherSnedecorNC(t,df1,df2,delta,coef,niid,tol)
%% cf_FisherSnedecorNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Fisher-Snedecor random variables, with
%  distributions F(df1_i,df2_i,delta_i). 
%
%  That is, cf_FisherSnedecorNC evaluates the characteristic function
%  cf(t) of  Y = coef_i*X_1 +...+ coef_N*X_N, where X_i ~
%  F(df1_i,df2_i,delta_i) are inedependent RVs, with df1_i and df2_i
%  degrees  of freedom, and the noncentrality parameters delta_i >0, for i
%  = 1,...,N.  
%
%  Random variable X has non-central F distribution with df1 and df2
%  degrees of freedom and the non-centrality parameter delta if X =
%  (X1/df1)/(X2/df2) where X1 ~ ChiSquare(df1,delta), i.e. with non-central
%  chi-square distribution, and X2 ~ ChiSquare(df2) is independent random
%  variable with central chi-square distriobution.
%
%  The characteristic function of X ~ F(df1,df2,delta) is Poisson  
%  mixture of the CFs of the scaled central F RVs of the form
%   cf(t) = cf_FisherSnedecorNC(t,df1,df2,delta) = 
%         = exp(-delta/2) sum_{j=1}^Inf (delta/2)^j/j! *
%           * cf_FisherSnedecor(t*(df1+2*j)/df1,df1+2*j,df2) 
%  where cf_FisherSnedecor(t,df1,df2) are the CFs of central F RVs with
%  parameters df1 and df2. Hence,the characteristic function of Y  =
%  coef(1)*Y1 + ... + coef(N)*YN is cf_Y(t) =  cf_Y1(coef(1)*t) * ... *
%  cf_YN(coef(N)*t), where cf_Yi(t) is evaluated with the parameters df1_i,
%  df2_i, and delta_i. 
%
% SYNTAX
%  cf = cf_FisherSnedecorNC(t,df1,df2,delta,coef,niid,tol)
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
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
%  tol   - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Noncentral_F-distribution
%
% EXAMPLE 1:
% % CF of the non-central F RV with delta = 1
%   df1   = 3;
%   df2   = 5;
%   delta = 1;
%   t     = linspace(-10,10,201);
%   cf    = cf_FisherSnedecorNC(t,df1,df2,delta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of non-central Fisher Snedecor RV')
%
% EXAMPLE 2:
% % CDF/PDF of the non-central F RV with delta = 1
%   df1   = 3;
%   df2   = 5;
%   delta = 1;
%   cf    = @(t) cf_FisherSnedecorNC(t,df1,df2,delta);
%   clear options;
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options)
%
% EXAMPLE 3:
% % CDF/PDF of the linear combination of non-central F RVs
%   df1   = [5 4 3];
%   df2   = [3 4 5];
%   delta = [0 1 2];
%   coef  = 1/3;
%   cf    = @(t) cf_FisherSnedecorNC(t,df1,df2,delta,coef);
%   clear options;
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],[],options)

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-Aug-2018 16:04:30
% Rev.: 28-Apr-2020 13:47:42

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
    cf = cf .* cf_ncFisherSnedecor(coef(i)*t,df1(i),df2(i),delta(i),tol);
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
function f = cf_ncFisherSnedecor(t,df1,df2,delta,tol)
% cf_ncFisherSnedecor Characteristic function of the distribution of
% the of the non-central Fisher-Snedecor distribution with df1 and df2
% degrees of freedom and the non-centrality parameter delta > 0.

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Feb-2018 13:53:40

%%
f = 0;
delta  = delta/2;
if delta == 0   % Deal with the central distribution
    f = cf_FisherSnedecor(t,df1,df2);
elseif delta > 0
    % Sum the Poisson series of CFs of independent central F RVs,
    % poisspdf(j,delta).*cf_FisherSnedecor(t*(df1+2*j)/df1,df1+2*j,df2)
    j0 = floor(delta/2);
    p0 = exp(-delta + j0 .* log(delta) - gammaln(j0 + 1));
    f  = f + p0 * cf_FisherSnedecor(t*(df1+2*j0)/df1,df1+2*j0,df2);
    p  = p0;
    j  = j0-1;
    while j >= 0 && p > tol
        p = p * (j+1) / delta;
        f = f + p * cf_FisherSnedecor(t*(df1+2*j)/df1,df1+2*j,df2);
        j = j - 1;
    end
    p  = p0;
    j  = j0+1;
    i  = 0;
    while p > tol && i <= 5000
        p = p * delta / j;
        f = f + p * cf_FisherSnedecor(t*(df1+2*j)/df1,df1+2*j,df2);
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