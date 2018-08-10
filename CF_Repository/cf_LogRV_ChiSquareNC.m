function cf = cf_LogRV_ChiSquareNC(t,df,delta,coef,niid,tol)
%% cf_LogRV_ChiSquareNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED non-central ChiSquare random variables, with
%  distributions ChiSquare(df_i,delta_i).
%
%  That is, cf_LogRV_ChiSquareNC evaluates the characteristic function cf(t)
%  of Y = coef_i*log(X_1) +...+ coef_N*log(X_N), where X_i ~
%  ChiSquare(df_i,delta_i) are inedependent RVs, with df_i degrees of
%  freedom and the noncentrality parameters delta_i >0, for i = 1,...,N. 
%
%  The characteristic function of Y = log(X) with X ~
%  ChiSquare(df,delta) is Poisson mixture of CFs of the central
%  log-transformed ChiSquare RVs of the form 
%   cf(t) = cf_LogRV_ChiSquareNC(t,df,delta) = 
%         = exp(-delta/2) sum_{j=1}^Inf (delta/2)^j/j! *
%           * cf_LogRV_ChiSquare(df+2*j) 
%  where cf_LogRV_ChiSquare(df+j) are the CFs of central log-transformed
%  ChiSquare RVs with df+j degrees of freedom. 
%  Alternatively, 
%   cf(t) = cf_LogRV_ChiSquareNC(t,df,delta) = 
%         = cf_LogRV_ChiSquare(t,df) * 1F1(-1i*t;df/2;-delta/2),
%  where 1F1(a;b;z) is hypergeometric function. 
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters df_i and delta_i.
%
% SYNTAX
%  cf = cf_LogRV_ChiSquareNC(t,df,delta,coef,niid,tol)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the  degrees of freedom df > 0. If empty, default
%          value is df = 1.  
%  delta - vector of the non-centrality parameters delta > 0. If empty,
%          default value is delta = 0.
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
% tol    - tolerance factor for selecting the Poisson weights, i.e. such
%          that PoissProb > tol. If empty, default value is tol = 1e-12.
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution
%
% EXAMPLE 1:
% % CF of the minus log-transformed non-central ChiSquare RV with delta = 1
%   df    = 3;
%   delta = 1;
%   coef  = -1;
%   t = linspace(-10,10,201);
%   cf = cf_LogRV_ChiSquareNC(t,df,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of minus log-transformed ChiSquare RV')
%
% EXAMPLE 2:
% % CDF/PDF of the minus log-transformed non-central ChiSquare RV
%   df    = 3;
%   delta = 1;
%   coef  = -1;
%   cf = @(t) cf_LogRV_ChiSquareNC(t,df,delta,coef);
%   clear options;
%   options.N = 2^10;
%   result = cf2DistGP(cf,[],[],options)
%
% EXAMPLE 3:
% % CDF/PDF of the linear combination of the minus log-transformed
% % non-central ChiSquare RVs 
%   df    = [3 4 5];
%   delta = [0 1 2];
%   coef  = [-1 -1 -1]/3;
%   cf = @(t) cf_LogRV_ChiSquareNC(t,df,delta,coef);
%   clear options;
%   options.N = 2^10;
%   result = cf2DistGP(cf,[],[],options)

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Feb-2018 13:53:40

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);

if nargin < 6, tol = []; end
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, delta = []; end
if nargin < 2, df = []; end

if isempty(df), df = 1; end
if isempty(delta), delta = 0; end
if isempty(coef), coef = 1; end
if isempty(tol), tol = 1e-12; end

%% SET THE COMMON SIZE of the parameters
[errorcode,coef,df,delta] = distchck(3,coef(:),df(:),delta(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent nc RVs
szt = size(t);
t   = t(:);
cf  = ones(length(t),1);
for i = 1:length(coef)
    cf = cf .* cf_ncLogRVChiSquare(coef(i)*t,df(i),delta(i),tol);
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
function f = cf_ncLogRVChiSquare(t,df,delta,tol)
% cf_ncLogRVChiSquare Characteristic function of the distribution of the
% non-central log-transformed ChiSquare RV with df degrees of freedom
% and the non-centrality parameter delta > 0. 

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Feb-2018 13:53:40

%%
f = 0;
delta  = delta/2;
if delta == 0   % Deal with the central distribution
    f = cf_LogRV_ChiSquare(t,df);
elseif delta > 0
    % Sum the Poisson series of CFs of central iid ChiSquare RVs,
    % poisspdf(j,delta) .* cf_LogRV_ChiSquare(t,alpha+j,beta)
    j0 = floor(delta/2);
    p0 = exp(-delta + j0 .* log(delta) - gammaln(j0 + 1));
    f  = f + p0 * cf_LogRV_ChiSquare(t,df+2*j0);
    p  = p0;
    j  = j0-1;
    while j >= 0 && p > tol
        p = p * (j+1) / delta;
        f = f + p * cf_LogRV_ChiSquare(t,df+2*j);
        j = j - 1;
    end
    p  = p0;
    j  = j0+1;
    i  = 0;
    while p > tol && i <= 5000
        p = p * delta / j;
        f = f + p * cf_LogRV_ChiSquare(t,df+2*j);
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