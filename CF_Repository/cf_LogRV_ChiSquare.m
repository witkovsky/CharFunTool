function cf = cf_LogRV_ChiSquare(t,df,coef,niid)
%% cf_LogRV_ChiSquare 
%  Characteristic function of a linear combination (resp.
%  convolution) of independent LOG-TRANSFORMED CHI-SQUARE random variables
%  (RVs) log(X), where X ~ ChiSquare(df) is central CHI-SQUARE distributed
%  RV with df degrees of freedom.
%  
%  That is, cf_LogRV_ChiSquare evaluates the characteristic function cf(t)
%  of Y = coef_1*log(X_1) +...+ coef_N*log(X_N), where X_i ~
%  ChiSquare(df_i), with df_i > 0 degrees of freedom, for i = 1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ ChiSquare(df) is
%  defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)). 
%  That is, the characteristic function can be derived from expression for
%  the r-th moment of X, E(X^r) by using (1i*t) instead of r. In
%  particular, the characteristic function of Y = log(X) is
%   cf_Y(t) = 2^(1i*t) * gamma(df/2 + 1i*t) / gamma(df/2).
%
%  Hence,the characteristic function of Y  = coef_1*X_1 +...+ coef_N*X_N
%  is  cf_Y(t) =  cf_1(coef_1*t) *...* cf_N(coef_N*t), where cf_i(t)
%  is the characteristic function of the ChiSquare distribution with df_i
%  degrees of freedom. 
%
% SYNTAX
%  cf = cf_LogRV_ChiSquare(t,df,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the degrees of freedom of the the chi-squared random
%          variables.  If df is scalar, it is assumed that all degrees of
%          freedom are equal. If empty, default value is df = 1. 
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed chi-square random variables. If coef is scalar,
%          it is assumed that all coefficients are equal. If empty, default
%          value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.    
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Chi-squared_distribution
%
% EXAMPLE 1:
% % CF of a weighted linear combination of independent log-ChiSquare RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   df     = [1 2 3 4 5];
%   t      = linspace(-10,10,501);
%   cf     = cf_LogRV_ChiSquare(t,df,weight);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a linear combination of minus log-ChiSquare RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a linear combination of independent log-ChiSquare RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   df     = [1 2 3 4 5];
%   cf     = @(t) cf_LogRV_ChiSquare(t,df,weight);
%   x      = linspace(1,6);
%   prob   = [0.9 0.95 0.99];
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 15:44:23

%% ALGORITHM
% cf = cf_LogRV_ChiSquare(t,df,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid = []; end
if nargin < 3, coef = []; end
if nargin < 2, df = []; end

%% SET the default values 
if isempty(df) && ~isempty(coef)
    df = 1;
end

if isempty(coef) && ~isempty(df)
    coef = 1;
end

if isempty(niid)
    niid = 1;
end

%% Check size of the parameters
[errorcode,coef,df] = distchck(2,coef(:)',df(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of linear combination 
szt = size(t);
t   = t(:);
aux = 1i*t*coef;
aux = GammaLog(bsxfun(@plus,aux,df/2))-ones(length(t),1)*GammaLog(df/2);
aux = aux + 1i*t*log(2);
cf  = prod(exp(aux),2);
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