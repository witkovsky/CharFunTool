function cf = cf4ChiSquaredConv(t,df,ncp,coef)
%cf4ChiSquaredConv Computes the characteristic function of the distribution
%  of a linear combination of independent chi-square random variables, i.e. of
%  the distribution of Y = sum( coef[i] * ChiSquareRV_{df[i],ncp[i]} ).
%  In particular, the characteristic function of Y is defined by
%   cf(t) = Prod ( (1-2*i*t*coef[i])^(-df[i]/2)
%                             * exp((i*t*ncp[i])/(1-2*i*t*coef[i])) )
%
%  For more details see, e.g., WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Chi-squared_distribution
%  https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution
%
% SYNTAX
%  cf = cf4ChiSquaredConv(t,df,ncp,coef)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the degrees of freedom of the the chi-squared random
%          variables.  If df is scalar, it is assumed that all degrees of
%          freedom are equal. If empty, default value is df = 1.
%  ncp   - vector of the non-centrality parameters of the the chi-squared
%          random variables. If ncp is scalar, it is assumed that all
%          non-centrality parameters are equal. If empty, default value is
%          ncp = 0. 
%  coef  - vector of the coefficients of the lienar combination of the
%          chi-squared random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%
% EXAMPLE1 (CF of the distribution of chi-square RV with DF = 1, NCP = 1)
%  df   = 1;
%  ncp  = 1;
%  coef = 1;
%  t    = linspace(-100,100,501);
%  cf   = cf4ChiSquaredConv(t,df,ncp,coef);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the \chi^2 RV with DF=1 and NCP=1')
%
% EXAMPLE2 (CF of a linear combination of K=100 independent chi2_1 RVs)
%  K    = 100;
%  idx  = 1:K;
%  coef = 1./((idx - 0.5)*pi).^2;
%  figure; plot(idx,coef,'.-')
%  title('Coefficients of the linear combination of \chi^2 RVs with DF=1')
%  df   = 1;
%  cf   = cf4ChiSquaredConv(t,df,[],coef);
%  t    = linspace(-100,100,501);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the linear combination of \chi^2_1 RVs')
%
% REFERENCES
%  Imhof, J.: Computing the distribution of quadratic forms in normal
%  variables. Biometrika 48, 419–426 (1961).

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-Apr-2017 18:10:40

%% Check and set the default input parameters
narginchk(1, 4);
if nargin < 4, coef = []; end
if nargin < 3, ncp  = []; end
if nargin < 2, df   = []; end

if isempty(ncp) && ~isempty(df)
    isNoncentral = false; ncp = 0;
elseif isempty(ncp) && ~isempty(coef)
    isNoncentral = false; ncp = 0;
elseif ~any(ncp)
    isNoncentral = false; ncp = 0;
else
    isNoncentral = true;
end

if isempty(df) && ~isempty(coef)
    df = 1;
elseif isempty(df) && ~isempty(ncp)
    df = 1;
end

if isempty(coef) && ~isempty(ncp)
    coef = 1;
elseif isempty(coef) && ~isempty(df)
    coef = 1;
end

%% Find the unique coefficients and their multiplicities
if ~isempty(coef) && isscalar(df) && ~isNoncentral 
    coef = sort(coef);
    m    = length(coef);
    [coef,idx] = unique(coef);
    df = df * diff([idx;m+1]);
end

% Check/set equal dimensions for the vectors coef, df, and ncp
[errorcode,coef,df,ncp] = distchck(3,coef(:)',df(:)',ncp(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% ALGORITHM 
% Characteristic function of a linear combination of noncentral chi-squares

szt   = size(t);
t     = t(:);
cf    = 1;
aux   = bsxfun(@times,t,coef);
if isNoncentral
    aux = bsxfun(@power,(1 -2i * aux),-df/2) .* ...
        exp(bsxfun(@times,aux,1i*ncp)./(1-2i*aux));
else
    aux = bsxfun(@power,(1 -2i * aux),-df/2);
end
cf   = cf .* prod(aux,2);
cf   = reshape(cf,szt);

end