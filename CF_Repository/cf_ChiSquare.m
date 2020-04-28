function cf = cf_ChiSquare(t,df,ncp,coef,niid)
%cf_ChiSquare
%  Characteristic function of a linear combination (resp.
%  convolution) of independent (possibly non-central) CHI-SQUARE random
%  variables.
%
%  That is, cf_ChiSquare evaluates the characteristic function cf(t)  of  Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ ChiSquare(df_i,ncp_i) are
%  inedependent RVs, with df_i > 0 degrees of freedom the 'non-centrality'
%  parameters ncp_i > 0, for i = 1,...,N.
%
%  The characteristic function of Y is defined by
%   cf(t) = Prod ( (1-2*i*t*coef(i))^(-df(i)/2)
%                   * exp((i*t*ncp(i))/(1-2*i*t*coef(i))) )
%
% SYNTAX:
%  cf = cf_ChiSquare(t,df,ncp,coef,niid)
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
%  coef  - vector of the coefficients of the linear combination of the
%          chi-squared random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Chi-squared_distribution
%  https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution
%
% EXAMPLE 1:
% % CF of the distribution of ChiSquare RV with DF = 1, NCP = 1
%   df   = 1;
%   ncp  = 1;
%   coef = 1;
%   t    = linspace(-100,100,501);
%   cf   = cf_ChiSquare(t,df,ncp,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the \chi^2 RV with DF=1 and NCP=1')
%
% EXAMPLE 2:
% % CF of a linear combination of K=100 independent ChiSquare RVs
%   coef = 1./(((1:50) - 0.5)*pi).^2;
%   figure; plot(coef,'.-')
%   title('Coefficients of the linear combination of \chi^2 RVs with DF=1')
%   df   = 1;
%   cf   = cf_ChiSquare(t,df,[],coef);
%   t    = linspace(-100,100,501);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the linear combination of \chi^2_1 RVs')
%
% EXAMPLE 3:
% % PDF/CDF from the CF by cf2DistGP
%   df   = 1;
%   coef = 1./(((1:50) - 0.5)*pi).^2;
%   cf = @(t) cf_ChiSquare(t,df,[],coef);
%   clear options
%   options.N = 2^12;
%   options.xMin = 0;
%   x = linspace(0,3.5,501);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% REFERENCES:
%  IMHOF J. (1961): Computing the distribution of quadratic forms in normal
%  variables. Biometrika 48, 419–426.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-May-2017 18:11:50
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf4Gamma(t,alpha,beta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);

if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, ncp  = []; end
if nargin < 2, df   = []; end

if isempty(ncp)
    ncp = 0;
end

if isempty(df)
    df = 1;
end

if isempty(coef)
    coef = 1;
end

if isempty(niid)
    niid = 1;
end

%% Find the unique coefficients and their multiplicities
%     if ~isempty(coef) && isscalar(df) && ~isNoncentral
%         coef = sort(coef);
%         m    = length(coef);
%         [coef,idx] = unique(coef);
%         df = df * diff([idx;m+1]);
%     end

% Check/set equal dimensions for the vectors coef, df, and ncp
if ncp~=0
    [errorcode,coef,df,ncp] = distchck(3,coef(:)',df(:)',ncp(:)');
else
    [errorcode,coef,df] = distchck(2,coef(:)',df(:)');
end
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt   = size(t);
t     = t(:);
cf    = 1;
aux   = bsxfun(@times,t,coef);
if ncp~=0
    aux = bsxfun(@power,(1 -2i * aux),-df/2) .* ...
        exp(bsxfun(@times,aux,1i*ncp)./(1-2i*aux));
else
    aux = bsxfun(@power,(1 -2i * aux),-df/2);
end
cf   = cf .* prod(aux,2);
cf   = reshape(cf,szt);

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end