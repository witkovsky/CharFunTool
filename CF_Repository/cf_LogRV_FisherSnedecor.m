function cf = cf_LogRV_FisherSnedecor(t,df1,df2,coef,niid)
%% cf_LogRV_FisherSnedecor 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent  LOG-TRANSFORMED FISHER-SNEDECOR F random variables (RVs)
%  log(X), where X ~ F(df1,df2) has the FISHER-SNEDECOR F distribution with
%  df1 > 0 and df2 > 0 degrees of freedom.
%  
%  That is, cf_LogRV_FisherSnedecor evaluates the characteristic function
%  cf(t) of  Y = coef_1*log(X_1) +...+ coef_N*log(X_N), where X_i ~
%  F(df1_i,df2_i), with degrees of freedom df1_i and df2_i, for i =
%  1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ F(df1,df2,lambda),
%  where lambda is the non-centrality parameter, is defined by cf_Y(t) =
%  E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)).  
%  That is, the characteristic function can be derived from expression for
%  the r-th moment of X, E(X^r) by using (1i*t) instead of r. In
%  particular, the characteristic function of Y = log(X) is defined by  
%   cf_Y(t) = (df2/df1)^(1i*t) .* gamma(df1/2 + 1i*t) / gamma(df1/2) .* ...
%             gamma(df2/2 - 1i*t) / gamma(df2/2).
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters df1(i), df2(i), and lambda(i).
%
% SYNTAX
%  cf = cf_LogRV_FisherSnedecor(t,df1,df2,coef,niid)
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
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
%
%  WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/F-distribution 
%
% EXAMPLE 1:
% % CF of a weighted linear combination of independent log-F RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   df1 = 5;
%   df2 = 3;
%   t   = linspace(-10,10,201);
%   cf  = cf_LogRV_FisherSnedecor(t,df1,df2,weight);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('Characteristic function of a linear combination of log-F RVs')
%
% EXAMPLE 2:
% % PDF/CDF from the CF by cf2DistGP
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   df1 = 5;
%   df2 = 3;
%   cf    = @(t) cf_LogRV_FisherSnedecor(t,df1,df2,weight);
%   clear options
%   options.N = 2^12;
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-Aug-2018 15:46:49

%% ALGORITHM
% cf = cf_LogRV_FisherSnedecor(t,df1,df2,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, df2 = []; end
if nargin < 2, df1 = []; end

%%
if isempty(df2) && ~isempty(df1)
    df2 = 1;
elseif isempty(df2) && ~isempty(coef)
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

if isempty(niid)
    niid = 1;
end


%% Check size of the parameters
[errorcode,coef,df1,df2] = distchck(3,coef(:)',df1(:)',df2(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);
aux = 1i*t*coef;
aux = GammaLog(bsxfun(@plus,aux,df1/2))  + ...
      GammaLog(bsxfun(@plus,-aux,df2/2)) - ...
      ones(length(t),1)*(GammaLog(df1/2) + GammaLog(df2/2)) + ...
      bsxfun(@times,aux,log(df2./df1));
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