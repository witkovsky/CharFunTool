function cf = cf_LogRV_Chi(t,df,coef,niid)
%% cf_LogRV_Chi 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED CHI random variables (RVs) log(X), where X ~
%  Chi(df) is central CHI distributed RV with df degrees of freedom.
%  
%  That is, cf_LogRV_Chi evaluates the characteristic function cf(t) of Y =
%  coef_1*log(X_1) +...+ coef_N*log(X_N), where X_i ~ Chi(df_i), with df_i
%  > 0 degrees of freedom, for i = 1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ Chi(df) is
%  defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)). 
%  That is, the characteristic function can be derived from expression for
%  the r-th moment of X, E(X^r) by using (1i*t) instead of r. In
%  particular, the characteristic function of Y = log(X) is
%   cf_Y(t) = 2^(1i*t/2) * gamma(df/2 + 1i*t/2) / gamma(df/2).
%
%  Hence,the characteristic function of Y  = coef_1*X_1 +...+ coef_N*X_N is
%  cf_Y(t) =  cf_1(coef_1*t) *...* cf_N(coef_N*t), where cf_i(t) is the
%  characteristic function of the Chi distribution with df_i degrees of
%  freedom.
%
% SYNTAX
%  cf = cf_LogRV_Chi(t,df,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the degrees of freedom of the the chi-distributes
%          random variables.  If df is scalar, it is assumed that all
%          degrees of freedom are equal. If empty, default value is df = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed chi-distributed random variables. If coef is
%          scalar, it is assumed that all coefficients are equal. If empty,
%          default value is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.    
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Chi_distribution
%
% EXAMPLE 1:
% % CF of a log-transformed Chi RVs
%   df     = 3;
%   t      = linspace(-20,20,501);
%   cf     = cf_LogRV_Chi(t,df);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a log-transformed Chi RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a log-transformed Chi RVs
%   df     = 3;
%   cf     = @(t) cf_LogRV_Chi(t,df);
%   x      = linspace(-3,2);
%   prob   = [0.9 0.95 0.99];
%   clear options
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 3:
% % CF of a linear combination of independent log-transformed Chi RVs
%   coef   = [1 2 3 4 5];
%   coef   = coef/sum(coef);
%   df     = [1 2 3 4 5];
%   t      = linspace(-20,20,501);
%   cf     = cf_LogRV_Chi(t,df,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a combination of independent log-transformed Chi RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent log-transformed Chi RVs
%   coef   = [1 2 3 4 5];
%   coef   = coef/sum(coef);
%   df     = [1 2 3 4 5];
%   cf     = @(t) cf_LogRV_Chi(t,df,coef);
%   x      = linspace(0.5,3);
%   prob   = [0.9 0.95 0.99];
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% SEE ALSO: cf_LogRV_ChiSquare

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 15:44:23

%% ALGORITHM
% cf = cf_LogRV_Chi(t,df,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid = []; end
if nargin < 3, coef = []; end
if nargin < 2, df = []; end

% CF of the linear combination of independent LOG-TRANSFORMED
% CHI-DISTIBUTED RVs expressed as a special case of the cf_LogRV_ChiSquare
cf = cf_LogRV_ChiSquare(t/2,df,coef,niid);

end