function cf = cf_LogRV_ChiNC(t,df,delta,coef,niid,tol)
%% cf_LogRV_ChiNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED NON-CENTAL CHI random variables, with
%  distributions Chi(df_i,delta_i), with df_i > 0 degrees of freedom and
%  the noncentrality parameters delta_i >= 0, for all i = 1,...,N.
% 
%  That is, cf_LogRV_ChiNC evaluates the characteristic function cf(t) of Y
%  = coef_i*log(X_1) +...+ coef_N*log(X_N), where X_i ~ Chi(df_i,delta_i)
%  are inedependent RVs, with df_i degrees of freedom and the noncentrality
%  parameters delta_i >0, for i = 1,...,N.
% 
%  The NON-CENTAL CHI distributions is related to the NON-CENTAL CHI-SQUARE
%  distribution. If X is a random variable with the non-central chi
%  distribution, X ~ Chi(df,delta), the random variable X^2 will have the
%  noncentral chi-square distribution, X^2 ~ Chi(df,delta^2). Hence, log(X)
%  = (1/2)*log(X^2).
%  
%  Notice that in definition of the non-central chi-distribution here the
%  noncentrality parameter delta is defined as a square root of the
%  non-centrality parameter of the non-central ChiSquare distribution, 
%  (delta = sqrt(sum_{i=1}^df mu_i^2)).
%
%  From that, the characteristic function of log(X) with X ~ Chi(df,delta)
%  is  
%   cf(t) = cf_LogRV_ChiNC(t,df,delta) = 
%         = cf_LogRV_ChiSquareNC(t/2,df,delta^2)
%  where cf_LogRV_ChiSquareNC is the CFs of log-transformed non-central
%  ChiSquare dtribution.
% 
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters df_i and delta_i.
%
% SYNTAX
%  cf = cf_LogRV_ChiNC(t,df,delta,coef,niid,tol)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - vector of the  degrees of freedom df > 0. If empty, default
%          value is df = 1.  
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0. Notice that the noncentrality
%          parameter delta can be interpreted as a square root of the sum
%          of squared means, delta = sqrt(sum_{i=1}^df mu_i^2).
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
%   https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%
% EXAMPLE 1:
% % CF of the minus log-transformed non-central Chi RV with delta = 1
%   df    = 3;
%   delta = 1;
%   coef  = -1;
%   t     = linspace(-10,10,201);
%   cf    = cf_LogRV_ChiNC(t,df,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of minus log-transformed Chi RV')
%
% EXAMPLE 2:
% % CDF/PDF of the minus log-transformed non-central Chi RV
%   df    = 3;
%   delta = 1;
%   coef  = -1;
%   cf    = @(t) cf_LogRV_ChiNC(t,df,delta,coef);
%   x     = linspace(-2,2);
%   prob  = [0.9,0.95,0.99];
%   clear options;
%   options.N = 2^10;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3:
% % CDF/PDF of the linear combination of the minus log-transformed
% % non-central Chi RVs 
%   df    = [3 4 5];
%   delta = [0 1 2];
%   coef  = [-1 -1 -1]/3;
%   cf = @(t) cf_LogRV_ChiNC(t,df,delta,coef);
%   x     = linspace(-1.5,0.5);
%   prob  = [0.9,0.95,0.99];
%   clear options;
%   options.N = 2^10;
%   result = cf2DistGP(cf,x,prob,options)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 16:48:55

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);

if nargin < 6, tol   = []; end
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 3, delta = []; end
if nargin < 2, df    = []; end

cf = cf_LogRV_ChiSquareNC(t/2,df,delta.^2,coef,niid,tol);

end
