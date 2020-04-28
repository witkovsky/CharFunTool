function cf = cfX_ChiSquare(t,df,ncp,coef,niid)
%% cfX_ChiSquare
%  Characteristic function of the CHI-SQUARE distribution with df > 0
%  degrees of freedom and the non-cetrality parameter ncp > 0.
%
%  cfX_ChiSquare is an ALIAS NAME of the more general function
%  cf_ChiSquare, used to evaluate the characteristic function of a linear 
%  combination of independent CHI-SQUARE distributed random variables.
%
%  The characteristic function of the CHI-SQUARE distribution is defined by
%   cf(t) = (1 - 2*i*t )^(df/2) * exp((i*t*ncp)/(1-2*i*t)).
%
% SYNTAX
%  cf = cfX_ChiSquare(t,df,ncp)
%  cf = cfX_ChiSquare(t,df,ncp,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - the degrees of freedom parameter df > 0. If empty, default value
%          is df = 1.   
%  ncp   - the non-centraloty parameter ncp > 0. If empty, default value
%          is ncp = 0. 
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1. 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Chi-squared_distribution.
%
% EXAMPLE 1: 
%  % CF of the ChiSquared distribution with df = 1
%  df = 1;
%  t = linspace(-50,50,501);
%  cf = cfX_ChiSquare(t,df);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Chi-square distribution with df = 1')
%
% EXAMPLE 2:
%  % PDF/CDF of the ChiSquare distribution with df = 3
%  df = 3;
%  x = linspace(0,15,101);
%  prob = [0.9 0.95 0.99];
%  cf = @(t) cfX_ChiSquare(t,df);
%  clear options
%  options.xMin = 0;
%  options.N = 2^14;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3:
%  % PDF/CDF of the compound Binomial-ChiSquared distribution
%  n = 25;  
%  p = 0.3;
%  df = 3;
%  cfX = @(t) cfX_ChiSquare(t,df);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,80,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jun-2017 10:07:43
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, ncp  = []; end
if nargin < 2, df   = []; end

cf = cf_ChiSquare(t,df,ncp,coef,niid);

end