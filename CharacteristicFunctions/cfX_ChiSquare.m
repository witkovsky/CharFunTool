function cf = cfX_ChiSquare(t,df,ncp)
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
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  df    - the degrees of freedom parameter df > 0. If empty, default value
%          is df = 1.   
%  ncp   - the non-centraloty parameter ncp > 0. If empty, default value
%          is ncp = 0.   
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

%% ALGORITHM
narginchk(1, 3);
if nargin < 3, ncp  = []; end
if nargin < 2, df   = []; end
if isempty(df),   df = 1; end
if isempty(ncp), ncp = 0; end

cf = cf_ChiSquare(t,df,ncp);

end