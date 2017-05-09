function cf = cfX_ChiSquared(t,df,ncp)
%cfX_ChiSquared evaluates the characteristic function cf(t) of the
% (non-central) CHI-SUQARED distribution with the parameters df (degrees of
% freedom, df >= 1) and ncp (noncentrality parameter), i.e. 
%   cf(t) = cfX_ChiSquared(t,df,ncp) 
%         = (1 - 2*i*t )^(df/2) * exp((i*t*ncp)/(1-2*i*t)).
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Chi-squared_distribution
%
% SYNTAX
%  cf = cfX_ChiSquared(t,df,ncp)
%
% EXAMPLE1 (CF of the ChiSquared distribution with df = 1)
%  df = 1;
%  t = linspace(-50,50,501);
%  cf = cfX_ChiSquared(t,df);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Chi-squared distribution with df = 1')
%
% EXAMPLE2 (PDF/CDF of the ChiSquared distribution with df = 3)
%  df = 3;
%  x = linspace(0,15,101);
%  prob = [0.9 0.95 0.99];
%  cf = @(t) cfX_ChiSquared(t,df);
%  clear options
%  options.xMin = 0;
%  options.N = 2^14;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-ChiSquared distribution)
%  n = 25;  
%  p = 0.3;
%  df = 3;
%  cfX = @(t) cfX_ChiSquared(t,df);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,80,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)

% Copyright (c) 2017, Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 9-May-2017 10:22:48

%% ALGORITHM
%cf = cfX_ChiSquared(t,df,ncp)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, ncp = []; end
if nargin < 2, df  = []; end

%%
if isempty(ncp)
    ncp = 0;
end

if isempty(df)
    df = 1;
end
%% Characteristic function of the ChiSquared distribution
szt = size(t);
t   = t(:);

if ncp ~= 0
    cf = exp((1i*t*ncp)./(1-2i*t));
else
    cf = 1;
end

cf  = cf .* (1-2i*t).^(-df/2);
cf  = reshape(cf,szt);
cf(t==0) = 1;

end