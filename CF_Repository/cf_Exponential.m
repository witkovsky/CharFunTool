function cf = cf_Exponential(t,lambda,coef,niid)
%% cf_Exponential 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent EXPONENTIAL random variables.
%
%  That is, cf_Exponential evaluates the characteristic function cf(t) of
%  Y = sum_{i=1}^N coef_i * X_i, where X_i ~ EXP(lambda_i) are inedependent
%  RVs, with the rate parameters lambda_i > 0, for i = 1,...,N.
%
%  The characteristic function of Y is defined by
%   cf(t) = Prod( lambda_i / (lambda_i - 1i*t) )
%
% SYNTAX:
%  cf = cf_Exponential(t,lambda,coef,niid)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  lambda - vector of the 'rate' parameters lambda > 0. If empty, default
%           value is lambda = 1.  
%  coef   - vector of the coefficients of the linear combination of the
%           GAMMA random variables. If coef is scalar, it is assumed
%           that all coefficients are equal. If empty, default value is
%           coef = 1.
%  niid   - scalar convolution coeficient niid, such that Z = Y +...+ Y is
%           sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%           coef(i) * X_i is independently and identically distributed
%           random variable. If empty, default value is niid = 1.   
% 
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Exponential_distribution.
%
% EXAMPLE 1:
%  % CF of the Exponential distribution with lambda = 5
%  lambda = 5;  
%  t = linspace(-50,50,501);
%  cf = cf_Exponential(t,lambda);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Exponential distribution with lambda = 5')
%
% EXAMPLE 2:
%  % CF of a linear combination of independent Exponential RVs
%  coef = 1./(((1:50) - 0.5)*pi).^2;
%  lambda = 5;  
%  t = linspace(-100,100,201);
%  cf = cf_Exponential(t,lambda,coef);
%  figure; plot(t,real(cf),t,imag(cf));grid on
%  title('CF of a linear combination of EXPONENTIAL RVs')
%
% EXAMPLE 3:
%  % PDF/CDF of the compound Binomial-Exponential distribution
%  n = 25;  
%  p = 0.3;
%  coef = 1./(((1:50) - 0.5)*pi).^2;
%  lambda = 5;  
%  cfX  = @(t) cf_Exponential(t,lambda,coef);
%  cf   = @(t) cfN_Binomial(t,n,p,cfX);
%  x    = linspace(0,5,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  options.N = 2^12;
%  options.SixSigmaRule = 15;
%  result = cf2DistGP(cf,x,prob,options)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-May-2017 18:11:50
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf_Exponential(t,lambda,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid = []; end
if nargin < 3, coef = []; end
if nargin < 2, lambda = []; end

%%
if isempty(lambda) 
    lambda = 1;
end

if isempty(coef) 
    coef = 1;
end

if isempty(niid)
    niid = 1;
end

%% Equal size of the parameters   
if ~isempty(coef) && isscalar(lambda) && isempty(niid)
    coef = sort(coef);
    m    = length(coef);
    [coef,idx] = unique(coef);
    lambda = lambda * diff([idx;m+1]);
end

[errorcode,coef,lambda] = distchck(2,coef(:)',lambda(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

% Special treatment for linear combinations with large number of RVs
szcoefs  = size(coef);
szcoefs  = szcoefs(1)*szcoefs(2);
szt      = size(t);
sz       = szt(1)*szt(2);
szcLimit = ceil(1e3 / (sz/2^16));
idc = 1:fix(szcoefs/szcLimit)+1;

%% Characteristic function
t        = t(:);
idx0     = 1;
cf       = 1;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux  = bsxfun(@times,t,coef(idx)./lambda(idx));
    cf   = cf .* prod(1./(1-1i*aux),2);
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