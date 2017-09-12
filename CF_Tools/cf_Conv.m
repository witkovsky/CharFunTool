function cf = cf_Conv(t,cfX,coef,n)
% cf_Conv Characteristic function of a linear combination (convolution)
%  of iid random variables X with given (common) characteristic function
%  cfX(t) and the coefficients coef, and such that Y = coef(1)*X + ... +
%  coef(N)*X, i.e. 
%     cf = cfX(coef(1)*t) * ... * cfX(coef(N)*t). 
%
%  If or all coef == 1, then cf = cfX(t)^N, . Moreover, with using the
%  optional parameter n, the algorithm evaluates CF of the convolution of n
%  independent copies of the random variable Y, i.e. Z = Y + ... + Y,
%  where Y = coef(1)*X + ... + coef(N)*X, i.e.
%    cf = (cfX(coef(1)*t) * ... * cfX(coef(N)*t))^n.
%
% SYNTAX:
%     cf = cf_Conv(t,cfX,coef,n)
%
% INPUT:
%     t        - vector (or array) of input values t where the cf_conv is
%                evaluated
%     cfX      - function handle to the given chracteristic function cfX(t)
%     coef     - vector of coeficients of the linear combination of the
%                iid random variables, such that Y = sum(coef(k) * X_k),
%                and cf = Prod_{k=1}^N cfX(coef(k)*t).
%     n        - optional power coeficient of additional convolution of
%                the combiend CF of Y. With using n (if empty, default
%                value is n = 1), cf = (Prod_{k=1}^N cfX(coef(k)*t)^n.
%
% OUTPUT:
%     cf       - The characteristic function of a linear combination of
%                iid RVs with characteristic function cf, evaluated at t.
%
% EXAMPLE:
% % CF of a linear combination of chi-square random variables
%     df = 1;
%     cfX = @(t) cfX_ChiSquared(t,df);
%     coef = 1./(1:100);
%     cf = @(t) cf_Conv(t,cfX,coef);
%     t = linspace(-10,10,501);
%     plot(t, real(cf(t)),t,imag(cf(t)));grid on
%     title('CF of a Linear Combination of iid Chi-Square RVs')
%     clear options;
%     options.xMin = 0;
%     result = cf2DistGP(cf,[],[],options)
%
% EXAMPLE:
% % CF of a linear combination of iid RVs sampled from the empirical
% % distribution function (atrtificially generated data)
%     n = 30;
%     p = [0.2 0.7 0.1];
%     data  = p(1) * chi2rnd(5,n,1) + p(2) * normrnd(10,1,n,1) ...
%             + p(3) * trnd(1,n,1);
%     cfE   = @(t) cfE_Empirical(t,data);
%     coef = 1./[ 1 2 3 4 5 6 7 8 9 10];
%     cf    = @(t) cf_Conv(t,cfE,coef);
%     t = linspace(-10,10,501);
%     plot(t, real(cf(t)),t,imag(cf(t)));grid on
%     title('CF of a Linear Combination of iid RVs')
%     result = cf2DistGP(cf)

% Copyright (c) 2017, Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 9-May-2017 10:22:48

%% CHECK THE INPUT PARAMETERS
narginchk(2, 4);
if nargin < 3, coef = []; end
if nargin < 4, n = []; end

if  isempty(coef)
    coef = 1;
end

%% Find the unique coefficients and their multiplicities
if ~isscalar(coef)
    coef = sort(coef);
    m     = length(coef);
    [coef,idx] = unique(coef);
    nums = diff([idx;m+1]);
else
    nums = 1;
end

[errorcode,coef,nums] = distchck(2,coef(:)',nums(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

% Special treatment for linear combinations with large number of RVs
szcoef  = size(coef);
szcoef  = szcoef(1)*szcoef(2);
szt      = size(t);
sz       = szt(1)*szt(2);
szcLimit = ceil(1e3 / (sz/2^16));
idc      = 1:fix(szcoef/szcLimit)+1;

%% ALGORITHM
t     = t(:);
idx0  = 1;
cf    = 1;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoef);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    aux  = bsxfun(@times,t,coef(idx));
    aux  = bsxfun(@power,cfX(aux),nums(idx));
    cf   = cf .* prod(aux,2);
end
cf = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(n)
    if isscalar(n)
        cf = cf .^ n;
    else
        error('n should be a scalar (positive integer) value');
    end
end
end