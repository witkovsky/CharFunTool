function cf = cf_FoldedDistribution(t,cfX,coef,n,tol)
% cf_FoldedDistribution Characteristic function of a linear combination
%  (convolution) of iid random variables X with given (common) FOLDED
%  DISTRIBUTION which is specified by its characteristic function. 
% 
%  The CF of a FOLDED DISTRIBUTION, say cf, is specified form the given
%  characteristic function cfX (from its real part) by convolving it with
%  the Fourier tranfsorm of the Heaviside function. 
%
%  Hence, CF of the  folded distribution is a characteristic function of a
%  distribution specified by cfX, which is folded about 0 to the positive
%  part of the real axis.   
%
% SYNTAX:
%     cf = cf_FoldedDistribution(t,cfX,coef,n)
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
%     tol      - relative tolerance for numerical integration. If empty,
%                the deafult value is tol = 1e-6.
%
% OUTPUT:
%     cf       - The characteristic function of a linear combination of
%                iid RVs with characteristic function cf, evaluated at t.
%
% EXAMPLE 1 (CF of the Half-Normal Distribution)
%  cfX = @(t) cf_Normal(t);
%  t   = linspace(-10,10,101);
%  cf  = cf_FoldedDistribution(t,cfX);
%  figure;plot(t,real(cf),t,imag(cf));grid
%  title('CF of the Half-Normal Distribution')
% 
% EXAMPLE 2 (CF of the Folded Normal Distribution)
%  mu = 2;
%  sigma = 1;
%  cfX = @(t) cf_Normal(t,mu,sigma);
%  t   = linspace(-10,10,101);
%  cf  = cf_FoldedDistribution(t,cfX);
%  figure;plot(t,real(cf),t,imag(cf));grid
%  title('CF of the Folded Normal Distribution')
%
% EXAMPLE 3 (PDF/CDF/QF of the Folded Normal Distribution)
%  mu = 2;
%  sigma = 1;
%  cfX = @(t) cf_Normal(t,mu,sigma);
%  cf  = @(t) cf_FoldedDistribution(t,cfX);
%  prob = [0.9 0.95 0.99];
%  clear options;
%  options.xMin = 0;
%  result = cf2DistGP(cf,[],prob,options)
%
% EXAMPLE 4 (PDF/CDF/QF of the Folded Distribution of a linear combination)
%  mu    = [2 2 3 4 5];
%  sigma = [1 0.5 2 3 4];
%  coef  = [1 -1 1 2 -3];
%  cfX = @(t) cf_Normal(t,mu,sigma,coef);
%  cf  = @(t) cf_FoldedDistribution(t,cfX);
%  prob = [0.9 0.95 0.99];
%  clear options;
%  options.xMin = 0;
%  result = cf2DistGP(cf,[],prob,options)

% (c) 2021, Viktor Witkovsky (witkovsky@savba.sk)
% Based on the idea by Tomy Duby, OAA Computing Ltd, 27-Mar-2021.
% Ver.: 11-Apr-2021 20:14:14

%% CHECK THE INPUT PARAMETERS
narginchk(2, 4);
if nargin < 3, coef = []; end
if nargin < 4, n = []; end
if nargin < 5, tol = []; end

if  isempty(coef)
    coef = 1;
end

if  isempty(tol)
    tol = 1e-6;
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

%% ALGORITHM
szt  = size(t);
t    = t(:);
cf   = 1;
m    = length(coef);
for j = 1:m
    cf = cf .* cfConv(coef(j)*t,cfX,tol).^nums(j);
end
cf       = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(n)
    if isscalar(n)
        cf = cf .^ n;
    else
        error('n should be a scalar (positive integer) value');
    end
end
end
%% Function cfConv
function cf = cfConv(t,cfX,tol)
% cfConv is an auxiliary function that calculates the convolution of a real
% part of  the characteristic function cfX with the Fourier tranfsorm of
% the Heaviside function. 
% 
% Based on the idea created by Tomy Duby, OAA Computing Ltd, 27-Mar-21.

% (c) 2021, Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 11-Apr-2021 20:14:14

%% ALGORITHM

cfReal   = @(omega) real(cfX(omega));
funConv  = @(u,omega) (-cfReal(-omega-u)+cfReal(omega-u))./u/pi;
fun      = @(x,omega) funConv((x./(1-x)).^2,omega) .* (2*x./(1-x).^3);
cfImag   = 1i*integral(@(x) fun(x,t),0,1,'ArrayValued',true,'RelTol',tol);
cf       = cfReal(t) + cfImag;

end
