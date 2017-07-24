function cf = cfX_Weibull(t,alpha,beta,tol)
%cfX_Weibull(t,alpha,beta) Computes the characteristic function cf(t)
%  of the Weibull distribution with the parameters alpha(scale parameter,
%  elsewhere denoted also as lambda, alpha > 0) and beta (shape parameter,
%  elsewhere denoted as k, beta > 0), for real (vector) argument t, i.e.
%  cf(t) = cfX_Weibull(t,alpha,beta);
%
%  The characteristic function of the Weibull distribution can be expressed
%  by a special Fox H-function, see Duby (2017). In particular, the CF can
%  be expressed by the Mellin-Barnes integral defining the Fox H-function.
%  The standard (known) represetation of the CF by the series expansion
%  based on the moments is numerically stable for abs(t*alpha) < 1 and beta
%  >= 1, only.
%
%  In particular,
%   cf(t) = 1/(2i*pi)* int_L gamma(u)*gamma(1-u/beta)*(-1i*t*alpha)^(-u) du
%         = H^{1,1}_{1,1}(-1i*t*alpha | (0,1/beta);(0,1))
%         = FoxWrightPsi(1,1/beta,[],[],-1i*alpha*t)
%         = sum_{n=0}^inf (-1i*alpha*t)^n * gamma(1 + n/beta) / gamma(1+n)
%  For  beta < 1 we can use the asymptotic expansion based on the
%  propeties of the H-functions, see Mathai et al. (2009) and Duby (2017),
%   cf(t) = sum_{n=0}^inf (-1)^n * beta * (1/(alpha*1i*t)^(beta*(1+n)) ...
%           * gamma(beta*(1+n)) / gamma(1+n)
%  
%  For more details see also WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Weibull_distribution
%
% SYNTAX:
%  cf = cfX_Weibull(t,alpha,beta)
%  cf = cfX_Weibull(t,alpha,beta,tol)
%
% EXAMPLE1 (CF of the Weibull distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta  = 1;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE2 (CF of the Weibull distribution with alpha=1, beta<1)
%  alpha = 1;
%  beta  = 0.5;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE3 (CF of the Weibull distribution with alpha=1, beta>1)
%  alpha = 1;
%  beta  = 5.5;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE4 (CDF/PDF of the  Weibull distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta  = 1;
%  x    = linspace(0,7,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^10;
%  cf = @(t) cfX_Weibull(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE5 (PDF/CDF of the compound Poisson-Weibull distribution)
%  alpha = 1;
%  beta = 1;
%  lambda = 10;
%  cfX = @(t)cfX_Weibull(t,alpha,beta);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,35,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] DUBY T. (2017). Characteristic function of Weibull distribution.
%     Work in Progress. 
% [2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models.  Working Paper.
%     Journal of Statistical Software.
% [3] MATHAI, A.M., SAXENA, R.K., HAUBOLD, H.J. (2009). The H-function:
%     theory and applications. Springer Science & Business Media. 
% [4] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [5] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  
% [6] WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. arXiv preprint
%     arXiv:1701.08299.     

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Jul-2017 11:09:09

%% ALGORITHM
%cf = cfX_Weibull(t,alpha,beta,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cfX_Weibull:TooFewInputs'));
end
if nargin < 2, alpha = []; end
if nargin < 3, beta = []; end
if nargin < 4, tol = []; end

if isempty(alpha), alpha = 1; end
if isempty(beta), beta = 1; end
if isempty(tol), tol = 1e-6; end

alpha(alpha <= 0) = NaN;
beta(beta <= 0) = NaN;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_Weibull(t,alpha,beta)
szt      = size(t);
t        = t(:);
cf       = zeros(size(t));

if beta > 1
    id     = abs(t*alpha) < 0.9;
    cf(id) = FoxWrightPsi(1,1/beta,[],[],1i*alpha*t(id));
    id     = abs(t*alpha) >= 0.9;
    cf(id) = cf_HIntegral(t(id),alpha,beta,tol);
elseif beta == 1
    cf     = alpha ./ (alpha - 1i*t);
else
    id     = abs(t)>0 & abs(t*alpha) < 1.1;
    cf(id) = cf_HIntegral(t(id),alpha,beta,tol);
    id     = abs(t*alpha) >= 1.1;
    cf(id) = cf_HAsymptotic(t(id),alpha,beta,100);
end

cf       = reshape(cf,szt);
cf(t==0) = 1;

end
%% Function cf_HAsymptotic
function cf = cf_HAsymptotic(t,alpha,beta,N)
% cf_HAsymptotic Computes the Weibull characteristic function by using the
% the aymptotic expansion (for large values t) of the Fox's H-function for
% 0< beta < 1, see Duby (2017).

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Jul-2017 12:58:35

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, N = []; end

if isempty(N)
    N = 500;
end

%% ALGORITHM
sz = size(t);
t  = t(:);
z  = -1i*t'*alpha;
n  = (0:N)';

cf = bsxfun(@plus,GammaLog((1+n)*beta)-GammaLog(1+n),-beta*(1+n)*log(z));
cf = sum(bsxfun(@times,beta*(-1).^n,exp(cf)));

cf(z==0) = 1;
cf  = reshape(cf,sz);

end
%% Function cf_HIntegral
function cf = cf_HIntegral(t,alpha,beta,reltol)
% cf_HIntegral Computes the Weibull characteristic function by using the
% integral representation of the H-function, see Duby (2017).  

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Jul-2017 12:58:35

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, reltol = []; end

if isempty(reltol)
    reltol = 1e-6;
end

%% ALGORITHM
sz = size(t);
cf = integral(@(x) integrandFun(x,t,alpha,beta),-Inf,+Inf, ...
    'ArrayValued',true,'RelTol',reltol)/(2*pi);
cf = reshape(cf,sz);

end
%%
function f = integrandFun(v,t,alpha,beta)
% integrandFun Integrand of the integral representation of the H-function
% based, see Duby (2017).   

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Jul-2017 12:58:35    
%%
g = beta/2;
v = v(:);
t = t(:)';
f = exp((GammaLog(g+1i*v)+GammaLog(1-(g+1i*v)/beta))*ones(size(t)) ...
     -(g+1i*v)*log(-1i*t*alpha));
end