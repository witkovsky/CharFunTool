function cf = cfX_Weibull(t,alpha,beta,tol)
%cfX_Weibull(t,alpha,beta) 
%  Characteristic function cf(t) of the Weibull distribution with the
%  parameters alpha(scale parameter, elsewhere denoted also as lambda,
%  alpha > 0) and beta (shape parameter, elsewhere denoted as k, beta > 0),
%  for real (vector) argument t, i.e. 
%   cf(t) = cfX_Weibull(t,alpha,beta);
%
%  NOTE: 
%  This is an experimental algorithm which is not fully documented and that
%  may provide unreliable results for a particular combination of
%  parameters and / or may cause an unexpected failure.
%
%  WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Weibull_distribution
%
%  If W ~ Weibull(alpha,beta) then W = alpha * X^(1/beta), where X ~
%  Exponential(1). In particular, if alpha=1 and beta=1 W ~ Exponential(1).
%  From that, the integral definition of the characteristic function of the
%  Weilbull random variable W, by using PDF of the Exponential(1)
%  distribution, is given by
%   cf(t) = int_0^Inf exp(1i*t * alpha*x^(1/beta)) * exp(-x) dx
%
%  Alternatively, the characteristic function of the Weibull distribution
%  can be expressed by a special Fox H-function, see Duby (2017). In
%  particular, the CF can be expressed by the Mellin-Barnes integral
%  defining the Fox H-function. The standard (known) represetation of the
%  CF by the series expansion based on the moments is numerically stable
%  for abs(t*alpha) < 1 and beta >= 1, only.
%
%  In summary,
%   cf(t) = int_0^Inf exp(1i*t * alpha*x^(1/beta)) * exp(-x) dx
%         = 1/(2i*pi)* int_L gamma(u)*gamma(1-u/beta)*(-1i*t*alpha)^(-u) du
%         = H^{1,1}_{1,1}(-1i*t*alpha | (0,1/beta);(0,1))
%         = FoxWrightPsi(1,1/beta,[],[],-1i*alpha*t)
%         = sum_{n=0}^inf (-1i*alpha*t)^n * gamma(1 + n/beta) / gamma(1+n)
%  For  beta < 1 we can use the asymptotic expansion based on the
%  propeties of the H-functions, see Mathai et al. (2009) and Duby (2017),
%   cf(t) = sum_{n=0}^inf (-1)^n * beta * (1/(alpha*1i*t)^(beta*(1+n)) ...
%           * gamma(beta*(1+n)) / gamma(1+n)
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
% EXAMPLE4 (CDF/PDF of the  Weibull distribution with alpha=1, beta=0.75)
%  alpha = 1;
%  beta  = 0.75;
%  x    = linspace(0,10,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.xMax = 20;
%  options.N = 2^10;
%  cf = @(t) cfX_Weibull(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE5 (PDF/CDF of the compound Poisson-Weibull distribution)
%  alpha = 1;
%  beta = 0.75;
%  lambda = 10;
%  cfX = @(t)cfX_Weibull(t,alpha,beta);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,40,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] DUBY T. (2017). Characteristic function of Weibull distribution.
%     Work in Progress.
% [2] MATHAI, A.M., SAXENA, R.K., HAUBOLD, H.J. (2009). The H-function:
%     theory and applications. Springer Science & Business Media.
% [3] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [4] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 13:25:21
%
% Revision history:
% Ver.: 29-Oct-2017 11:25:57

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
szt = size(t);
t   = t(:);
cf  = zeros(size(t));

if beta > 1
    T = 25*(1+beta^0.75);
    % CF by using direct integration with exponential PDF
    id = abs(t*alpha) < T;
    if any(id)
        cf(id) = cfWintegral(t(id),alpha,beta,tol);
    end
    % CF by using linear approximation of log(cf) for large t
    id = abs(t*alpha) >= T;
    if any(id)
        cf(id) = cfFitFun(t(id),alpha,beta,T);
    end
elseif beta == 1
    % CF by using the exact CF of the exponential distribution
    cf = alpha ./ (alpha - 1i*t);
else
    % CF by using the inegral representation of the H-function
    T  = 1.1;
    id = abs(t)>0 & abs(t*alpha) <  T;
    if any(id)
        cf(id) = cfHintegral(t(id),alpha,beta,tol);
    end
    % CF by using the asymptotic expansion of the H-function
    id = abs(t*alpha) >= T;
    if any(id)
        cf(id) = cfHasymptotic(t(id),alpha,beta,100);
    end
end

cf = reshape(cf,szt);
cf(t==0) = 1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function cfWintegral
function [cf,f] = cfWintegral(t,alpha,beta,reltol)
%  cfWintegral Computes the Weibull characteristic function by using the
%  direct integration. If W ~ Weibull(alpha,beta) then W = alpha *
%  X^(1/beta), where X ~ Exponential(1). From that we get the integral
%  representation of the characteristic function of the Weilbull random
%  variable W, using PDF of the Exponential(1) distribution, given by
%  cf(t) = int_0^Inf exp(1i*t*alpha*x^(1/beta)) * exp(-x) dx

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 18:01:45

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, reltol = []; end

if isempty(reltol)
    reltol = 1e-6;
end

%% ALGORITHM
sz = size(t);

% f = @(x) integrandWfun(x,t,alpha,beta);
% cf = integral(@(x) integrandWfun(x,t,alpha,beta),0,+Inf, ...
%     'ArrayValued',true,'RelTol',reltol);

f = @(x) bsxfun(@times,integrandW2fun((x./(1-x)).^2,t,alpha,beta), ...
    2*x./(1-x).^3);
cf = integral(@(x) bsxfun(@times,integrandWfun((x./(1-x)).^2, ...
    t,alpha,beta), 2*x./(1-x).^3),0,1,'ArrayValued',true,'RelTol',reltol);

cf = reshape(cf,sz);

end
%% Function integrandWfun
function f = integrandWfun(x,t,alpha,beta)
% integrandWfun Integrand of the integral representation of characteristic
% function  cf(t) = int_0^Inf exp(1i*alpha*t*x^(1/beta)) * exp(-x) dx

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 18:01:45
%%
x = x(:);
t = t(:)';

f =  exp(1i*alpha*x.^(1/beta)*t) .* exp(-x*ones(size(t)));
end

%% Function cfHintegral
function cf = cfHintegral(t,alpha,beta,reltol)
% cfHintegral Computes the Weibull characteristic function by using the
% integral representation of the H-function, see Duby (2017).

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 18:01:45

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, reltol = []; end

if isempty(reltol)
    reltol = 1e-6;
end

%% ALGORITHM
sz = size(t);
cf = integral(@(x) integrandHfun(x,t,alpha,beta),-Inf,+Inf, ...
    'ArrayValued',true,'RelTol',reltol)/(2*pi);
cf = reshape(cf,sz);

end
%% Function integrandHfun
function f = integrandHfun(v,t,alpha,beta)
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
%% Function cfHasymptotic
function cf = cfHasymptotic(t,alpha,beta,N)
% cfHasymptotic Computes the Weibull characteristic function by using the
% the aymptotic expansion (for large values t) of the Fox's H-function for
% 0 < beta < 1, see Duby (2017).

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 18:01:45

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

% Optimum stopping rule for the (divergent) asymptotic expansion
cf = bsxfun(@plus,GammaLog((1+n)*beta)-GammaLog(1+n),-beta*(1+n)*log(z));
[~,id] = min(max(-325,real(cf)));

% CF by using the asymptotic expansion of the H-function
cf = cumsum(bsxfun(@times,beta*(-1).^n,exp(cf)));
cf = diag(cf(id,:));

cf(z==0) = 1;
cf = reshape(cf,sz);

end
%% Function cfFitFun
function f = cfFitFun(t,alpha,beta,T)
% cfFitFun estimates the asymptotic behaviour of CF for large abs(t) by
% using linear regression fitted for log(cf)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Oct-2017 11:14:14

%%
tt = linspace(T/100,T,7)';
cf = cfWintegral(tt,alpha,beta);
W  = diag(tt/norm(tt));
X  = [tt.^0 tt];
b  = (X'*W*X)\(X'*W*log(cf));

id = t >= 0;
if any(id)
    f(id) = exp([t(id).^0 t(id)] * b);
end

id = t < 0;
if any(id)
    f(id) = conj(exp([t(id).^0 -t(id)] * b));
end

end