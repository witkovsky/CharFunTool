function cf = cfX_Pareto(t,alpha,sigma,type,tol)
%cfX_Pareto(t,alpha,sigma,type,tol)
% Characteristic function of the Pareto distribution with parameters alpha
% > 0 (shape), sigma > 0 (scale, sometimes denoted by x_m - which is the
% minimum value of the distribution support, i.e. x >= x_m), and the
% parameter type ('typeI' or 'typeII'), computed for real vector argument
% t, i.e.  
%   cf(t) = cfX_Pareto(t,alpha,sigma,type);
%
% NOTE: 
% This is an experimental algorithm which is not fully documented and that
% may provide unreliable results for a particular combination of parameters
% and / or may cause an unexpected failure.
%
% WIKIPEDIA:
% https://en.wikipedia.org/wiki/Pareto_distribution
%
% SYNTAX:
%  cf = cfX_Pareto(t,alpha)
%  cf = cfX_Pareto(t,alpha,sigma,type)
%
% INPUTS:
%  t      -  vector of real values where the CF is evaluated, i.e. CF(t),
%  alpha  -  scalar shape parameter of the Pareto distribution, alpha > 0,
%  sigma  -  scalar scale parameter of the Pareto distribution, sigma > 0,
%  type   -  type of the Pareto distribution: 'typeI' or 'typeII',
%  tol    -  relative tolerance of the integral function (tol = 1e-6).
%
% OUTPUT:
% cf     -   calculated values of the characteristic function CF(t)
%
% FURTHER DETAILS:
% The Pareto Type I distribution (also known as European Pareto) is
% characterized by a scale parameter sigma and a shape parameter alpha,
% which is known as the tail index. The pdf of X is defined by pdf_X(x) =
% alpha * sigma^alpha / x^(alpha+1), for x >= sigma, otherwise it is zero.
% When this distribution is used to model the distribution of wealth, then
% the parameter alpha is called the Pareto index.
%
% The Type II Pareto distribution (also known as Lomax or American Pareto
% distribution, type = 'typeII' or type = 'Amer') is shifted to zero.
% The default type is 'TypeI' (type = 'typeI' or type = 'Eur').
% The CFs of the Pareto type distributions are defined by
% cf_{typeI}(t)  = alpha * exp(1i*sigma*t) * U(1,1-alpha,-1i*sigma*t)
% cf_{typeII}(t) = alpha * U(1,1-alpha,-1i*sigma*t), where U(a,b,z) denotes
% the confluent hypergeometric function of the second kind. 
%
% A special case of Pareto TypeI distribution is the Pareto 1 distribution
% with sigma = 1, cf_{1}(t) = alpha * exp(1i*t) * U(1,1-alpha,-1i*t).
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Pareto_distribution
%
% EXAMPLE1 (CF of the Pareto 1 distribution)
%  alpha = 3/2;  
%  t = linspace(-10,10,1001);
%  cf = cfX_Pareto(t,alpha);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Pareto 1 distribution with alpha = 3/2')
%
% EXAMPLE2 (CF of the Pareto Type I/EUR distribution)
%  alpha = 3/2;  
%  sigma  = 2;
%  type = 'eur';
%  t = linspace(-10,10,1001);
%  cf = cfX_Pareto(t,alpha,sigma,type);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Pareto EUR distribution with alpha = 3/2, sigma = 2')
%
% EXAMPLE3 (CF of the Pareto Type II/AMER distribution)
%  alpha = 3/2;  
%  sigma  = 2;
%  type = 'amer';
%  t = linspace(-10,10,1001);
%  cf = cfX_Pareto(t,alpha,sigma,type);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Pareto AMER distribution with alpha = 3/2, sigma = 2')
%
% EXAMPLE4 (PDF/CDF of the compound Poisson-Pareto distribution)
%  alpha = 3/2;  
%  sigma  = 2;
%  type = 'eur';
%  lambda = 10;
%  cfX = @(t) cfX_Pareto(t,alpha,sigma,type);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,500,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [2] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 13:25:21
%
% Revision history:
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfX_Pareto(t,alpha,sigma,type,tol);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, tol = []; end
if nargin < 4, type = []; end
if nargin < 3, sigma = []; end
if nargin < 2, alpha = []; end
if isempty(tol), tol = 1e-6; end
if isempty(type), type = 'TypeI'; end
if isempty(sigma), sigma = 1; end
if isempty(alpha), alpha = 1; end

%% Characteristic function of the Pareto distribution
szt  = size(t);
t    = t(:);
switch lower(type)
    case lower({'eur','typei'})
        cf  = alpha .* exp(1i*sigma.*t) .* auxFun(1,1-alpha,-1i*sigma.*t,tol);         
    case lower({'amer','lomax','typeii'})
        cf  = alpha .* auxFun(1,1-alpha,-1i*sigma.*t,tol);
    otherwise
        cf  = alpha .* exp(1i*sigma.*t) .* auxFun(1,1-alpha,-1i*sigma.*t,tol); 
end
cf = reshape(cf,szt);
cf(t==0) = 1;

end
%%
function f = auxFun(a,b,z,tol)
% Computes the required auxiliary function.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 05-Apr-2016 17:34:57

%%
narginchk(2, 4);
if nargin < 4, tol = 1e-6; end
if nargin < 3, z = []; end

if any(real(z)~=0)
     error('Error. \nInput z must be purely imaginary complex vector')
end
%
sz = size(z);
z  = z(:);
f  = integral(@(x) bsxfun(@times,integrand(a,b,z,(x/(1-x))^2), ...
    2*x/(1-x)^3),0,1,'ArrayValued',true,'RelTol',tol);
f  = reshape(f,sz);

end

%%
function f = integrand(a,b,z,x)
% Computes the required integrand

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 05-Apr-2016 17:34:57

%%
z  = imag(z(:));
nz = length(z);
oz = ones(nz,1);
x  = x(:)';
nx = length(x);
ox = ones(1,nx);
c  = - gammaln(a);
f  = zeros(nz,nx);
id = (z==0);
if any(id)
    f(id,:) = exp(c + (a-1).*log(oz(id)*x) + (b-a-1).*log(1+oz(id)*x));
end
id = (abs(z)>=1);
if any(id)
    zi = -1i./z(id);
    f(id,:) = exp(c + log(zi*ox) + (a-1)*log(zi*x) + ...
        (b-a-1)*log(1+zi*x) - oz(id)*x);
end
id = (z>0 & z<1);
if any(id)
    f(id,:) = exp(c + log(-1i) + (a-1)*log(-1i*oz(id)*x) + ...
        (b-a-1)*log(1-1i*oz(id)*x) - (z(id))*x);
end
id = (z<0 & z>-1);
if any(id)
    f(id,:) = exp(c + log(1i) + (a-1)*log(1i*oz(id)*x) + ...
        (b-a-1)*log(1+1i*oz(id)*x) + (z(id))*x);
end

end