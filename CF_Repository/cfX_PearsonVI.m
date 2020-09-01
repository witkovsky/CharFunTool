function cf = cfX_PearsonVI(t,alpha,beta)
%cfX_PearsonVI(t,alpha,beta) 
% Characteristic function of the Pearson type VI distribution with
% parameters alpha > 0 (shape) and beta > 0 (scale), computed for real
% vector argument t, i.e.
%   cf(t) = cfX_PearsonVI(t,alpha,beta)
%         = (gamma(alpha+beta)./gamma(beta)) .* U(alpha,1-beta,-1i*t),
% where U(a,b,z) denotes the confluent hypergeometric function of the
% second kind. 
%
% NOTE: 
% This is an experimental algorithm which is not fully documented and that
% may provide unreliable results for a particular combination of parameters
% and / or may cause an unexpected failure.
% 
% WIKIPEDIA:
% https://en.wikipedia.org/wiki/Pearson_distribution
%
% SYNTAX:
%  cf = cfX_PearsonVI(t,alpha,beta)
%
% INPUTS:
%  t      -  vector of real values where the CF is evaluated, i.e. CF(t),
%  alpha  -  scalar shape parameter of the PearsonVI distribution, alpha > 0,
%  beta   -  scalar scale parameter of the PearsonVI distribution, beta > 0,
%
% OUTPUT:
% cf     -   calculated values of the characteristic function CF(t)
%
% EXAMPLE1: (CF of the PearsonVI distribution)
%  alpha = 3/2;
%  beta = 2/3;
%  t = linspace(-10,10,1001);
%  cf = cfX_PearsonVI(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the PearsonVI distribution with alpha = 3/2, beta = 2/3')
%
% EXAMPLE2 (PDF/CDF of the PearsonVI distribution)
%  alpha = 3/2;
%  beta = 2/3;
%  x = linspace(0,200,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^14;
%  cf = @(t) cfX_PearsonVI(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-PearsonVI distribution)
%  n = 25;  
%  p = 0.3;
%  alpha = 3/2;
%  beta = 2/3;
%  cfX = @(t) cfX_PearsonVI(t,alpha,beta);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,10000,101);
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
%cf = cfX_PearsonVI(t,alpha,beta);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end
if isempty(beta), beta = 1; end
if isempty(alpha), alpha = 1; end

%% Characteristic function of the Pearson VI distribution
szt  = size(t);
t    = t(:);
cf   = (gamma(alpha+beta)./gamma(beta)) .* U(alpha,1-beta,-1i*t); 
cf   = reshape(cf,szt);
cf(t==0) = 1;

end
%%
function f = U(a,b,z,tol)
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