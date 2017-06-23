function cf = cf_Stable(t,alpha,mu,sigma)
%cf_Stable(t,alpha,mu,sigma) evaluates the characteristic function
% cf(t) of the symmetric alpha stable (location–scale symmetric stable
% Paretian) distribution with the stability parameters 0 < alpha <= 2,
% location parameter mu in (??, ?), and the scale parameter sigma > 0, 
% i.e.
%   cf(t) = cf_Stable(t,alpha,mu,sigma) = 
%         = exp( 1i*mu*t - |sigma*t|^alpha).
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Stable_distribution
%
% SYNTAX
%  cf = cf_Stable(t,alpha,mu,sigma)
%
% EXAMPLE1 (CF of the symmetric Alpha Stable distribution with alpha=1)
%  alpha = 1;
%  t = linspace(-10,10,501);
%  cf = cf_Stable(t,alpha);
%  figure; plot(t,cf),grid
%  title('CF of the symmetric Alpha Stable distribution with alpha=1')
%
% EXAMPLE2 (PDF/CDF the symmetric Alpha Stable distribution with alpha=1)
%  alpha = 1;
%  x = linspace(-5,5,501);
%  cf = @(t) cf_Stable(t,alpha);
%  clear options
%  options.N = 2^12;
%  options.SixSigmaRule = 3;
%  result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Feb-2017 14:35:54

%% ALGORITHM
%cf = cf_Stable(t,alpha,mu,sigma);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, alpha = []; end
if nargin < 3, mu = []; end
if nargin < 4, sigma = []; end

if isempty(alpha), alpha = 2; end
if isempty(mu), mu = 0; end
if isempty(sigma), sigma = 1; end

%% Characteristic function of the Exponential distribution
szt = size(t);
t   = t(:);

cf  = exp(1i*mu*t - abs(sigma*t).^alpha);
cf = reshape(cf,szt);
cf(t==0) = 1;

end