function cf = cf_Beta(t,alpha,beta)
%cf_Beta(t,alpha,beta) evaluates the characteristic function cf(t) of
% the Beta distribution with the shape parameters alpha > 0 and beta >0,
% defined on the interval (0,1), i.e. beta distribution with the mean Mean =
% alpha / (alpha + beta) and the Variance = (alpha*beta) /
% ((alpha+beta)^2*(alpha+beta+1)).  Then, the standard deviation is given
% by STD = sqrt(Variance), i.e.  
%   cf(t) = cf_Beta(t,alpha,beta) = 1F1(alpha ;alpha +beta ; i*t),
% where 1F1(.;.;.) is the Confluent hypergeometric function.
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Beta_distribution
%
% SYNTAX
%  cf = cf_Beta(t,alpha,beta)
%
% EXAMPLE1 (CF of the Beta distribution with alpha = 1/2, beta = 3/2)
%  alpha = 1/2;
%  beta = 3/2;
%  t = linspace(-100,100,501);
%  cf = cf_Beta(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the the Beta distribution')
%
% EXAMPLE2 (PDF/CDF of the the Beta with alpha = 1/2, beta = 3/2)
%  alpha = 1/2;
%  beta = 3/2;
%  cf = @(t) cf_Beta(t,alpha,beta);
%  x = linspace(0,1,101);
%  xRange = 1;
%  clear options
%  options.dt = 2*pi/xRange;
%  result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Feb-2017 14:58:48

%% ALGORITHM
%cf = cf_Beta(t,alpha,beta);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end
if isempty(alpha), alpha = 1; end
if isempty(beta), beta = 1; end


%% Characteristic function of the Beta distribution
szt = size(t);
t   = t(:);

cf = min(1,hypergeom1F1(alpha,alpha+beta,1i*t));
cf = reshape(cf,szt);
cf(t==0) = 1;

end