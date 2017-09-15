function cf = cf_vonMises(t,mu,kappa,coef,niid)
%% cf_vonMises(t) 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent VON MISES random variables. 
%
%  The VON MISES distribution is circular distribution on the interval of
%  length 2*pi, here we consider (-pi,pi), equivalent of the normal
%  distribution with the real parameter mu and rate parameter kappa > 0 (mu
%  and 1/kappa are analogous to mu and sigma^2, the mean and variance in
%  the normal distribution), on a whole circle, i.e. the interval of angles
%  (-pi,pi). 
%
%  cf_vonMises evaluates the characteristic function cf(t) of  Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ vonMises(mu_i,kappa_i) are
%  inedependent RVs, with the locarion parameters mu_i in Real and the rate
%  parameters kappa_i > 0, for i = 1,...,N.
%
%  The characteristic function of the vonMises(mu,kappa) distribution is 
%   cf(t) = cf_vonMises(t,mu,kappa) 
%         = besseli(t,kappa)/besseli(0,kappa) .* exp(1i*t*mu).
%
% SYNTAX
%  cf = cf_vonMises(t,mu,kappa,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of location parameters, mu in (-pi,pi). If empty, default
%          value is mu = 0. 
%  kappa - vector of rate (1/scale) parameters, kappa_i > 0. If empty,
%          default value is kappa = 0, i.e. uniform distribution on
%          (-pi,pi).  
%  coef  - vector of the coefficients of the linear combination of the
%          IGamma random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y +...+ Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef_i * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1. 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Von_Mises_distribution.
%
% EXAMPLE 1:
% % CF of the weighted linear combinantion of the von Mises RVs
%   mu    = [0 0 -pi/2 pi/2 0];
%   kappa = [1 2 3 4 5];
%   coef  = [1 2 3 4 5]/15;
%   t     = linspace(-20,20,201);
%   cf    = cf_vonMises(t,mu,kappa,coef);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of the weighted linear combinantion of the von Mises RVs')
%
% EXAMPLE 2:
% % CDR/PDF of the weighted linear combinantion of the von Mises RVs
%   mu    = [0 0 -pi/2 pi/2 0];
%   kappa = [1 2 3 4 5];
%   coef  = [1 2 3 4 5]/15;
%   t     = linspace(-20,20,201);
%   cf    = @(t) cf_vonMises(t,mu,kappa,coef);
%   x     = linspace(-pi,pi,201);
%   prob  = [0.9 0.95 0.99];
%   clear options
%   options.xMin = -pi;
%   options.xMax = pi;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%   angle  = result.x;
%   radius = result.pdf;
%   figure; polarplot(angle,radius);
%   ax = gca; ax.ThetaAxisUnits = 'radians';
%
% EXAMPLE 3:
% % CF of the mixture of the von Mises distribution on (-pi,pi)
%   mu1 = 0; kappa1 = 5;
%   mu2 = 1; kappa2 = 15;
%   cf = @(t) 0.25*cf_vonMises(t,mu1,kappa1) + 0.75*cf_vonMises(t,mu2,kappa2);
%   clear options
%   options.xMin = -pi;
%   options.xMax = pi;
%   result = cf2DistGP(cf,[],[],options)
%   angle  = result.x;
%   radius = result.pdf;
%   figure; polarplot(angle,radius);
%   ax = gca; ax.ThetaAxisUnits = 'radians';
%
% EXAMPLE 4:
% % PDF/CDF of the mixture of the von Mises distribution on (0,2*pi)
%   mu1 = 0;  kappa1 = 5;
%   mu2 = 1;  kappa2 = 15;
%   mu3 = pi; kappa3 = 10;
%   cf  = @(t) 0.25*cf_vonMises(t,mu1,kappa1) + ...
%         0.25*cf_vonMises(t,mu2,kappa2) + 0.5*cf_vonMises(t,mu3,kappa3);
%  clear options
%  options.isCircular   = true;
%  options.correctedCDF = true;
%  options.xMin = 0;
%  options.xMax = 2*pi;
%  x = linspace(0,2*pi);
%  result = cf2DistGP(cf,x,[],options)
%  angle  = result.x;
%  radius = result.pdf;
%  figure; polarplot(angle,radius);
%  ax = gca; ax.ThetaAxisUnits = 'radians';

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jun-2017 18:25:56

%% ALGORITHM
%cf = cf_vonMises(t,mu,kappa,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, kappa = []; end
if nargin < 2, mu = []; end

%%
if isempty(mu), mu = 0; end
if isempty(kappa), kappa = 0; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

[errorcode,coef,mu,kappa] = distchck(3,coef(:)',mu(:)',kappa(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function 
szt = size(t);
t   = t(:);

if length(coef)==1
    cf  = (besseli(abs(t*coef),kappa,1) ...
        ./ besseli(0,kappa,1)) .* exp(1i*t*mu*coef);
else
    cf  = prod((besseli(abs(t*coef),ones(size(t))*kappa,1) ...
        ./ besseli(0,ones(size(t))*kappa,1)) .* exp(1i*t*(mu.*coef)),2);
end

cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end