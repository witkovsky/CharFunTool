function cf = cfS_Trapezoidal(t,lambda)
%cfS_Trapezoidal(t,lambda) evaluates the characteristic function cf(t) of
% the symmetric zero-mean Trapezoidal distribution on the interval
% (-1,1) with halfwidth of the constant pdf region given by the parameter
% lambda (offset, 0<=lambda<=1),i.e., it is  
% symmetric trapezoidal distribution with zero mean and variance V =
% (1+lambda^2)/6, i.e.     
%   cf(t) = cfS_Trapezoidal(t) 
%         = cfX_Rectangular(w*t))*cfX_Rectangular((1-w)*t);
%         = (sin(w*t)./(w*t)).*(sin((1-w)*t)./((1-w)*t)),
% where w  = (1+lambda)/2.
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Triangular_distribution
%
% SYNTAX
%  cf = cfS_Trapezoidal(t,lambda)
%
% EXAMPLE1 (CF of the symmetric Trapezoidal distribution with lambda = 0.5)
%  lambda = 0.5;
%  t = linspace(-50,50,501);
%  cf = cfS_Trapezoidal(t,lambda)
%  figure; plot(t,cf),grid
%  title('CF of the symmetric Trapezoidal distribution with lambda = 0.5')
%
% EXAMPLE2 (PDF/CDF of the symmetric Trapezoidal distribution, lambda = 0.5)
%  lambda = 0.5;
%  cf = @(t) cfS_Trapezoidal(t,lambda)
%  x = linspace(-1,1,101);
%  xRange = 2;
%  clear options
%  options.N = 2^10;
%  options.dt = 2*pi/xRange;
%  result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfS_Trapezoidal(t,lambda);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);

if nargin < 2, lambda = []; end

%%
if isempty(lambda)
    lambda = 0;
end

if (lambda<0 || lambda>1)
    error('Parametr d is out of [0,1]');
end

%% Characteristic function of the Exponential distribution
szt = size(t);
t   = t(:);

w  = (1+lambda)/2;
%cf = cfX_Rectangular(w*t).*cfX_Rectangular((1-w)*t);
cf = min(1,(sin(w*t)./(w*t)).*(sin((1-w)*t)./((1-w)*t)));
cf = reshape(cf,szt);
cf(t==0) = 1;

end