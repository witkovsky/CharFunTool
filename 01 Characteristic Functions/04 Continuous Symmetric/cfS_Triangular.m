function cf = cfS_Triangular(t)
%cfS_Triangular(t) evaluates the characteristic function cf(t) of
% the symmetric zero-mean Triangular distribution on the interval
% (-1,1)(symmetric triangular distribution with zero mean and variance
% 1/6), i.e.    
%   cf(t) = cfS_Triangular(t) = (2-2*cos(t))/t^2;
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Triangular_distribution
%
% SYNTAX
%  cf = cfS_Triangular(t)
%
% EXAMPLE1 (CF of the symmetric Triangular distribution on (-1,1))
%  t = linspace(-50,50,501);
%  cf = cfS_Triangular(t);
%  figure; plot(t,cf),grid
%  title('CF of the the symmetric Triangular distribution on (-1,1)')
%
% EXAMPLE2 (PDF/CDF of the the symmetric Triangular distribution on (-1,1))
%  cf = @(t) cfS_Triangular(t);
%  x = linspace(-2,2,101);
%  xRange = 2;
%  clear options
%  options.N = 2^5;
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
%cf = cfS_Triangular(t);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 1);

%% Characteristic function of the Exponential distribution
szt = size(t);
t   = t(:);

cf = min(1,(2-2*cos(t))./t.^2);
cf = reshape(cf,szt);
cf(t==0) = 1;

end