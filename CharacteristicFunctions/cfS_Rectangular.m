function cf = cfS_Rectangular(t)
%cfS_Rectangular(t) evaluates the characteristic function cf(t) of
% the symmetric zero-mean Rectangular distribution on the interval
% (-1,1)(rectangular distribution with zero mean and variance 1/3), i.e.   
%   cf(t) = cfS_Rectangular(t) = sinc(t) = sin(t)/t;
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
%
% SYNTAX
%  cf = cfS_Rectangular(t)
%
% EXAMPLE1 (CF of the Rectangular distribution on (-1,1))
%  t = linspace(-50,50,501);
%  cf = cfS_Rectangular(t);
%  figure; plot(t,cf),grid
%  title('CF of the Rectangular distribution on (-1,1)')
%
% EXAMPLE2 (PDF/CDF of the Rectangular distribution on (-1,1))
%  cf = @(t) cfS_Rectangular(t);
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
%cf = cfS_Rectangular(t);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 1);

%% Characteristic function of the Exponential distribution
szt = size(t);
t   = t(:);

cf = min(1,sin(t)./t);
cf = reshape(cf,szt);
cf(t==0) = 1;

end