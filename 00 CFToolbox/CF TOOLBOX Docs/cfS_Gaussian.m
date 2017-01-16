%function cf = cfS_Gaussian(t)
%cfS_Gaussian(t) evaluates the characteristic function cf(t) of
% the symmetric zero-mean standard Gaussian distribution (i.e. the standard
% normal distribution with zero mean and unit variance: N(0,1)), i.e.   
%   cf(t) = cfS_Gaussian(t) = exp(-t^2/2)
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Normal_distribution
%
% SYNTAX
%  cf = cfS_Gaussian(t)
%
% EXAMPLE1 (CF of the Gaussian distribution N(0,1))
%  t = linspace(-5,5,501);
%  cf = cfS_Gaussian(t);
%  figure; plot(t,cf),grid
%  title('CF of the Gaussian distribution N(0,1)')
%
% EXAMPLE2 (PDF/CDF of the Gaussian distribution N(0,1))
%  cf = @(t) cfS_Gaussian(t);
%  x = linspace(-4,4,101);
%  clear options
%  options.N = 2^5;
%  options.SixSigmaRule = 8;
%  result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfS_Gaussian(t);