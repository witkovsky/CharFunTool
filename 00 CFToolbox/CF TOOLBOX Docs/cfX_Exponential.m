%function cf = cfX_Exponential(t,lambda)
%cfX_Exponential(t,lambda) evaluates the characteristic function cf(t) of
% the Exponential distribution with the parameter lambda (rate, lambda > 0),
% i.e.   
%   cf(t) = cfX_Exponential(t,lambda) = lambda / (lambda - 1i*t)
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Exponential_distribution
%
% SYNTAX
%  cf = cfX_Exponential(t,lambda)
%
% EXAMPLE1 (CF of the Exponential distribution with lambda = 5)
%  lambda = 5;  
%  t = linspace(-50,50,501);
%  cf = cfX_Exponential(t,lambda);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Exponential distribution with lambda = 5')
%
% EXAMPLE2 (PDF/CDF of the Exponential distribution with lambda = 5)
%  lambda = 5;  
%  cf = @(t) cfX_Exponential(t,lambda);
%  x  = linspace(0,1.5,101);
%  clear options
%  options.xMin = 0;
%  options.SixSigmaRule = 8;
%  result = cf2DistGP(cf,x,[],options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-Exponential distribution)
%  n = 25;  
%  p = 0.3;
%  lambda = 5;
%  cfX  = @(t) cfX_Exponential(t,lambda);
%  cf   = @(t) cfN_Binomial(t,n,p,cfX);
%  x    = linspace(0,5,101);
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
% [3] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Working Paper.
%     Insurance: Mathematics and Economics. 
% [4] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models.  Working Paper.
%     Journal of Statistical Software.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfX_Exponential(t,lambda);