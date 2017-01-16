%function cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol)
%cfX_GeneralizedPareto(t,xi,sigma,theta) Computes the characteristic
% function cf(t) of the Generalized Pareto distribution with parameters xi
% (shape, here xi >= 0), sigma (scale, sigma > 0), and theta (threshold,
% theta >= 0), for real (vector) argument t, i.e.  
%   cf(t) = cfX_GeneralizedPareto(t,xi,sigma,theta);
% The closed-form analytic expression of the characteristic function of the
% Generalized Pareto distribution is unknown. Thus, the characteristic
% function is numerically evalueted from its definition as suggested in
% [3]. For more details see [3], and also WIKIPEDIA:    
% https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
%  
% SYNTAX:
%  cf = cfX_GeneralizedPareto(t,xi,sigma,theta)
%  cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol)
%
% EXAMPLE1 (CF of the Generalized Pareto distribution)
%  xi = 1;
%  sigma = 1;
%  theta = 1;
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_GeneralizedPareto(t,xi,sigma,theta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Generalized Pareto distribution')
%
% EXAMPLE2 (CDF/PDF of the  Generalized Pareto distribution)
%  xi = 1;
%  sigma = 1;
%  theta = 0;
%  x = linspace(0,300,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.N = 2^15;
%  options.SixSigmaRule = 15;
%  cf = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Generalized Pareto distribution)
%  xi = 1;
%  sigma = 1;
%  theta = 0;
%  lambda = 10;
%  cfX = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,5000,101);
%  prob = [0.9 0.95 0.99];
%  clear options
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
%cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol);