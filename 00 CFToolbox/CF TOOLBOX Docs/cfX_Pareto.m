%function cf = cfX_Pareto(t,alpha,sigma,type,tol)
%cfX_Pareto(t,alpha,sigma,type,tol) Computes the characteristic function of
% the Pareto distribution with parameters alpha > 0 (shape), sigma > 0
% (scale, sometimes denoted by x_m - which is the minimum value of the
% distribution support, i.e. x >= x_m), and the parameter type ('typeI' or
% 'typeII'), computed for real vector argument t, i.e.
%  cf(t) = cfX_Pareto(t,alpha,sigma,type);
% For more details see [1], and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Pareto_distribution
%
% SYNTAX:
%  cf = cfX_Pareto(t,alpha)
%  cf = cfX_Pareto(t,alpha,sigma,type)
%
% INPUTS:
%  t      -  vector of real values where the CF is evaluated, i.e. CF(t),
%  alpha  -  scalar shape parameter of the Pareto distribution, alpha > 0,
%  sigma  -  scalar scale parameter of the Pareto distribution, sigma > 0,
%  type   -  type of the Pareto distribution: 'typeI' or 'typeII',
%  tol    -  relative tolerance of the integral function (tol = 1e-6).
%
% OUTPUT:
% cf     -   calculated values of the characteristic function CF(t)
%
% FURTHER DETAILS:
% The Pareto Type I distribution (also known as European Pareto) is
% characterized by a scale parameter sigma and a shape parameter alpha,
% which is known as the tail index. The pdf of X is defined by pdf_X(x) =
% alpha * sigma^alpha / x^(alpha+1), for x >= sigma, otherwise it is zero.
% When this distribution is used to model the distribution of wealth, then
% the parameter alpha is called the Pareto index.
%
% The Type II Pareto distribution (also known as Lomax or American Pareto
% distribution, type = 'typeII' or type = 'Amer') is shifted to zero.
% The default type is 'TypeI' (type = 'typeI' or type = 'Eur').
% The CFs of the Pareto type distributions are defined by
% cf_{typeI}(t)  = alpha * exp(1i*sigma*t) * U(1,1-alpha,-1i*sigma*t)
% cf_{typeII}(t) = alpha * U(1,1-alpha,-1i*sigma*t), where U(a,b,z) denotes
% the confluent hypergeometric function of the second kind. 
%
% A special case of Pareto TypeI distribution is the Pareto 1 distribution
% with sigma = 1, cf_{1}(t) = alpha * exp(1i*t) * U(1,1-alpha,-1i*t).
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Pareto_distribution
%
% EXAMPLE1 (CF of the Pareto 1 distribution)
%  alpha = 3/2;  
%  t = linspace(-10,10,1001);
%  cf = cfX_Pareto(t,alpha);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Pareto 1 distribution with alpha = 3/2')
%
% EXAMPLE2 (CF of the Pareto Type I/EUR distribution)
%  alpha = 3/2;  
%  sigma  = 2;
%  type = 'eur';
%  t = linspace(-10,10,1001);
%  cf = cfX_Pareto(t,alpha,sigma,type);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Pareto EUR distribution with alpha = 3/2, sigma = 2')
%
% EXAMPLE3 (CF of the Pareto Type II/AMER distribution)
%  alpha = 3/2;  
%  sigma  = 2;
%  type = 'amer';
%  t = linspace(-10,10,1001);
%  cf = cfX_Pareto(t,alpha,sigma,type);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Pareto AMER distribution with alpha = 3/2, sigma = 2')
%
% EXAMPLE4 (PDF/CDF of the compound Poisson-Pareto distribution)
%  alpha = 3/2;  
%  sigma  = 2;
%  type = 'eur';
%  lambda = 10;
%  cfX = @(t) cfX_Pareto(t,alpha,sigma,type);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,500,101);
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
%cf = cfX_Pareto(t,alpha,sigma,type,tol);