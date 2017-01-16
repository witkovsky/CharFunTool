%function cf = cfX_LogNormal(t,mu,sigma,tol)
%cfX_LogNormal(t,mu,sigma) Computes the characteristic function cf(t) of
% the LogNormal distribution with parameters mu (real) and sigma > 0,
% computed for real (vector) argument t, i.e.
%   cf(t) = cfX_LogNormal(t,mu,sigma);
% For more details see also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Log-normal_distribution
%  
% SYNTAX:
%  cf = cfX_LogNormal(t,mu,sigma)
%  cf = cfX_LogNormal(t,mu,sigma,tol)
%
% EXAMPLE1 (CF of the LogNormal distribution with mu = 0,sigma = 1)
%  mu = 0;
%  sigma = 1;
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_LogNormal(t,mu,sigma);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogNormal distribution')
%
% EXAMPLE2 (CDF/PDF of the LogNormal distribution with mu = 0,sigma = 1)
%  mu = 0;
%  sigma = 1;
%  x = linspace(0,15,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^13;
%  options.SixSigmaRule = 8;
%  cf = @(t) cfX_LogNormal(t,mu,sigma);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-LogNormal distribution)
%  mu = 0;
%  sigma = 1;
%  lambda = 10;
%  cfX = @(t)cfX_LogNormal(t,mu,sigma);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,70,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% FURTHER DETAILS:
% In probability theory, a log-normal (or lognormal) distribution is a
% continuous probability distribution of a random variable whose logarithm
% is normally distributed. The lognormal distribution is defined for x in
% (0,+inf) by its PDF/CDF/CF, as follows
%  pdf(x) = \frac{1}{x\sigma\sqrt{2\pi}}\ 
%            e^{-\frac{\left(\ln x-\mu\right)^2}{2\sigma^2}},
%  cdf(x) = \frac12 + \frac12\,\mathrm{erf}\Big[\frac{\ln x-\mu}
%            {\sqrt{2}\sigma}\Big],
%  cf(t)  = \sum_{n=0}^{\infty}\frac{(it)^n}{n!}e^{n\mu+n^2\sigma^2/2}.
% As noted, this representation is asymptotically divergent but sufficient
% for numerical purposes. 
%
% cfX_LogNormal is based on the standard integral representation of the
% characteristic function of the lognormal distribution, i.e.
%  cf(t) = Integral_0^inf exp(i*t*x) * PDF(x) dx.
% By using the half-space Fourier integral transformation we get
%  cf(t) = Integral_0^inf (i/t) * exp(-x) * PDF(i*x/t) dx.
% If we define the integrand as funCF(t,x) = (i/t) * exp(-x) * PDF(i*x/t),
% then by using a stabilizing transformation from [0,inf] to [0,1], we can
% evaluate the CF by the following (well behaved) integral:
%  cf(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx.
%
% cfX_LogNormal evaluates this integral by using the MATLAB built in
% function 'integral', with precission specified by tolerance tol (default
% value is tol = 1e-6). 
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
%cf = cfX_LogNormal(t,mu,sigma,tol);