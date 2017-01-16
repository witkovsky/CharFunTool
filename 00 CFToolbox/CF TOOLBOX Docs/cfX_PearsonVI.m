%function cf = cfX_PearsonVI(t,alpha,beta)
%cfX_PearsonVI(t,alpha,beta) Computes the characteristic function of the
% Pearson type VI distribution with parameters alpha > 0 (shape) and beta >
% 0 (scale), computed for real vector argument t, i.e. 
%   cf(t) = cfX_PearsonVI(t,alpha,beta)
%         = (gamma(alpha+beta)./gamma(beta)) .* U(alpha,1-beta,-1i*t),
% where U(a,b,z) denotes the confluent hypergeometric function of the
% second kind. See also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Pearson_distribution
%
% SYNTAX:
%  cf = cfX_PearsonVI(t,alpha,beta)
%
% INPUTS:
%  t      -  vector of real values where the CF is evaluated, i.e. CF(t),
%  alpha  -  scalar shape parameter of the PearsonVI distribution, alpha > 0,
%  beta   -  scalar scale parameter of the PearsonVI distribution, beta > 0,
%
% OUTPUT:
% cf     -   calculated values of the characteristic function CF(t)
%
% EXAMPLE1: (CF of the PearsonVI distribution)
%  alpha = 3/2;
%  beta = 2/3;
%  t = linspace(-10,10,1001);
%  cf = cfX_PearsonVI(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the PearsonVI distribution with alpha = 3/2, beta = 2/3')
%
% EXAMPLE2 (PDF/CDF of the PearsonVI distribution)
%  alpha = 3/2;
%  beta = 2/3;
%  x = linspace(0,200,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^14;
%  cf = @(t) cfX_PearsonVI(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-PearsonVI distribution)
%  n = 25;  
%  p = 0.3;
%  alpha = 3/2;
%  beta = 2/3;
%  cfX = @(t) cfX_PearsonVI(t,alpha,beta);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,10000,101);
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
%cf = cfX_PearsonVI(t,alpha,beta);
