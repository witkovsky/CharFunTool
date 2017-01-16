%function cf = cfN_PolyaEggenberger(t,a,b,m,cfX)
%cfN_PolyaEggenberger(t,a,b,m) evaluates the characteristic function cf(t)
% of the Polya-Eggenberger distribution, with the parameters a (a real), b
% (b real), and m (m integer), i.e.  
%   cf(t) = cfN_PolyaEggenberger(t,a,b,m) = 2F1(-m,a,a+b,1-e^(1i*t));
% where 2F1 denotes the Gauss hypergeometric function. For more details see
% [4], p. 525. 
%
% SYNTAX
%  cf = cfN_PolyaEggenberger(t,a,b,m)
%  cf = cfN_PolyaEggenberger(t,a,b,m,cfX)
%
% EXAMPLE1 (CF of Polya-Eggenberger distribution with a=2.2, b=3.3, m=4)
%  a = 2.2;
%  b = 3.3;
%  m = 4; 
%  t = linspace(-15,15,1001);
%  cf = cfN_PolyaEggenberger(t,a,b,m);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Polya-Eggenberger distribution with a=2.2, b=3.3, m=4')
%
% EXAMPLE2 (CF of the compound Polya-Eggenberger-Exponential distribution)
%  a = 2.2;
%  b = 3.3;
%  m = 4; 
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  t = linspace(-50,50,501);
%  cf = cfN_PolyaEggenberger(t,a,b,m,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the compound Polya-Eggenberger-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Polya-EggenbergerExponential distr.)
%  a = 2.2;
%  b = 3.3;
%  m = 4; 
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_PolyaEggenberger(t,a,b,m,cfX);
%  x = linspace(0,2.5,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Preprint submitted
%     to Insurance: Mathematics and Economics.
% [2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models. Preprint submitted
%     to Journal of Statistical Software.
% [3] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44. 
% [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
%     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN
%     3-87773-025-6. 

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfN_PolyaEggenberger(t,a,b,m,cfX);