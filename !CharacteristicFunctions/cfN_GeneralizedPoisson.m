function [cf,N] = cfN_GeneralizedPoisson(t,a,p,cfX)
%cfN_GeneralizedPoisson(t,a,p) evaluates the characteristic function cf(t)
% of the Generalized-Poisson distribution, with the parameters a (parameter
% of variable mean, a > 0), and p (success probability, p in [0,1]), i.e. 
%   cf(t) = cfN_GeneralizedPoisson(t,a,p) 
%         = exp(a*(sum_{j=1}^Inf ((p*j)^(j-1)*e^(-p*j)/j!)*e^(1i*t*j)-1));
% For more details see [4], p. 93.
%
%  The Generalized-Poisson distribution is equivalent with the Borel-Tanner
%  distribution with parameters (p,m), see WIMMER & ALTMANN (1999), where 
%  m ~ Poisson(a).
%
% cfN_GeneralizedPoisson(t,a,p,cfX) evaluates the compound characteristic
% function 
%   cf(t) = cfN_GeneralizedPoisson(-1i*log(cfX(t)),a,p),
% where cfX is function handle of the characteristic function cfX(t) of a
% continuous distribution and/or random variable X.  
% 
% Note that such CF is characteristic function of the compound distribution,
% i.e. distribution of the random variable Y = X_1 + ... + X_N, where X_i ~
% F_X are i.i.d. random variables with common CF cfX(t), and N ~ F_N is
% independent RV with its CF given by cfN(t). 
%
% SYNTAX
%  cf = cfN_GeneralizedPoisson(t,a,p)
%  cf = cfN_GeneralizedPoisson(t,a,p,cfX)
%
% EXAMPLE1 (CF of the Generalized-Poisson distribution with a=10, p=0.5)
%  a = 10;
%  p = 0.5;
%  t = linspace(-10,10,501);
%  cf = cfN_GeneralizedPoisson(t,a,p);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Generalized-Poisson distribution with a=10, p=0.5')
%
% EXAMPLE2 (CF compound Generalized-Poisson-Exponential distribution)
%  a = 10;
%  p = 0.5;
%  lambda = 5;
%  t = linspace(-5,5,501);
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = cfN_GeneralizedPoisson(t,a,p,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF compound Generalized-Poisson-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Generalized-Poisson-Exponential distr.)
%  a = 10;
%  p = 0.5;
%  lambda = 5;
%  x = linspace(0,15,101);
%  prob = [0.9 0.95 0.99];
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_GeneralizedPoisson(t,a,p,cfX);
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
%     output quantity in linear measurement models. Acta IMEKO, 5(3),32-44. 
% [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
%     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN
%     3-87773-025-6. 

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%[cf,N] = cfN_GeneralizedPoisson(t,a,p,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, cfX = []; end

%% Characteristic function of the Generalized-Poisson distribution
szt = size(t);
t   = t(:);
 
if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

exppt = exp(-p) * expit;
c     = exppt;
cf    = c;
N     = ceil(-36.84/(1+log(p)-p));
for j = 1:N
    c  = p * ((1+j)/j)^(j-1) .* exppt .* c;
    cf = cf + c;
end
cf = exp(a * (cf-1));
cf = reshape(cf,szt);
cf(t==0) = 1;

end