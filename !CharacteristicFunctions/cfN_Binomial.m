function cf = cfN_Binomial(t,n,p,cfX)
%cfN_Binomial(t,n,p) evaluates the characteristic function cf(t) of the
% Binomial distribution with the parameters n (number of trials, n in N)
% and p (success probability, p in [0,1]), i.e. 
%   cf(t) = cfN_Binomial(t,n,p) = (1 - p + p*e^(1i*t))^n
% For more details see [4], and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Binomial_distribution
%
% cfN_Binomial(t,n,p,cfX) evaluates the compound characteristic function
%   cf(t) = cfN_Binomial(-1i*log(cfX(t)),n,p),
% where cfX is function handle of the characteristic function cfX(t) of a
% continuous distribution and/or random variable X.  
% 
% Note that such CF is characteristic function of the compound distribution,
% i.e. distribution of the random variable Y = X_1 + ... + X_N, where X_i ~
% F_X are i.i.d. random variables with common CF cfX(t), and N ~ F_N is
% independent RV with its CF given by cfN(t). 
% 
% SYNTAX
%  cf = cfN_Binomial(t,n,p)
%  cf = cfN_Binomial(t,n,p,cf_X)
%
% EXAMPLE1 (CF of the Binomial distribution with n = 25, p = 0.3)
%  n = 25;  
%  p = 0.3;
%  t = linspace(-15,15,1001);
%  cf = cfN_Binomial(t,n,p);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF the Binomial distribution with n = 25, p = 0.3')
%
% EXAMPLE2 (CF of the compound Binomial-Exponential distribution)
%  n = 25;  
%  p = 0.3;
%  lambda = 10;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  t = linspace(-10,10,501);
%  cf = cfN_Binomial(t,n,p,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the compound Binomial-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-Exponential distribution)
%  n = 25;  
%  p = 0.3;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,5,101);
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
%cf = cfN_Binomial(t,n,p,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, n = []; end
if nargin < 3, p = []; end
if nargin < 4, cfX = []; end

%%
if isempty(n)
    n = 10;
end

if isempty(p)
    p = 0.5;
end
%% Characteristic function of the (compound) Binomial distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf  = (1 - p + p.*expit).^n;
cf = reshape(cf,szt);
cf(t==0) = 1;

end