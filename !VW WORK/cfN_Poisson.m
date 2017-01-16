function cf = cfN_Poisson(t,lambda,cfX)
%cfN_Poisson(t,lambda) evaluates the characteristic function cf(t) of the
% Poisson distribution with the rate parameter lambda > 0, i.e. 
%   cf(t) = cfN_Poisson(t,lambda) = exp(lambda*(e^(1i*t)-1))
% For more details see [4], and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Poisson_distribution.
%
% cfN_Poisson(t,lambda,cfX) evaluates the compound characteristic function
%   cf(t) = cfN_Poisson(-1i*log(cfX(t)),lambda),
% where cfX is function handle of the characteristic function cfX(t) of a
% continuous distribution and/or random variable X.  
% 
% Note that such CF is characteristic function of the compound distribution,
% i.e. distribution of the random variable Y = X_1 + ... + X_N, where X_i ~
% F_X are i.i.d. random variables with common CF cfX(t), and N ~ F_N is
% independent RV with its CF given by cfN(t). 
% 
% SYNTAX
%  cf = cfN_Poisson(t,lambda)
%  cf = cfN_Poisson(t,lambda,cfX)
%
% EXAMPLE1 (CF of the Poisson distribution with the parameter lambda = 10)
%  lambda = 10;  
%  t = linspace(-10,10,501);
%  cf = cfN_Poisson(t,lambda);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('Characteristic function of the Poisson distribution')
%
% EXAMPLE2 (CF compound Poisson & Exponential distribution)
%  lambda1 = 10;
%  lambda2 = 5;
%  t = linspace(-10,10,501);
%  cfX = @(t) cfX_Exponential(t,lambda2);
%  cf = cfN_Poisson(t,lambda1,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF compound Poisson & Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Exponential distribution)
%  lambda1 = 10;
%  lambda2 = 5;
%  cfX = @(t) cfX_Exponential(t,lambda2);
%  cf = @(t) cfN_Poisson(t,lambda1,cfX);
%  x = linspace(0,8,101);
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
%     output quantity in linear measurement models. Acta IMEKO, 5(3),32-44.  
% [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
%     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN
%     3-87773-025-6. 

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfN_Poisson(t,lambda,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 2, lambda = []; end
if nargin < 3, cfX = []; end

%%
if isempty(lambda)
    lambda = 1;
end

%% Characteristic function of the (compound) Poisson distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf = lambda * (expit-1);
cf = exp(cf);
cf = reshape(cf,szt);
cf(t==0) = 1;

end