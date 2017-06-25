function cf = cfX_LogLogistic(t,alpha,beta,tol)
%cfX_LogLogistic(t,alpha,beta) Computes the characteristic function cf(t)
% of the LogLogistic distribution with parameters alpha (scale, alpha > 0)
% and beta (shape, beta > 0), for real (vector) argument t, i.e. 
%   cf(t) = cfX_LogLogistic(t,alpha,beta);
% The closed-form analytic expression of the characteristic function of the
% LogLogistic distribution is unknown. Thus, the characteristic function is
% numerically evalueted from its definition as suggested in [1]. For more
% details see [1], and also WIKIPEDIA:    
% https://en.wikipedia.org/wiki/Log-logistic_distribution
%  
% SYNTAX:
%  cf = cfX_LogLogistic(t,alpha,beta)
%  cf = cfX_LogLogistic(t,alpha,beta,tol)
%
% EXAMPLE1 (CF of the LogLogistic distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta = 1;
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_LogLogistic(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogLogistic distribution')
%
% EXAMPLE2 (CDF/PDF of the  Loglogistic distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta = 1;
%  x = linspace(0,700,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.xMax = 1000;
%  options.N = 2^14;
%  cf = @(t) cfX_LogLogistic(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Loglogistic distribution)
%  alpha = 1;
%  beta = 1;
%  lambda = 10;
%  cfX = @(t)cfX_LogLogistic(t,alpha,beta);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,3500,101);
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
%cf = cfX_LogLogistic(t,alpha,beta,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cfX_LogLogistic:TooFewInputs'));
end
if nargin < 2, alpha = []; end
if nargin < 3, beta = []; end
if nargin < 4, tol = []; end

if isempty(alpha), alpha = 1; end
if isempty(beta), beta = 1; end
if isempty(tol), tol = 1e-6; end

alpha(alpha <= 0) = NaN;
beta(beta <= 0) = NaN;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_LogLogistic(t,alpha,beta)

pdfFun = @(x)(beta./alpha) .* (x./alpha).^(beta-1) ./ (1 + (x./alpha).^beta).^2;
cf = cfX_PDF(t,pdfFun,tol);

end