function cf = cfX_Weibull(t,alpha,beta,tol)
%cfX_Weibull(t,alpha,beta) Computes the characteristic function cf(t)
%  of the Weibull distribution with the parameters alpha(scale parameter,
%  elsewhere denoted also as lambda, alpha > 0) and beta (shape parameter,
%  elsewhere denoted as k, beta > 0), for real (vector) argument t, i.e.
%  cf(t) = cfX_Weibull(t,alpha,beta);
%
%  The closed-form analytic expression of the characteristic function of
%  the Weibull distribution is unknown. The series expansion of the CF,
%  numerically stable for abs(t)<1, is known. For more details see also
%  WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Weibull_distribution
%  
% SYNTAX:
%  cf = cfX_Weibull(t,alpha,beta)
%  cf = cfX_Weibull(t,alpha,beta,tol)
%
% EXAMPLE1 (CF of the Weibull distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta  = 1;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE2 (CF of the Weibull distribution with alpha=1, beta<1)
%  alpha = 1;
%  beta  = 0.5;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE3 (CF of the Weibull distribution with alpha=1, beta>1)
%  alpha = 1;
%  beta  = 5.5;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE4 (CDF/PDF of the  Weibull distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta  = 1;
%  x    = linspace(0,7,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^10;
%  cf = @(t) cfX_Weibull(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE5 (PDF/CDF of the compound Poisson-Weibull distribution)
%  alpha = 1;
%  beta = 1;
%  lambda = 10;
%  cfX = @(t)cfX_Weibull(t,alpha,beta);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,35,101);
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
%cf = cfX_Weibull(t,alpha,beta,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cfX_Weibull:TooFewInputs'));
end
if nargin < 2, alpha = []; end
if nargin < 3, beta = []; end
if nargin < 4, tol = []; end

if isempty(alpha), alpha = 1; end
if isempty(beta), beta = 1; end
if isempty(tol), tol = 1e-6; end

alpha(alpha <= 0) = NaN;
beta(beta <= 0) = NaN;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_Weibull(t,alpha,beta)

pdfFun   = @(x) (x./alpha).^(beta-1) .* ...
           exp(-((x./alpha).^beta)) .* beta./alpha;
szt      = size(t);
t        = t(:);
cf       = zeros(size(t));

if beta > 1
    method = 'def';
    cf     = cfX_PDF(t,pdfFun,method,tol);
elseif beta == 1
    method = 'fit';
    id     = abs(t/alpha) <= 0.99;
    cf(id) = FoxWrightPsi(1,1/beta,[],[],-1i*alpha*t(id));
    id     = abs(t/alpha) > 0.99;
    cf(id) = cfX_PDF(t(id),pdfFun,method,tol);
else
    method = 'fit';
    cf     = cfX_PDF(t,pdfFun,method,tol);
end

cf       = reshape(cf,szt);
cf(t==0) = 1;

end