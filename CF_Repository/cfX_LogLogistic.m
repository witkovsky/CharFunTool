function cf = cfX_LogLogistic(t,alpha,beta,tol)
%cfX_LogLogistic(t,alpha,beta) 
% Characteristic function cf(t) of the LogLogistic distribution with
% parameters alpha (scale, alpha > 0) and beta (shape, beta > 0), for real
% (vector) argument t, i.e.  
%   cf(t) = cfX_LogLogistic(t,alpha,beta);
%
% NOTE: 
% This is an experimental algorithm which is not fully documented and that
% may provide unreliable results for a particular combination of parameters
% and / or may cause an unexpected failure.
%
% The closed-form analytic expression of the characteristic function of the
% LogLogistic distribution is unknown. Thus, the characteristic function is
% numerically evalueted from its definition. 
% 
% WIKIPEDIA:    
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
%  tol = 1e-3;
%  cf = cfX_LogLogistic(t,alpha,beta,tol);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogLogistic distribution')
%
% EXAMPLE2 (CDF/PDF of the Loglogistic distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta = 1;
%  tol = 1e-3;
%  cf = @(t) cfX_LogLogistic(t,alpha,beta,tol);
%  x = linspace(0,70,101);
%  prob = [];
%  clear options
%  options.xMin = 0;
%  options.SixSigmaRule = 15;
%  options.N = 2^14;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (CDF/PDF of the Loglogistic distribution with alpha=1, beta=0.75)
% % NOTE: Here cf2DistBV returns more reliable result for the heavy tailed
% %       distribution with beta < 1
%  alpha = 1;
%  beta = 0.75;
%  tol = 1e-3;
%  cf = @(t) cfX_LogLogistic(t,alpha,beta,tol);
%  x = linspace(0,70,101);
%  prob = [];
%  clear options
%  options.xMin = 0;
%  result = cf2DistBV(cf,x,prob,options)
%
% EXAMPLE4 (PDF/CDF of the compound Poisson-Loglogistic distribution)
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

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 12:25:48

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
if isempty(tol), tol = 1e-4; end

alpha(alpha <= 0) = NaN;
beta(beta <= 0) = NaN;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_LogLogistic(t,alpha,beta)

pdfFun = @(x)(beta./alpha) .* (x./alpha).^(beta-1) ./ (1 + (x./alpha).^beta).^2;

sz = size(t);
t  = t(:);

if beta < 2
    method = 'fit';
else
    method = 'def';
end

cf = cfX_PDF(t,pdfFun,method,tol);
cf = reshape(cf,sz);
cf(t==0) = 1;

end