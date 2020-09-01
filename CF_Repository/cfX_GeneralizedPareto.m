function cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol)
%cfX_GeneralizedPareto(t,xi,sigma,theta) 
% Characteristic function cf(t) of the Generalized Pareto distribution with
% the shape parameter xi (any real value, however here we assume 0 <= xi <
% 10), scale parameter sigma (sigma > 0), and the threshold parameter theta
% (real). cf(t) is evaluated for real (vector) argument t, i.e.   
%   cf(t) = cfX_GeneralizedPareto(t,xi,sigma,theta);
%
% NOTE: 
% This is an experimental algorithm which is not fully documented and that
% may provide unreliable results for a particular combination of parameters
% and / or may cause an unexpected failure.
%
% The closed-form analytic expression of the characteristic function of the
% Generalized Pareto distribution is unknown. Thus, the characteristic
% function is numerically evalueted from its definition. 
%
% WIKIPEDIA:
% https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
%  
% SYNTAX:
%  cf = cfX_GeneralizedPareto(t,xi,sigma,theta)
%  cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol)
%
% EXAMPLE1 (CF of the Generalized Pareto distribution)
%  xi = 1;
%  sigma = 1;
%  theta = 1;
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_GeneralizedPareto(t,xi,sigma,theta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Generalized Pareto distribution')
%
% EXAMPLE2 (CDF/PDF of the Generalized Pareto distribution by cf2DistGP)
% % NOTE: Here cf2DistGP returns ureliable result for the heavy tailed
% %       distribution with xi > 1
%  xi = 2;
%  sigma = 1;
%  theta = 0;
%  x = linspace(0,5,101);
%  prob = [];
%  clear options
%  options.xMin = 0;
%  options.N = 2^12;
%  options.SixSigmaRule = 15;
%  cf = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (CDF/PDF of the Generalized Pareto distribution by cf2DistBV)
% % NOTE: Here cf2DistBV returns  more reliable result for the heavy tailed
% %       distribution with xi > 1
%  xi = 2;
%  sigma = 1;
%  theta = 0;
%  cf = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  x = linspace(0,5,101);
%  prob = [];
%  clear options
%  options.xMin = 0;
%  result = cf2DistBV(cf,x,prob,options)
%
% EXAMPLE4 (CDF/PDF of the Generalized Pareto distribution)
%  xi = 0.5;
%  sigma = 1;
%  theta = 0;
%  x = linspace(0,5,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^12;
%  cf = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE5 (PDF/CDF of the compound Poisson-Generalized Pareto distribution)
%  xi = 1;
%  sigma = 1;
%  theta = 0;
%  lambda = 10;
%  cfX = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,3000,101);
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
%cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cfX_GeneralizedPareto:TooFewInputs'));
end
if nargin < 2, xi = []; end
if nargin < 3, sigma = []; end
if nargin < 4, theta = []; end
if nargin < 5, tol = []; end

if isempty(xi), xi = 1; end
if isempty(sigma), sigma = 1; end
if isempty(theta), theta = 0; end
if isempty(tol), tol = 1e-6; end

sigma(sigma <= 0) = NaN;
%xi(xi <= 0) = NaN;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_GeneralizedPareto(t,alpha,beta)

if xi == 0
    pdfFun = @(x)(1./sigma) .* exp(-x./sigma);
elseif xi > 0
    pdfFun = @(x)(1./sigma) .* (1 + (xi./sigma) .* x).^(-(1./xi)-1) ;
else
    pdfFun = @(x)(1./sigma) .* (1 + (xi./sigma) .* x).^(-(1./xi)-1) .* ...
        (x >= 0) .* (x <= -1./xi);
end

sz = size(t);
t  = t(:);

if xi > 0
    method = 'fit';
else
    method = 'def';
end

if theta == 0
    cf = cfX_PDF(t,pdfFun,method,tol);
else
    cf = exp(1i*t*theta) .* cfX_PDF(t,pdfFun,method,tol);
end
cf = reshape(cf,sz);
cf(t==0) = 1;

end