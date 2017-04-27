function cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol)
%cfX_GeneralizedPareto(t,xi,sigma,theta) Computes the characteristic
% function cf(t) of the Generalized Pareto distribution with parameters xi
% (shape, here xi >= 0), sigma (scale, sigma > 0), and theta (threshold,
% theta >= 0), for real (vector) argument t, i.e.  
%   cf(t) = cfX_GeneralizedPareto(t,xi,sigma,theta);
% The closed-form analytic expression of the characteristic function of the
% Generalized Pareto distribution is unknown. Thus, the characteristic
% function is % numerically evalueted from its definition as suggested in
% [1]. For more details see [1], and also WIKIPEDIA:    
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
% EXAMPLE2 (CDF/PDF of the  Generalized Pareto distribution)
%  xi = 1;
%  sigma = 1;
%  theta = 0;
%  x = linspace(0,300,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.N = 2^15;
%  options.SixSigmaRule = 15;
%  cf = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Generalized Pareto distribution)
%  xi = 1;
%  sigma = 1;
%  theta = 0;
%  lambda = 10;
%  cfX = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
%  cf = @(t) cfN_Poisson(t,lambda,cfX);
%  x = linspace(0,3000,101);
%  prob = [0.9 0.95 0.99];
%  clear options
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
if isempty(xi), theta = 0; end
if isempty(tol), tol = 1e-6; end

sigma(sigma <= 0) = NaN;
xi(xi <= 0) = NaN;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_GeneralizedPareto(t,alpha,beta)

pdfFun = @(x)(1./sigma) .* (1 + (xi./sigma) .* x).^(-(1./xi)-1);
sz = size(t);
t  = t(:);
cf = exp(1i*t*theta) .* cfX_PDF(t,pdfFun,tol);
cf = reshape(cf,sz);
cf(t==0) = 1;

end

%% Function cfX_PDF
function cf = cfX_PDF(t,pdfFun,tol)
%cfX_PDF Characteristic function of the continuos distribution defined by 
%  its PDF function, @(x) pdfFun(x), here we assume x>=0, computed for real 
%  vector argument t.
%
% DEFINITION:
%  cfX_PDF is based on the standard integral representation of the
%  characteristic function of the continuous distribution defined by its 
%  PDF (here PDF is represented by the function handle pdfFun(x) defined 
%  for x >= 0). Then, 
%    CF(t) = Integral_0^inf exp(i*t*x) * pdfFun(x) dx.
%  By using the half-space Fourier integral transformation we get
%    CF(t) = Integral_0^inf (i/t) * exp(-x) * pdfFun(i*x/t) dx.
%  If we define the integrand as 
%    funCF(t,x) = (i/t) * exp(-x) * pdfFun(i*x/t),
%  then by using a stabilizing transformation from [0,inf] to [0,1], we can
%  evaluate the CF by the following (well behaved) integral:
%    CF(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx.
%
%  cfX_PDF evaluates this integral by using the MATLAB built in function
%  'integral', with precission specified by tolerance tol (default value is 
%  tol = 1e-6). For more details see WITKOVSKY (2016).
%
% SYNTAX:
%  cf = cfX_PDF(t,pdfFun,tol)
%
% INPUTS:
%  t      - real vector, where the characteristic function CF(t) will
%           be evaluated,
%  pdfFun - function handle used as the PDF function with the argument x,
%           for x >= 0. However, pdfFun(z) should accept as an input value 
%           any complex matrix z,
%  tol    - relative tolerance parameter used in the built-in Matlab
%           numerical integration algorithm 'integral'. Default value is
%           tol = 1e-6.
%
% OUTPUT:
%  cf     - (complex) vector of the characteristic function values,
%            evalated at the required t, i.e. CF(t).
%
% EXAMPLE1 (CF of the Lognormal distribution with mu = 0,sigma=1)
%  mu = 0;
%  sigma = 1;
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_PDF(t,pdfFun);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Lognormal distribution')
%
% EXAMPLE2 (CDF/PDF of the Lognormal distribution with mu = 0,sigma=1)
%  mu = 0;
%  sigma = 1; 
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  cf = @(t) cfX_PDF(t,pdfFun);
%  clear options
%  options.isForcedSymmetric = true;
%  result = cf2DIST(cf,options);
%
% EXAMPLE3 (CF of the Exponential distribution with mu = 1)
%  pdfFun = @(x) exp(-x);
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_PDF(t,pdfFun);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Exponential distribution')
%
% EXAMPLE4 (CF of the Weibull distribution with a = 1.5, b = 1)
%  a = 1.5;
%  b = 1;
%  pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;;
%  t = linspace(-20,20,2^10+1)';
%  tol = 1e-10;
%  cf = cfX_PDF(t,pdfFun,tol);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE5 (CF of the Log-logistic distribution with a = 1, b = 1)
%  a = 1;
%  b = 8;     
%  pdfFun = @(x)(b./a) .* (x./a).^(b-1) ./ (1 + (x./a).^b).^2;
%  t = linspace(-20,20,2^10+1)';
%  tol = 1e-10;
%  cf = cfX_PDF(t,pdfFun,tol);
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Log-logistic distribution')
%  clear options
%  options.isForcedSymmetric = true;
%  result = cf2DIST(cf,options);
%
% REFERENCES:
%  WITKOVSKY V. (2016). On computing the characteristic functions of 
%  lognormal distribution and its applications. Working Paper.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Apr-2016 18:03:14

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cfX_PDF:TooFewInputs'));
end
if nargin < 2, pdfFun = []; end
if nargin < 3, tol = []; end

if isempty(pdfFun)
    pdfFun = @(x) exp(-x);
end

if isempty(tol) 
    tol = 1e-6; 
end
reltol = tol;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfX_PDF(t,pdfFun)

sz = size(t);
t  = t(:);
cf = ones(size(t));
id = t~=0;
cf(id) = integral(@(x) bsxfun(@times,funCF(pdfFun,t(id),(x/(1-x))^2), ...
    2*x/(1-x)^3),0,1,'ArrayValued',true,'RelTol',reltol);
cf = reshape(cf,sz);

end

%% Function funCF
function f = funCF(pdfFun,t,x)
%funCF Integrand function of the integral representation of the
%  characteristic function CF of the distribution defined by its pdfFun(x),
%  for x >=0, and the real (vector) argument t.
%
% SYNTAX:
%   f = funCF(pdfFun,t,x)
%
% EXAMPLE (Integrand function of the Lognormal(0,1) distribution)
%  mu     = 0;
%  sigma  = 1;
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  t   = 1:5;
%  x   = linspace(0,1);
%  f   = funCF(pdfFun,t,x)
%  plot(x,real(f),x,imag(f))
%
% REFERENCES:
%  WITKOVSKY V. (2016). On computing the characteristic functions of 
%  lognormal distribution and its applications. Working Paper.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Apr-2016 18:03:14

%% ALGORITHM FUNCF
ti  = 1./t(:);
x  = x(:)';
ot = ones(size(ti));
ox = ones(size(x));

funT    = @(x,ti)(1i*ti*ox) .* pdfFun(1i*ti*x) .* exp(-ot*x);
f       = funT(x,ti);

end