function cf = cfX_PDF(t,pdfFun,method,tol)
%cfX_PDF Computes the characteristic function of the continuos
% distribution defined by its PDF function, pdfFun = @(x) pdf(x), here we
% assume x>=0, computed for real vector argument t.
%
% DEFINITION:
%  cfX_PDF is based on the standard integral representation of the
%  characteristic function of the continuous distribution defined by its 
%  PDF (here PDF is represented by the function handle pdfFun(x) defined 
%  for x >= 0). Then, 
%    CF(t) = Integral_0^inf exp(i*t*x) * pdfFun(x) dx.
%  Alternatively, by using the half-space Fourier Integral Transformation
%  (FIT), which could improve the highly oscillatory behaviour of the
%  integrand function, we get    
%    CF(t) = Integral_0^inf (i/t) * exp(-x) * pdfFun(i*x/t) dx.
%  If we define the integrand as 
%    funCF(t,x) = (i/t) * exp(-x) * pdfFun(i*x/t),
%  then by using a stabilizing transformation from [0,inf] to [0,1], we can
%  evaluate the CF by the following (well behaved) integral:
%    CF(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx.
%
%  Selection of the proper method (standard definition or the integral
%  transformation FIT) depends on the distribution and the particular form
%  of its PDF.   
%
%  cfX_PDF evaluates this integral by using the MATLAB built in function
%  'integral', with precission specified by tolerance tol (default value is 
%  tol = 1e-6).
%
% SYNTAX:
%  cf = cfX_PDF(t,pdfFun,method,tol)
%
% INPUTS:
%  t      - real vector, where the characteristic function CF(t) will
%           be evaluated,
%  pdfFun - function handle used as the PDF function with the argument x,
%           for x >= 0. However, pdfFun(z) should accept as an input value 
%           any complex matrix z,
%  method - set the method used to compute the characteristic function,
%           either method = 'def' (or the standard definition of CF by using
%           the pdf) or method = 'fit' (by using the half-space Fourier
%           integral transform of the pdf). Default value is method =
%           'def'. 
%  tol    - relative tolerance parameter used in the built-in Matlab
%           numerical integration algorithm 'integral'. Default value is
%           tol = 1e-6.
%
% OUTPUT:
%  cf     - (complex) vector of the characteristic function values,
%            evalated at the required t, i.e. CF(t).
%
% EXAMPLE1 (CF of the Exponential distribution with mu = 1)
%  pdfFun = @(x) exp(-x);
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_PDF(t,pdfFun);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Exponential distribution')
%
% EXAMPLE2 (CF of the LogNormal distribution with mu = 0, sigma = 1)
%  mu     = 0;
%  sigma  = 1;
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  method = 'fit';
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_PDF(t,pdfFun,method);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogNormal distribution')
%
% EXAMPLE3 (CDF/PDF of the LogNormal distribution with mu = 0, sigma = 1)
%  mu     = 0;
%  sigma  = 1; 
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  method = 'fit';
%  cf = @(t) cfX_PDF(t,pdfFun,method);
%  clear options
%  options.xMin = 0;
%  result = cf2DistGP(cf,[],[],options);
%
% EXAMPLE4 (CF of the Weibull distribution with a = 1.5, and small b<1)
%  a      = 1.5;
%  b      = 0.5;
%  pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  method = 'fit';
%  t  = linspace(-20,20,2^10+1)';
%  cf = cfX_PDF(t,pdfFun,method);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE5 (CF of the Weibull distribution with a = 1.5, and large b > 1)
%  a      = 1.5;
%  b      = 3.5;
%  pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  method = 'def';
%  t  = linspace(-10,10,2^10+1)';
%  cf = cfX_PDF(t,pdfFun,method);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
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
% Ver.: 21-Jul-2017 20:53:22

%% ALGORITHM
%cf = cfX_PDF(t,pdfFun,method,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cfX_PDF:TooFewInputs'));
end
if nargin < 2, pdfFun = []; end
if nargin < 3, method = []; end
if nargin < 4, tol = []; end

if isempty(pdfFun)
    pdfFun = @(x) exp(-x);
end

if isempty(method) 
    method = 'definition'; 
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

switch lower(method)
    case {'def','definition','standard','direct'}
        cf(id) = integral(@(x) bsxfun(@times,funCFdef(pdfFun,t(id), ...
            (x/(1-x))^2), 2*x/(1-x)^3),0,1,'ArrayValued',true, ...
            'RelTol',reltol);
    case {'fit','fourier','transform'}
        cf(id) = integral(@(x) bsxfun(@times,funCFfit(pdfFun,t(id), ...
            (x/(1-x))^2), 2*x/(1-x)^3),0,1,'ArrayValued',true, ...
            'RelTol',reltol);
    otherwise
        cf(id) = integral(@(x) bsxfun(@times,funCFdef(pdfFun,t(id), ...
            (x/(1-x))^2), 2*x/(1-x)^3),0,1,'ArrayValued',true, ...
            'RelTol',reltol);
end
cf = reshape(cf,sz);

end
%% Function funCFdef
function f = funCFdef(pdfFun,t,x)
%funCFdef Integrand function of the integral representation of the
%  characteristic function CF defined by using the pdfFun(x), 
%  for x >=0, and the real (vector) argument t.
%
% SYNTAX:
%  f = funCFdef(pdfFun,t,x)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Jul-2017 20:53:22

%% ALGORITHM 
x  = x(:)';
f  = pdfFun(x) .* exp(1i*t*x);
end
%% Function funCFfit
function f = funCFfit(pdfFun,t,x)
%funCFfit Integrand function of the integral representation of the
%  characteristic function CF of the distribution defined by using the
%  half-space Fourier Integral Transform (FIT) of its pdfFun(x), 
%  for x >=0, and the real (vector) argument t.
%
% SYNTAX:
%    f = funCFfit(pdfFun,t,x)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Jul-2017 20:53:22

%% ALGORITHM 
ti  = 1./t(:);
x  = x(:)';
ot = ones(size(ti));
ox = ones(size(x));
f  = (1i*ti*ox) .* pdfFun(1i*ti*x) .* exp(-ot*x);
end