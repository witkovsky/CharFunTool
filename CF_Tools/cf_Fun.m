function cf = cf_Fun(t,fun,pdf,A,B,tol)
%cf_Fun Computes the characteristic function CF of the TRANSFORMED RANDOM
%  VARIABLE, Y = fun(X), where X is a continuous random variable specified
%  by its PDF function, pdf = @(x) pdf(x) on the interval (A,B). By
%  default here we assume that the transformation function is identity,
%  i.e. fun(x) = x.
%
% DEFINITION:
%  cf_Fun is based on the standard integral representation of the
%  characteristic function of the continuous distribution defined by its 
%  PDF (here PDF is represented by the function handle pdf(x)). Then, 
%    CF(t) = Integral_A^B exp(i*t*fun(x)) * pdf(x) dx. 
%
%  cf_Fun evaluates this integral by using the MATLAB built in function
%  'integral', with precission specified by tolerance tol (default value is 
%  tol = 1e-6).
%
% SYNTAX:
%  cf = cf_Fun(t,fun,pdf,A,B,tol)
%
% INPUTS:
%  t      - real vector, where the characteristic function CF(t) will
%           be evaluated,
%  fun    - function handle used as the transformation function, such that
%           y = fun(x). If empty, default value is fun = @(x) x.
%  pdf    - function handle used as the PDF function with the argument x.
%           If empty, default value is pdf = @(x) exp(-x) (i.e. the
%           exponential distribution on (0,Inf)).  
%  A      - lower limit of the X distribution support. If empty, default
%           value is A = 0.
%  B      - upper limit of the X distribution support. If empty, default
%           value is B = Inf.
%  tol    - relative tolerance parameter used in the built-in Matlab
%           numerical integration algorithm 'integral'. Default value is
%           tol = 1e-6.
%
% OUTPUT:
%  cf     - (complex) vector of the characteristic function values,
%            evalated at the required t, i.e. CF(t).
%
% EXAMPLE1 (CF of Y = X-log(X) with X being exponentially distributed)
%  fun  = @(x) x - log(x);
%  pdf  = @(x) exp(-x);
%  t    = linspace(-100,100,2^10+1)';
%  cf   = cf_Fun(t,fun,pdf);
%  plot(t,real(cf),t,imag(cf));grid
%  title('CF of Y = X-log(X) with X being exponentially distributed')
%
% EXAMPLE2 (PDF/CDF of Y = X-log(X) with X being uniformly distributed)
%  fun  = @(x) x - log(x);
%  pdf  = @(x) 1;
%  A    = 0;
%  B    = 1;
%  cf   = @(t) cf_Fun(t,fun,pdf,A,B);
%  x    = linspace(1,7);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 1;
%  options.SixSigmaRule = 15;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of Y = (X-df)-df*log(X/df) with chi-square distributed X)
%  df   = 7;
%  pdf  = @(x) chi2pdf(x,df);
%  fun  = @(x) (x-df) - df*log(x/df);
%  cf   = @(t) cf_Fun(t,fun,pdf);
%  x    = linspace(0,10);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.N = 2^9;
%  options.xMin = 0;
%  options.SixSigmaRule = 10;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 06-Sep-2018 12:55:53

%% ALGORITHM
%cf = cf_Fun(t,fun,pdf,A,B,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 1
    error(message('VW:cf_Fun:TooFewInputs'));
end
if nargin < 2, fun = []; end
if nargin < 3, pdf = []; end
if nargin < 4, A = []; end
if nargin < 5, B = []; end
if nargin < 6, tol = []; end

if isempty(fun)
    fun = @(x) x;
end

if isempty(pdf)
    pdf = @(x) exp(-x);
end

if isempty(A)
    A = 0;
end

if isempty(B)
    B = Inf;
end

if isempty(tol) 
    tol = 1e-6; 
end
reltol = tol;

%% EVALUATE THE CHARACTERISTIC FUNCTION: cf_Fun(t,pdf)

sz = size(t);
t  = t(:);

if A==0 && B==Inf
cf = integral(@(x) bsxfun(@times,funCF(fun,pdf,t, ...
    (x/(1-x))^2), 2*x/(1-x)^3),0,1,'ArrayValued',true, ...
    'RelTol',reltol);
else
    cf = integral(@(x) funCF(fun,pdf,t,x),A,B, ...
        'ArrayValued',true,'RelTol',reltol);
end

cf(t==0) = 1;
cf = reshape(cf,sz);

end
%% function funCF
function f = funCF(fun,pdf,t,x)
%funCF Integrand function of the integral representation of the
%  characteristic function CF defined by using the fun(x) and pdf(x), 
%  for x and the real (vector) argument t.
%
% SYNTAX:
%  f = funCF(fun,pdf,t,x)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 06-Sep-2018 12:55:53

%% ALGORITHM 
x  = x(:)';
f  = pdf(x) .* exp(1i*t*fun(x));

end