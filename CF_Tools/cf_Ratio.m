function cf = cf_Ratio(t,cfX1,pdfX2,a,b,tol)
%% cf_Ratio 
%  Characteristic function of a ratio of two independent (continuous)
%  random variables, i.e. Y = X1 / X2. The distribution of the first random
%  variable X1 is given by its characteristic function (CF), say cf_X1(t),
%  and the distribution of the second random variable is given by its
%  probability density function (PDF), say pdf_X2(x). 
% 
%  Then, the characteristic function of the ratio of two independent
%  random variables, say cf_Y(t), is given by
%   cf_Y(t) = integral(@(x) cf_X1(t/x) .* pdf_X2(x),a,b),
%  where a and b are limits of integration (interval representing the
%  support of the distribution of X2). 
%
% SYNTAX:
%  cf = cf_Ratio(t,cfX1,pdfX2,a,b)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  cfX1   - function handle @(t) of the chcaracteristic function of the
%           random variable X1.   
%  pdfX2  - function handle @(x) of the probability density function of the
%           random variable X2. If empty, the default value is PDF of the
%           rectangular distribution over interval the interval [-1,1],
%           i.e. pdf = @(x) 1/(b-a), where b = 1 and a = -1.
%  a      - lower limit of the support interval [a,b] of the distribution
%           of the random variable X2. If empty,the default value a = -1;
%  b      - upper limit of the support interval [a,b] of the distribution
%           of the random variable X2. If empty, the default value b = 1;
%  tol    - relative tolerance of integration. If empty, the default value
%           is tol = 1e-6; 
%
% OUTPUT:
%  cf     - characteristic function of the product of two independent
%           random variables, i.e. CF of Y = X1/X2, evaluated at the vector
%           of real values t.
%
% EXAMPLE 1:
%  % CF of a ratio of independent RVs with central symmetric triangular
%  % and rectangular distributions
%  t = linspace(-20,20,2^10+1);
%  cfX1 = @(t) cfS_Triangular(t);
%  cf   = cf_Ratio(t,cfX1);
%  figure; plot(t,real(cf),t,imag(cf));grid on
%  title('CF of a ratio of independent RVs ')
%
% EXAMPLE 2:
%  % CF of a ratio of independent RVs with chi-squared (df = 5) and
%  % rectangular distributions 
%  t = linspace(-5,5,2^10+1);
%  cfX1 = @(t) cfX_ChiSquare(t,5);
%  cf   = cf_Ratio(t,cfX1);
%  figure; plot(t,real(cf),t,imag(cf));grid on
%  title('CF of a ratio of independent RVs ')
%
% EXAMPLE 3:
%  % CF of a ratio of independent RVs with chi-squared (df = 5) and
%  % standard normal distributions 
%  t = linspace(-5,5,2^10+1);
%  cfX1  = @(t) cfX_ChiSquare(t,5);
%  pdfX2 = @(x) normpdf(x);
%  a     = -8;
%  b     =  8;
%  cf    = cf_Ratio(t,cfX1,pdfX2,a,b);
%  figure; plot(t,real(cf),t,imag(cf));grid on
%  title('CF of a ratio of independent RVs ')
%
% EXAMPLE 4:
%  % PDF/CDF of a ratio of independent RVs with symmetric standard
%  % triangular and rectangular distributions
%  cfX1  = @(t) cfS_Triangular(t);
%  pdfX2 = @(x) 0.5;
%  a     = -1;
%  b     = 1;
%  tol   = 1e-4;
%  cf    = @(t) cf_Ratio(t,cfX1,pdfX2,a,b,tol);
%  x     = linspace(-30,30,1001);
%  prob  = [0.9 0.95 0.99];
%  clear options;
%  options.N = 2^10;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 5:
%  % PDF/CDF of a ratio of independent RVs with chi-squared (df = 5) and
%  % standard normal distributions 
%  cfX1  = @(t) cfX_ChiSquare(t,5);
%  pdfX2 = @(x) normpdf(x);
%  a     = -8;
%  b     = 8;
%  tol   = 1e-4;
%  cf    = @(t) cf_Ratio(t,cfX1,pdfX2,a,b,tol);
%  prob  = [0.9 0.95 0.99];
%  x = linspace(-100,100,1001);
%  clear options;
%  options.N = 2^12;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 6:
%  % PDF/CDF of a ratio of independent RVs with standard normal
%  % distribution and chi-square distribution with 5 degrees of freedom 
%  cfX1  = @(t) cfS_Gaussian(t);
%  pdfX2 = @(x) chi2pdf(x,5);
%  a     = 0;
%  b     = Inf;
%  tol   = 1e-4;
%  cf    = @(t) cf_Ratio(t,cfX1,pdfX2,a,b,tol);
%  x     = linspace(-3,3,501);
%  prob  = [0.9 0.95 0.99];
%  clear options;
%  options.N = 2^10;
%  result = cf2DistGP(cf,x,prob,options)

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Feb-2018 09:10:33

%% ALGORITHM
% cf = cf_Ratio(t,cfX1,pdfX2,a,b)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 6);

if nargin < 6, tol   = []; end
if nargin < 5, b     = []; end
if nargin < 4, a     = []; end
if nargin < 3, pdfX2 = []; end

if isempty(b),   b   =  1; end
if isempty(a),   a   = -b; end
if isempty(tol), tol = 1e-6; end

if isempty(pdfX2) 
    pdfX2 = @(x) 1/(b-a);
end

%% 
szt = size(t);
t   = t(:);

fun = @(x) cfX1(t/x) .* pdfX2(x);
cf  = integral(fun,a,b,'ArrayValued',true,'RelTol',tol);

cf  = reshape(cf,szt);

end