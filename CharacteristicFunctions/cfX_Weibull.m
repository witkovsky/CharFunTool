function cf = cfX_Weibull(t,alpha,beta,tol)
%cfX_Weibull(t,alpha,beta) Computes the characteristic function cf(t)
% of the Weibull distribution with the parameters alpha(scale parameter,
% elsewhere denoted also as lambda, alpha > 0) and beta (shape parameter,
% elsewhere denoted as k, beta > 0), for real (vector) argument t, i.e.
%  cf(t) = cfX_Weibull(t,alpha,beta);
%
% REMARK (Caution Notice)
% This implementation of the algorithm is experimental. It DOES NOT WORK
% PROPERLY for all parameters!!!! 
%
% The closed-form analytic expression of the characteristic function of the
% Weibull distribution is unknown. The series expansion of the CF, valid
% for abs(t)<1, is known and given, e.g. in WIKIPEDIA. 
% 
% Thus, for abs(t)<1 the CF is evalueted from its series expansion, and for
% abs(t)>=1 from its definition, as suggested in [3]. 
% For more details see also WIKIPEDIA:  
% https://en.wikipedia.org/wiki/Weibull_distribution
%  
% SYNTAX:
%  cf = cfX_Weibull(t,alpha,beta)
%  cf = cfX_Weibull(t,alpha,beta,tol)
%
% EXAMPLE1 (CF of the Weibull distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta = 1;
%  t = linspace(-20,20,2^10+1)';
%  cf = cfX_Weibull(t,alpha,beta);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% EXAMPLE2 (CDF/PDF of the  Weibull distribution with alpha=1, beta=1)
%  alpha = 1;
%  beta = 1;
%  x = linspace(0,7,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.xMin = 0;
%  options.N = 2^10;
%  cf = @(t) cfX_Weibull(t,alpha,beta);
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Weibull distribution)
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
id       = abs(t/alpha) <= 0.99;
cf(id)   = FoxWrightPsi(1,1/beta,[],[],-1i*alpha*t(id));
id       = abs(t/alpha) > 0.99;
cf(id)   = cfX_PDF(t(id),pdfFun,tol);
cf       = reshape(cf,szt);
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
%% FoxWrightPsi
function [pPSIq,c] = FoxWrightPsi(a,A,b,B,z,n)
%FoxWrightPsi  The Fox-Wright Psi function pPSIq(a,A,b,B,z)
%  The Fox-Wright Psi function pPSIq(a,A,b,B,z) is defined for complex
%  arguments ai, Ai, for i = 1,...,p, and bj, Bj for j = 1,...,q, and z.
%  With a = [a1, a2, …, ap], A = [A1, A2, …, Ap], b = [b1, b2, …, bq], and
%  B = [B1, B2, …, Bq], the The Fox-Wright Psi function of order p, q is
%  defined as
%   pPSIq(a,A,b,B,z) = \sum_{n=0}^\infty 
%      \frac{\Gamma( a_1 + A_1 n )\cdots\Gamma( a_p + A_p n )}
%           {\Gamma( b_1 + B_1 n )\cdots\Gamma( b_q + B_q n )} 
%      \frac {z^n} {n!}.
%
%  This is a very simple implementation of the function pPSIq(a,A,b,B,z)
%  which is computed directly from the definition - as a truncated power
%  series. The sum is truncated to n terms (by default n = 500).
% 
%  An important property is the convergence criteria of the  functions
%  depending on the values of p and q (the radius of convergence of a
%  series as a variable of z).
% 
%  Numerical stability of the present algorithm could be strongly affected
%  at z values close to the border of the convergence region.
%
% SYNTAX
%    pPSIq = FoxWrightPsi(a,A,b,B,z)
%    pPSIq = FoxWrightPsi(a,A,b,B,z,n)
%
% EXAMPLE1
%  a = 3; A = 1.5;
%  b = 1.5; B = 2.5;
%  z = 1i*(0:0.05:1)';
%  pPSIq = FoxWrightPsi(a,A,b,B,z)
%
% EXAMPLE2 (Characteristic function of Weibull distribution)
%  lambda = 1;
%  k = 5
%  a = 1; A = 1/k;
%  b = []; B = [];
%  t = linspace(-20,20,501)';
%  z = lambda*1i*t;
%  cf = FoxWrightPsi(a,A,b,B,z)
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the Weibull(1,5) distribution')
%  xlabel('t')
%  ylabel('CF')

% REFERENCES
% https://en.wikipedia.org/wiki/Fox%E2%80%93Wright_function

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Aug-2015 11:00:07

%% CHECK THE INPUT PARAMETERS
narginchk(5, 6);
if nargin < 6, n = []; end

if isempty(n)
    n = 500;
end

%% ALGORITHM
sz = size(z);
z  = z(:);
a  = a(:)';
A  = A(:)';
b  = b(:)';
B  = B(:)';
k  = (0:n)';

c  = sum(gammalog(bsxfun(@plus,a,bsxfun(@times,A,k))),2) - ...
     sum(gammalog(bsxfun(@plus,b,bsxfun(@times,B,k))),2);

pPSIq  = sum(exp(bsxfun(@plus,c-gammalog(k+1),bsxfun(@times,k,log(z')))));

pPSIq(z==0) = complex(exp(c(1)));
pPSIq  = reshape(pPSIq,sz);

end
%%
function [f] = gammalog(z)
% GAMMALOG  Natural Log of the Gamma function valid in the entire complex
%           plane. This routine uses an excellent Lanczos series
%           approximation for the complex ln(Gamma) function.
%
%usage: [f] = gammalog(z)
%             z may be complex and of any size.
%             Also  n! = prod(1:n) = exp(gammalog(n+1))
%
%tested under version 5.3.1
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%
%see also:   GAMMA GAMMALN GAMMAINC PSI
%see also:   mhelp GAMMA
%see also:   mhelp lnGAMMA

%Paul Godfrey
%pgodfrey@conexant.com
%07-13-01

siz = size(z);
z   = z(:);
zz  = z;

%f = 0.*z; % reserve space in advance

p = find(real(z)<0);
if ~isempty(p)
    z(p) = -z(p);
end

%Lanczos approximation for the complex plane

g=607/128; % best results when 4<=g<=5

c = [0.99999999999999709182;
    57.156235665862923517;
    -59.597960355475491248;
    14.136097974741747174;
    -0.49191381609762019978;
    0.33994649984811888699e-4;
    0.46523628927048575665e-4;
    -0.98374475304879564677e-4;
    0.15808870322491248884e-3;
    -0.21026444172410488319e-3;
    0.21743961811521264320e-3;
    -0.16431810653676389022e-3;
    0.84418223983852743293e-4;
    -0.26190838401581408670e-4;
    0.36899182659531622704e-5];

s = 0;
for k = size(c,1):-1:2
    s = s + c(k)./(z+(k-2));
end

zg   = z+g-0.5;
s2pi = 0.9189385332046727417803297;

f = (s2pi + log(c(1)+s)) - zg + (z-0.5).*log(zg);

f(z==1 | z==2) = 0.0;

if ~isempty(p)
    lpi  = 1.14472988584940017414342735 + 1i*pi;
    f(p) = lpi-log(zz(p))-f(p)-log(sin(pi*zz(p)));
end

p = find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
    f(p) = Inf;
end

f = reshape(f,siz);
end