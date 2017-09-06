function cf = cf_PDF(t,pdfFun,pdfMin,pdfMax,nPts)
%cf_PDF Computes the characteristic function of the continuos
% distribution defined by its PDF function, pdfFun = @(x) pdf(x). 
%
% DEFINITION:
%  cfX_PDF is based on the standard integral representation of the
%  characteristic function of the continuous distribution defined by its 
%  PDF (here PDF is represented by the function handle pdfFun(x) defined 
%  for x >= 0). Then, 
%    CF(t) = Integral_0^inf exp(i*t*x) * pdfFun(x) dx.
%
% SYNTAX:
%  cf = cf_PDF(t,pdfFun,pdfMin,pdfMax,nPts)
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
% EXAMPLE1 (CF of the Exponential distribution with lambda = 1)
%  pdfFun = @(x) exp(-x);
%  pdfMin = 0;
%  pdfMax = 100;
%  nPts   = 25;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cf_PDF(t,pdfFun,pdfMin,pdfMax,nPts);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Exponential distribution')
%
% EXAMPLE2 (CF of the LogNormal distribution with mu = 0, sigma = 1)
%  mu     = 0;
%  sigma  = 1;
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  pdfMin = 0;
%  pdfMax = 100;
%  nPts   = 100;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cf_PDF(t,pdfFun,pdfMin,pdfMax,nPts);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogNormal distribution')
%
% EXAMPLE3 (CF of the Weibull distribution with a = 1.5, and large b > 1)
%  a      = 1.5;
%  b      = 3.5;
%  pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  pdfMin = 0;
%  pdfMax = 100;
%  t  = linspace(-10,10,2^10+1)';
%  cf = cf_PDF(t,pdfFun,pdfMin,pdfMax);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%
% REFERENCES:
% [1] BAKHVALOV, N.S., VASILEVA, L.G. Evaluation of the integrals of
%     oscillating functions by interpolation at nodes of Gaussian
%     quadratures. USSR Computational Mathematics and Mathematical Physics,
%     1968, 8(1): 241-249.
% [2] EVANS, G.A.; WEBSTER, J.R. A comparison of some methods for the
%     evaluation of highly oscillatory integrals. Journal of Computational
%     and Applied Mathematics, 1999, 112(1): 55-69.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 06-Sep-2017 17:28:00

%% ALGORITHM
%cf = cf_PDF(t,pdfFun,pdfMin,pdfMax,nPts)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, nPts = []; end
if nargin < 4, pdfMax = []; end
if nargin < 3, pdfMax = []; end
if nargin < 2, pdfFun = []; end

if isempty(nPts)
    nPts = 100;
end

if isempty(pdfFun)
    pdfFun = @(x) exp(-x);
end

if isempty(pdfMax)
    pdfMax = Inf;
end

if isempty(pdfMin)
    pdfMin = 0;
end

%%

[x,w]=legzo(nPts+1);
P = zeros(nPts+1);
P(1,:) = 1;
p1 = 1;
p2 = 0;
for jj = 1:nPts
    p3 = p2;
    p2 = p1;
    p1 = ((2*jj-1).*x.*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
    P(jj+1,:) = p1;
end

m = (pdfMax-pdfMin)/2;
c = (pdfMax+pdfMin)/2;
f = pdfFun(m*x+c);
fw = f.*w;
Pfw = sum(bsxfun(@times,P,fw),2);

szt = size(t);
t   = t(:)';
id    = t < 0;
t(id) = -t(id); 
nt  = length(t);
K = m*exp(1i*t*c);
t = m*t;
B = zeros(nPts+1,nt);
for k = 0:nPts
    B(k+1,:) = 1i^k * (2*k+1) * (pi ./ (2*t)).^0.5 .* besselj(k+0.5,t);
end

cf       = K .* sum(bsxfun(@times,B,Pfw));
cf(t==0) = 1;
cf(id)   = conj(cf(id));
cf       = reshape(cf,szt);


%% Function legzo
function [x,w]=legzo(n, a, b)
%       =========================================================
%       Purpose : Compute the zeros of Legendre polynomial Pn(x)
%                 in the interval [a,b], and the corresponding
%                 weighting coefficients for Gauss-Legendre
%                 integration
%       Input :   n    --- Order of the Legendre polynomial
%                 a    --- Lower boundary (optional)
%                 b    --- Upper boundary (optional)
%       Output:   x(n) --- Zeros of the Legendre polynomial
%                 w(n) --- Corresponding weighting coefficients
%       =========================================================
if nargin == 1
    a = -1;
    b =  1;
end
x = zeros(1, n);
w = zeros(1, n);
m = (n+1)/2;
h = b-a;

for ii=1:m
    z = cos(pi*(ii-.25)/(n+.5)); % Initial estimate.
    z1 = z+1;
    while abs(z-z1)>eps
        p1 = 1;
        p2 = 0;
        for jj = 1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
        end
        pp = n*(z*p1-p2)/(z^2-1); % The L.P. derivative.
        z1 = z;
        z = z1-p1/pp;
    end
    x(ii) = z; % Build up the abscissas.
    x(n+1-ii) = -z;
    w(ii) = h/((1-z^2)*(pp^2)); % Build up the weights.
    w(n+1-ii) = w(ii);
end

if a ~= -1 || b ~= 1
    x = (x+1)*(h/2) + a;
end