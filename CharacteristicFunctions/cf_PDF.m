function cf = cf_PDF(t,pdfFun,A,B,nPts)
%cf_PDF Computes the characteristic function of the continuos
% distribution defined by its PDF function.
%
% DEFINITION:
%  cf_PDF is evaluated from the standard integral representation of the
%  characteristic function of the continuous distribution defined by its
%  PDF (here represented by the function handle pdfFun), i.e.
%    CF(t) = Integral_A^B exp(i*t*x) * pdfFun(x) dx,
%  using the efficient BAKHVALOV-VASILEVA algorithm suggested for computing
%  the oscillatory Fourier integrals, which is based on approximation of
%  the PDF function by the Legendre polynomials and observation that
%  Fourier transform of the Legendre polynomials is related to the Belssel
%  J functions.  
%
% For more details see Evans and Webster (1999).
%
% SYNTAX:
%  cf = cf_PDF(t,pdfFun,A,B,nPts)
%
% INPUTS:
%  t      - real vector, where the characteristic function CF(t) will
%           be evaluated,
%  pdfFun - function handle used as the PDF function with the argument x.
%  A      - finite minimum value of the distribution support. If the true
%           minimum is -Inf, than A should be set as a reasonable finite
%           approximation of the minimum value of the support. 
%  B      - finite maximum value of the distribution support. If the true
%           maximum is Inf, than B should be set as a reasonable finite
%           approximation of the maximum value of the support. 
%  nPts   - Order of Legendre polynomial approximation. Default value is
%           nPts = 100.
%
% OUTPUT:
%  cf     - (complex) vector of the characteristic function values,
%            evalated at the required t, i.e. CF(t).
%
% EXAMPLE1 (CF of the Exponential distribution with lambda = 1)
%  pdfFun = @(x) exp(-x);
%  A = 0;
%  B = 100;
%  nPts   = 25;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cf_PDF(t,pdfFun,A,B,nPts);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Exponential distribution')
%
% EXAMPLE2 (CF of the LogNormal distribution with mu = 0, sigma = 1)
%  mu     = 0;
%  sigma  = 1;
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  A = 0;
%  B = 100;
%  nPts   = 100;
%  t  = linspace(-20,20,2^10+1)';
%  cf = cf_PDF(t,pdfFun,A,B,nPts);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogNormal distribution')
%
% EXAMPLE3 (CF of the Weibull distribution with a = 1.5, and large b > 1)
%  a      = 1.5;
%  b      = 3.5;
%  pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  A = 0;
%  B = 100;
%  t  = linspace(-10,10,2^10+1)';
%  cf = cf_PDF(t,pdfFun,A,B);
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
% Ver.: 07-Sep-2017 17:28:00

%% ALGORITHM
%cf = cf_PDF(t,pdfFun,A,B,nPts)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, nPts = []; end
if nargin < 4, B = []; end
if nargin < 3, A = []; end
if nargin < 2, pdfFun = []; end

if isempty(nPts)
    nPts = 100;
end

if isempty(pdfFun)
    pdfFun = @(x) exp(-x);
end

if isempty(B)
    B = Inf;
end

if isempty(A)
    A = -Inf;
end

% % Set proper rules for dealing with infinite A and/or B
% if ~isfinite(A)
%     A = -1e2;
% end
% 
% if ~isfinite(B)
%     B = 1e2;
% end

%% BAKHVALOV-VASILEVA ALGORITHM

% Nodes and weights of the n-th order Gauss-Legendre quadrature on [-1,1]
[x,w]  = LegendrePts(nPts+1);
P      = zeros(nPts+1);
P(1,:) = 1;

% The Legendre polynomials evaluated at x
p1 = 1;
p2 = 0;
for k = 1:nPts
    p3 = p2;
    p2 = p1;
    p1 = ((2*k-1) .* x .* p2 - (k-1) * p3) / k;
    P(k+1,:) = p1;
end

% PDF function (weighted, shifted and scaled) evaluated at x
scale   = (B - A) / 2;
shift   = (B + A) / 2;
f   = pdfFun(scale*x + shift);
fw  = f.*w;
Pfw = sum(bsxfun(@times,P,fw),2);

% Bessel J functions evaluated at required values
szt   = size(t);
t     = t(:)';
id    = t < 0;
t(id) = -t(id);
nt    = length(t);
K     = scale * exp(1i*t*shift);
t     = scale * t;
B     = zeros(nPts+1,nt);
for k = 0:nPts
    B(k+1,:) = 1i^k * (2*k+1) * (pi ./ (2*t)).^0.5 .* besselj(k+0.5,t);
end

% Characteristic function evaluated at required values t
cf       = K .* sum(bsxfun(@times,B,Pfw));
cf(t==0) = 1;
cf(id)   = conj(cf(id));
cf       = reshape(cf,szt);
end