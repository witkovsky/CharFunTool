function [FI,coefs] = FourierIntegral_BV(t,fun,A,B,nPts)
% FourierIntegral_BV approximates the Fourier integral of function FUN from
%   A to B using the BAKHVALOV-VASILEVA method suggested for computing the
%   oscillatory Fourier integrals based on approximation of the integrand
%   function by the FOURIER-LEGENDRE SERIES EXPANSION, and observation that
%   Fourier transform of the Legendre polynomials is related to the Belssel
%   J functions. For more details see Bakhlanov and Vasileva (1968) and
%   Evans and Webster (1999).  
%
%   The FOURIER INTEGRAL, evaluated at given vector omega, is defined by
%     FI(omega) = Integral_A^B  fun(x) * exp(i*omega*x) dx.
%
%   FUN must be a function handle. In the current version, A and B must be
%   finite. Function Y = FUN(X) must accept a vector argument X and return
%   a vector result Y.
%
% REMARK:
%  FourierIntegral_BV is suggested for situations when the PDF can be well
%  approximated by the nth order polynomial over known (given) support
%  interval [A,B], as e.g. the uniform distribution PDF = 1 over [0,1].
%  Otherwise the computed result could be misleading!
%
% SYNTAX:
%  [FI,coefs] = FourierIntegral_BV(t,fun,A,B,nPts)
%
% INPUTS:
%  t      - real vector, where the characteristic function CF(t) will
%           be evaluated,
%  fun    - function handle used as the PDF function with the argument x.
%  A      - finite minimum value of the distribution support. If the true
%           minimum is -Inf, than A should be set as a reasonable finite
%           approximation of the minimum value of the support. Default
%           value is A = -1.
%  B      - finite maximum value of the distribution support. If the true
%           maximum is Inf, than B should be set as a reasonable finite
%           approximation of the maximum value of the support.  Default
%           value is B = 1.
%  nPts   - Order of Legendre polynomial approximation. Default value is
%           nPts = 100.
%
% OUTPUT:
%  FI     - (complex) vector of Fourier integral values evalated at the
%           required t, i.e. FI(t).
%  coefs  - vector of the polynomial expansion coefficiens of the PDF.
%           This allows to check the quality of the polynomial
%           approaximation over the interval [A,B].
%
% EXAMPLE1 (CF of the Uniform distribution on the interval [0,1])
%  fun = @(x) 1;
%  A = 0;
%  B = 1;
%  nPts   = 1;
%  t  = linspace(-50,50,2^10+1)';
%  cf = FourierIntegral_BV(t,fun,A,B,nPts);
%  cf(t==0) = 1;
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Uniform distribution')
%
% EXAMPLE2 (CF of the Normal distribution on the interval [-8,8])
%  % !!PDF cannot be approximated by the 20th degree Legendre expansion!
%  fun = @(x) exp(-x.^2/2)/sqrt(2*pi);
%  A = -8;
%  B = 8;
%  nPts   = 20;
%  t  = linspace(-10,10,2^10+1)';
%  [cf,coefs] = FourierIntegral_BV(t,fun,A,B,nPts);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  cf(t==0) = 1;
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Normal distribution')
%
% EXAMPLE3 (CF of the Exponential distribution with lambda = 1)
%  fun = @(x) exp(-x);
%  A = 0;
%  B = 100;
%  nPts   = 25;
%  t  = linspace(-20,20,2^10+1)';
%  [cf,coefs] = FourierIntegral_BV(t,fun,A,B,nPts);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  cf(t==0) = 1;
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Exponential distribution')
%
% EXAMPLE4 (CF of the LogNormal distribution with mu = 0, sigma = 1)
%  mu     = 0;
%  sigma  = 1;
%  fun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  A = 0;
%  B = 100;
%  nPts   = 100;
%  t  = linspace(-20,20,2^10+1)';
%  [cf,coefs] = FourierIntegral_BV(t,fun,A,B,nPts);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  cf(t==0) = 1;
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogNormal distribution')
%
% EXAMPLE5 (CF of the Weibull distribution with a = 1.5, and large b > 1)
%  a      = 1.5;
%  b      = 3.5;
%  fun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  A = 0;
%  B = 100;
%  t  = linspace(-10,10,2^10+1)';
%  [cf,coefs] = FourierIntegral_BV(t,fun,A,B);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  cf(t==0) = 1;
%  figure
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

%% ALGORITHM CALL
%[FI,coefs] = FourierIntegral_BV(t,fun,A,B,nPts)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, nPts = []; end
if nargin < 4, B = []; end
if nargin < 3, A = []; end
if nargin < 2, fun = []; end

if isempty(nPts)
    nPts = 100;
end

if isempty(fun)
    fun = @(x) exp(-x);
end

if isempty(B)
    B = Inf;
end

if isempty(A)
    A = -Inf;
end

%% BAKHVALOV-VASILEVA ALGORITHM

% Nodes and weights of the n-th order Gauss-Legendre quadrature on [-1,1]
[x,w]   = LegendrePoints(nPts+1);
P       = zeros(nPts+1);
P(1,:)  = 1;

% The Legendre polynomials evaluated at nodes x
p1 = 1;
p2 = 0;
for k   = 1:nPts
    p3  = p2;
    p2  = p1;
    p1  = ((2*k-1) .* x .* p2 - (k-1) * p3) / k;
    P(k+1,:) = p1;
end

% Function fun (weighted, shifted and scaled) evaluated at nodes x
scale  = (B - A) / 2;
shift  = (B + A) / 2;
F      = fun(scale*x + shift);
coefs  = F.*w;
coefs  = sum(bsxfun(@times,P,coefs),2);

% Bessel J functions evaluated at required values t
szt    = size(t);
t      = t(:)';
id     = t < 0;
t(id)  = -t(id);
nt     = length(t);
K      = scale * exp(1i*t*shift);
t      = scale * t;
B      = zeros(nPts+1,nt);
for k  = 0:nPts
    B(k+1,:) = 1i^k * (2*k+1) * besselj(k+0.5,t);
end

% Fourier integral evaluated at required values t
FI     = K .* sum(bsxfun(@times,B,coefs)) ./ sqrt(2*t/pi);
FI(id) = conj(FI(id));
FI     = reshape(FI,szt);
end