function [cf,coefs] = cf_PdfBV(t,pdfFun,A,B,nPts)
%cf_PdfBV 
%  Computes the characteristic function of the continuos distribution from
%  its probability density function (PDF) by using the BAKHVALOV-VASILEVA
%  method for computing the Fourier Integral over the interval (A,B).
%  pdfFun is function handle for computing the PDF at arbitrary x in (A,B).
%
%  Here, the interval (A,B) denotes the support of the distribution or
%  the given truncation points with A < B. In the current version, A and B
%  must be finite. 
%
%  cf_PdfBV is evaluated from the standard integral representation of the
%  characteristic function of the continuous distribution defined by its
%  PDF (here represented by the function handle pdfFun), i.e.
%    CF(t) = Integral_A^B exp(i*t*x) * pdfFun(x) dx,
%  using the algorithms suggested for computing the oscillatory Fourier
%  integrals, based on approximation of the PDF function by the polynomials
%  and observation that Fourier transform of the Legendre polynomials is
%  related to the Belssel J functions. For more details see Evans and
%  Webster (1999). 
%
%  pdfFun must be a function handle. In the current version, A and B must
%  be finite. 
%
% REMARK:
%  cf_PdfBV is suggested for situations when the PDF can be well approximated
%  by the nth order polynomial over known (given) support interval [A,B],
%  as e.g. the uniform distribution PDF = 1 over [0,1].
%  Otherwise the computed result could be misleading!
%
% SYNTAX:
%  [cf,coefs] = cf_PdfBV(t,pdfFun,A,B,nPts)
%
% INPUTS:
%  t      - real vector, where the characteristic function CF(t) will
%           be evaluated.
%  pdfFun - function handle used as the PDF function with the argument x.
%  A      - finite minimum value of the distribution support or the lower
%           truncation point. If the true minimum is -Inf, than A should be
%           set as a reasonable finite approximation of the minimum value
%           of the support. Default value is A = -100.
%  B      - finite maximum value of the distribution support or the upper
%           truncation point. If the true maximum is Inf, than B should be
%           set as a reasonable finite approximation of the maximum value
%           of the support. Default value is B = 100.
%  nPts   - Order of Legendre polynomial approximation. Default value is
%           nPts = 100.
%
% OUTPUT:
%  cf     - (complex) vector of the characteristic function values,
%            evalated at the required t, i.e. CF(t).
%  coefs  - vector of the polynomial expansion coefficiens of the PDF.
%           This allows to check the quality of the polynomial
%           approaximation over the interval [A,B].
%
% EXAMPLE 1 (CF of the Uniform distribution on the interval [0,1])
%  pdfFun = @(x) ones(size(x));
%  A = 0;
%  B = 1;
%  nPts = 1;
%  t  = linspace(-50,50,2^10+1)';
%  cf = cf_PdfBV(t,pdfFun,A,B,nPts);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Uniform distribution')
%
% EXAMPLE 2 (CF of the truncated Normal distribution on the interval [0,5])
%  pdfFun = @(x) exp(-x.^2/2)/sqrt(2*pi);
%  A = 0;
%  B = 5;
%  nPts = 30;
%  t  = linspace(-20,20,2^10+1)';
%  [cf,coefs] = cf_PdfBV(t,pdfFun,A,B,nPts);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the truncated Normal distribution')
%
% EXAMPLE 3 (PDF/CDF of the truncated Normal distribution)
%  pdfFun = @(x) exp(-x.^2/2)/sqrt(2*pi);
%  A = -2;
%  B = 5;
%  t  = linspace(-20,20,2^10+1)';
%  cf = @(t) cf_PdfBV(t,pdfFun,A,B);
%  x  = linspace(-2,5);
%  prob = [0.9, 0.95 0.99];
%  clear options
%  options.xMin = -2;
%  options.xMax = 5;
%  figure
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 4 (CF of the truncated Exponential distribution with lambda = 1)
%  pdfFun = @(x) exp(-x);
%  A = 0;
%  B = 15;
%  nPts   = 25;
%  t  = linspace(-20,20,2^10+1)';
%  [cf,coefs] = cf_PdfBV(t,pdfFun,A,B,nPts);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the truncated Exponential distribution')
%
% EXAMPLE 5 (CF of the LogNormal distribution with mu = 0, sigma = 1)
%  mu     = 0;
%  sigma  = 1;
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  A = 1e-8;
%  B = 100;
%  t  = linspace(-20,20,2^10+1)';
%  [cf,coefs] = cf_PdfBV(t,pdfFun,A,B);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogNormal distribution')
%
% EXAMPLE 6 (CF of the Weibull distribution with a = 1.5, and large b > 1)
%  a      = 1.5;
%  b      = 3.5;
%  pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  A = 1e-8;
%  B = 100;
%  t  = linspace(-10,10,2^10+1)';
%  [cf,coefs] = cf_PdfBV(t,pdfFun,A,B);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Weibull distribution')
%   
% REFERENCES: 
% [1] BAKHVALOV, N.S., VASILEVA, L.G. Evaluation of the integrals of
%     oscillating functions by interpolation at nodes of Gaussian
%     quadratures. USSR Computational Mathematics and Mathematical Physics,
%     1968, 8(1): 241-249.
% [2] PATTERSON, T. N. L. On high precision methods for the evaluation of
%     Fourier integrals with finite and infinite limits. Numerische
%     Mathematik, 1976, 27(1): 41-52.  
% [3] EVANS, G.A.; WEBSTER, J.R. A comparison of some methods for the
%     evaluation of highly oscillatory integrals. Journal of Computational
%     and Applied Mathematics, 1999, 112(1): 55-69.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Sep-2018 18:04:33

%% ALGORITHM CALL
%[cf,coefs] = cf_PdfBV(t,pdfFun,A,B,nPts)

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
    pdfFun = @(x) exp(-x.^2/2)/sqrt(2*pi);
end

if isempty(B)
    B = 100;
end

if isempty(A)
    A = -100;
end

%% ALGORITHM

[cf,coefs] = FourierIntegral_BV(t,pdfFun,A,B,nPts);

pr = (abs(FourierIntegral_BV(-1e-8,pdfFun,A,B,nPts)) + ...
    abs(FourierIntegral_BV(1e-8,pdfFun,A,B,nPts)))/2;

if pr < 1
    cf = cf / pr;
end

cf(t==0) = 1;

end