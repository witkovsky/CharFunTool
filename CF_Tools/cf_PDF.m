function [cf,coefs,method] = cf_PDF(t,pdfFun,A,B,method,nPts)
%cf_PDF Computes the characteristic function of the continuos
%  distribution defined by its PDF function.
%
%  cf_PDF is evaluated from the standard integral representation of the
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
%  cf_PDF is suggested for situations when the PDF can be well approximated
%  by the nth order polynomial over known (given) support interval [A,B],
%  as e.g. the uniform distribution PDF = 1 over [0,1].
%  Otherwise the computed result could be misleading!
%
% SYNTAX:
%  [cf,coefs,method] = cf_PDF(t,pdfFun,A,B,method,nPts)
%
% INPUTS:
%  t      - real vector, where the characteristic function CF(t) will
%           be evaluated.
%  pdfFun - function handle used as the PDF function with the argument x.
%  A      - finite minimum value of the distribution support. If the true
%           minimum is -Inf, than A should be set as a reasonable finite
%           approximation of the minimum value of the support. Default
%           value is A = -100.
%  B      - finite maximum value of the distribution support. If the true
%           maximum is Inf, than B should be set as a reasonable finite
%           approximation of the maximum value of the support. Default
%           value is B = 100.
%  method - select method for computing the required Fourier Integral.
%           Currently, there are two possible algorithms: (i)the
%           BAKHVALOV-VASILEVA method (method = 'bv') and (ii) the
%           PATTERSON method (method = 'pat'). Deafault value is  method =
%           'pat'. 
%  nPts   - Order of Legendre polynomial approximation. Default value is
%           nPts = 100.
%
% OUTPUT:
%  cf     - (complex) vector of the characteristic function values,
%            evalated at the required t, i.e. CF(t).
%  coefs  - vector of the polynomial expansion coefficiens of the PDF.
%           This allows to check the quality of the polynomial
%           approaximation over the interval [A,B].
%  method - selected method for computing the required Fourier Integral.
%
% EXAMPLE1 (CF of the Uniform distribution on the interval [0,1])
%  pdfFun = @(x) 1;
%  A = 0;
%  B = 1;
%  method = 'bv';
%  nPts   = 1;
%  t  = linspace(-50,50,2^10+1)';
%  cf = cf_PDF(t,pdfFun,A,B,method,nPts);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Uniform distribution')
%
% EXAMPLE2 (CF of the Normal distribution on the interval [-8,8])
%  % !!PDF cannot be approximated by the 20th degree Legendre expansion!
%  pdfFun = @(x) exp(-x.^2/2)/sqrt(2*pi);
%  A = -8;
%  B = 8;
%  method = 'bv';
%  nPts   = 20;
%  t  = linspace(-10,10,2^10+1)';
%  [cf,coefs] = cf_PDF(t,pdfFun,A,B,method,nPts);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Normal distribution')
%
% EXAMPLE3 (CF of the Exponential distribution with lambda = 1)
%  pdfFun = @(x) exp(-x);
%  A = 0;
%  B = 100;
%  method = 'patterson';
%  nPts   = 25;
%  t  = linspace(-20,20,2^10+1)';
%  [cf,coefs] = cf_PDF(t,pdfFun,A,B,method,nPts);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the Exponential distribution')
%
% EXAMPLE4 (CF of the LogNormal distribution with mu = 0, sigma = 1)
%  mu     = 0;
%  sigma  = 1;
%  pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
%  A = 1e-8;
%  B = 100;
%  t  = linspace(-20,20,2^10+1)';
%  [cf,coefs] = cf_PDF(t,pdfFun,A,B);
%  plot(coefs,'o-');grid
%  title('Expansion coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the LogNormal distribution')
%
% EXAMPLE5 (CF of the Weibull distribution with a = 1.5, and large b > 1)
%  a      = 1.5;
%  b      = 3.5;
%  pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  A = 1e-8;
%  B = 100;
%  t  = linspace(-10,10,2^10+1)';
%  [cf,coefs] = cf_PDF(t,pdfFun,A,B);
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
% Ver.: 8-Sep-2017 15:27:01

%% ALGORITHM CALL
%[cf,coefs,method] = cf_PDF(t,pdfFun,A,B,method,nPts)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, nPts = []; end
if nargin < 5, method = []; end
if nargin < 4, B = []; end
if nargin < 3, A = []; end
if nargin < 2, pdfFun = []; end

if isempty(nPts)
    nPts = 100;
end

if isempty(method)
    method = 'patterson';
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

switch lower(method)
    case {'patterson','pat','p'}
        [cf,coefs] = FourierIntegral_P(t,pdfFun,A,B,nPts);
    case {'bakhvalov-vasileva','bv','b'}
        [cf,coefs] = FourierIntegral_BV(t,pdfFun,A,B,nPts);
    otherwise
        [cf,coefs] = FourierIntegral_P(t,pdfFun,A,B,nPts);
end
cf(t==0) = 1;
end