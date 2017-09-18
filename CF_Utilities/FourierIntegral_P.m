function [FI,coefs] = FourierIntegral_P(omega,fun,A,B,nPts)
% FourierIntegral_P approximates the Fourier integral of function FUN
%   from A to B by using the PATTERSON method (extended Clenshaw-Curtis
%   quadrature) suggested in Patterson (1976), based on the Bakhlanov
%   and Vasileva (1968) method. For more details see also Evans and Webster
%   (1999).
%
%   The Fourier integral, evaluated at given vector omega, is defined by
%     FI(omega) = Integral_A^B  fun(x) * exp(i*omega*x) dx.
%
%   FUN must be a function handle. In the current version, A and B must be
%   finite. Function Y = FUN(X) must accept a vector argument X and return
%   a vector result Y.
%
% REMARK:
%  FourierIntegral_P is suggested for situations when the function FUN
%  can be well approximated by the nth order polynomial over known finite
%  support interval [A,B], as e.g. the uniform distribution PDF = 1 over
%  [0,1]. Otherwise the computed result could be misleading!
%  
% SYNTAX:
%  [FI,coefs] = FourierIntegral_P(omega,fun,A,B,nPts)
%
% INPUTS:
%  omega  - real vector, where the Fourier integral will be evaluated,
%  fun    - function handle of (complex) function defined on real finite
%           interval [A,B].  
%  A      - finite minimum value of the function support. If the true
%           minimum is -Inf, than A should be set as a reasonable finite
%           approximation of the minimum value of the support. Default
%           value is A = -1.
%  B      - finite maximum value of the function support. If the true
%           maximum is Inf, than B should be set as a reasonable finite
%           approximation of the maximum value of the support. Default
%           value is B = 1.
%  nPts   - Order of Chebyshev polynomial approximation. Default value is
%           nPts = 100.
%
% OUTPUT:
%  FI     - (complex) vector of Fourier integral values evalated at the
%           required t, i.e. FI(t). 
%  coefs  - vector of the polynomial expansion coefficiens of the PDF.
%           This allows to check the quality of the polynomial
%           approaximation over the interval [A,B].
%
% EXAMPLE 1 (Fourier transform of Gaussian PDF)
%  t   = linspace(-5,5);
%  pdf = @(x) exp(-x.^2/2)/sqrt(2*pi);
%  [cf,coefs] = FourierIntegral_P(t,pdf,-8,8);
%  plot(coefs,'o-');grid
%  title('Chebyshev coefficients')
%  figure
%  plot(t,real(cf))
%  title('Characteristic function: Fourier transform of Gaussian PDF')
%
% EXAMPLE 2 (Inverse Fourier transform of Gaussian characteristic function)
%  x   = linspace(-5,5);
%  cf  = @(t) cfS_Gaussian(t);
%  [Q,coefs] = FourierIntegral_P(x,cf,-10,10);
%  pdf = real(Q/(2*pi));
%  plot(coefs,'o-');grid
%  title('Chebyshev coefficients')
%  figure
%  plot(x,pdf)
%  title('PDF: Inverse Fourier transform of Gaussian CF')
%
% EXAMPLE 3 (Fourier transform of Beta PDF)
%  t   = linspace(-100,100,1001);
%  pdf = @(x) betapdf(x,2.5,1.5);
%  [cf,coefs] = FourierIntegral_P(t,pdf,0,1);
%  plot(coefs,'o-');grid
%  title('Chebyshev coefficients')
%  figure
%  plot(t,real(cf),t,imag(cf))
%  title('Characteristic function: Fourier transform of Beta PDF')
%
% EXAMPLE 4 (CF of the Weibull distribution with a = 1.5, and b = 0.5)
%  a   = 1.5;
%  b   = 3.5;
%  pdf = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
%  A   = 0;
%  B   = 15;
%  t   = linspace(-10,10,2^10+1)';
%  [cf,coefs] = FourierIntegral_P(t,pdf,A,B);
%  plot(coefs,'o-');grid
%  title('Chebyshev coefficients')
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
%
% CREDENTIALS:
%  'FourierIntegral_P' was developed based on ideas of the algorithm
%  'integral_oscillatory_Chebyshev' created by Tomy Duby, OAA Computing
%  Ltd (tomy@oaacomputing.co.uk).

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Sep-2017 14:17:51

%% ALGORITHM CALL
%[FI,coefs] = FourierIntegral_P(omega,fun,A,B,nPts)

%% CHECK THE INPUT PARAMETERS
narginchk(2,5);
if nargin < 5, nPts = []; end
if nargin < 4, B = []; end
if nargin < 3, A = []; end

if isempty(A)
    A = -1;
end

if isempty(B)
    B = 1; 
end

if isempty(nPts)
    nPts = 100; 
end

%% ALGORITHM

% Transform Fun from [A,B] to [-1,1]
sz          = size(omega);
scale       = (B - A) / 2;
shift       = (B + A) / 2;
omega       = omega(:)';
id          = omega < 0;
omega(id)   = -omega(id);
K           = scale * exp(1i*omega*shift);
Omega       = scale * omega;
s           = 0:nPts;
k           = s'; 
F           = fun(scale * cos(pi*s/nPts) + shift);

% Coefficients of Chebyshev series expansion of the function F
w           = cos(pi * k * s / nPts);
F(1)        = F(1)/2;
F(nPts+1)   = F(nPts+1)/2;
coefs       = 2 * sum(bsxfun(@times,w,F),2) / nPts;

% Fourier integrals of the Chebyshev polynomials
I           = ChebFourierIntegral(nPts,Omega);
I(1,:)      = I(1,:)/2;
I(nPts+1,:) = I(nPts+1,:)/2;

% Fourier integral Q = Integral_A^B  fun(x) * exp(i*omega*x) dx
FI          = sum(bsxfun(@times,I,coefs));
FI          = K .* FI;
FI(id)      = conj(FI(id));
FI          = reshape(FI,sz);
end