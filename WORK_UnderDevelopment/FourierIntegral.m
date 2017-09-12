function Q = FourierIntegral(omega,Fun,A,B,options)
%   Q = FourierIntegral(FUN,A,B) approximates the Fourier integral of
%   function FUN from A to B using extended Clenshaw-Curtis quadrature
%   suggested by Patterson (1976) and based on the Bakhlanov and Vasileva
%   (1968) method, see also Evans and Webster (1999). The Fourier integral
%   is defined as Q =  Integral_A^B  fun(x) * exp(i*omega*x) dx. 
%
%   FUN must be a function handle, A and B can be -Inf or Inf. Function Y =
%   FUN(X) must accept a vector argument X and return a vector result Y. 
%  
% SYNTAX:
%  Q = FourierIntegral(omega,fun,A,B,options)
%
% EXAMPLES:
%  % Integrate exp(-x^2)* epx(i*omega*x) from -Inf to Inf for given omega:
%  omega = linspace(-5,5);
%  f = @(x) exp(-x.^2)
%  Q = FourierIntegral(f,-Inf,Inf)
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
% Ver.: 11-Sep-2017 13:41:45

%% ALGORITHM
%Q = FourierIntegral(omega,fun,A,B,options)

%% CHECK THE INPUT PARAMETERS
narginchk(2,5);
if nargin < 5, options = []; end
if nargin < 4, B = []; end
if nargin < 3, A = []; end

if isempty(A)
    A = -Inf;
end

if isempty(B)
    B = Inf; 
end

if ~isfield(options, 'AbsTol')
    options.AbsTol = 1e-8;
end

if ~isfield(options, 'RelTol')
    options.RelTol = 1e-8;
end

if ~isfield(options, 'nPts')
    options.nPts = 100;
end

nPts   = options.nPts;

%% ALGORITHM

sz      = size(omega);
scale   = (B - A) / 2;
shift   = (B + A) / 2;
omega   = omega(:)';
K       = scale * exp(1i*omega*shift);
Omega   = scale * omega;
s       = 0:nPts;
k       = s'; 
w       = cos(pi * k * s / nPts);
F       = Fun(scale * cos(pi*s/nPts) + shift);
F(1)    = F(1)/2;
F(nPts) = F(nPts)/2;
Fw      = sum(bsxfun(@times,w,F),2);
I       = ChebyIntegral(Omega,nPts);
I(1)    = I(1)/2;
I(nPts) = I(nPts)/2;

% Fourier integral Q = Integral_A^B  fun(x) * exp(i*omega*x) dx
Q       = sum(bsxfun(@times,I,Fw));
Q       = (2/nPts) * K .* Q;
Q       = reshape(Q,sz);

end
%% Function ChebyIntegrals
function I = ChebyIntegral(omega,nPts)
%   ChebyIntegral evaluates the required Fourier integrals of Chebyshev
%   polynomials used in the Pattersom Fourier integral method,
%     I(omega) = Integral_{-1}^1 Tk(x) * exp(1i*omega*x) dx,
%   for given vector of values omega and all k = 0,...,nPts.
%
% SYNTAX:
%  I = ChebyIntegral(omega,nPts)
%
% EXAMPLES:
%  omega = linspace(-5,5)';
%  nPts  = 10;
%  I = ChebyIntegral(omega,nPts)
%
% REFERENCES: 
%   EVANS, G.A.; WEBSTER, J.R. A comparison of some methods for the
%   evaluation of highly oscillatory integrals. Journal of Computational
%   and Applied Mathematics, 1999, 112(1): 55-69.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 11-Sep-2017 13:41:45

%% ALGORITHM

omega   = omega(:)';
nOmega  = length(omega);
I       = zeros(nPts+1,nOmega);
kE      = 0:2:nPts;
kO      = 1:2:nPts;
R       = zeros(1,nPts+1);
R(kE+1) = ((-1).^kE+1)./(1-kE.^2);
R(kO+1) = -2./((kO-2).*(kO+2));
R(2)    = 2/3;

for k  = 1:2:(nPts+1)
    J  = 0;
    Rk = R(k); 
    for j  = 1:2:k
        Rk = Rk *((k-1)^2-(j-1)^2) / ((k-1)^2-((j-1)+3)^2);       
        c  = Rk * (2*(j-1)+1) / 2;
        J  = J + 1i^(j-1) * c * besselj((j-1)+0.5,omega);
    end
    I(k,:) = sqrt(2*pi./omega) .* J;    
end

for k  = 2:2:(nPts+1)
    J  = 0;
    Rk = R(k);
    for j  = 2:2:k
        Rk = Rk*((k-1)^2 - (j-1)^2)/((k-1)^2-((j-1)+3)^2);      
        c  = Rk * (2*(j-1)+1)/2;
        J  = J + 1i^(j-1) * c * besselj((j-1)+0.5,omega);
    end
    I(k,:) = sqrt(2*pi./omega) .* J;     
end
I(:,omega==0) = 0;

end