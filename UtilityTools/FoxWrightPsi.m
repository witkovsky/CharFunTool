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

c  = sum(GammaLog(bsxfun(@plus,a,bsxfun(@times,A,k))),2) - ...
     sum(GammaLog(bsxfun(@plus,b,bsxfun(@times,B,k))),2);

pPSIq  = sum(exp(bsxfun(@plus,c-GammaLog(k+1),bsxfun(@times,k,log(z')))));

pPSIq(z==0) = complex(exp(c(1)));
pPSIq  = reshape(pPSIq,sz);

end