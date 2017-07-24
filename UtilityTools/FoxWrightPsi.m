function [pPSIq,c] = FoxWrightPsi(a,A,b,B,z,N)
%FoxWrightPsi  The Fox-Wright Psi function pPSIq(a,A,b,B,z) for complex
%  arguments ai, Ai, for i = 1,...,p, and bj, Bj for j = 1,...,q, and z.
%  With a = [a1, a2, …, ap], A = [A1, A2, …, Ap], b = [b1, b2, …, bq], 
%  and B = [B1, B2, …, Bq], the Fox-Wright Psi function of order p, q
%  is defined as
%   pPSIq(a,A,b,B,z) = \sum_{n=0}^\infty 
%      \frac{\Gamma( a_1 + A_1 n )\cdots\Gamma( a_p + A_p n )}
%           {\Gamma( b_1 + B_1 n )\cdots\Gamma( b_q + B_q n )} 
%      \frac {z^n} {n!}.
%
%  This is a very simple implementation of the function pPSIq(a,A,b,B,z)
%  which is computed directly from the definition - as a truncated power
%  series. The sum is truncated to n terms (by default N = 500).
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
%    pPSIq = FoxWrightPsi(a,A,b,B,z,N)
% 
% For more details see WIKIPEDIA:
% https://en.wikipedia.org/wiki/Fox%E2%80%93Wright_function
%
% EXAMPLE1
%  a = 3; A = 1.5;
%  b = 1.5; B = 2.5;
%  z = 1i*(0:0.05:1)';
%  pPSIq = FoxWrightPsi(a,A,b,B,z)
%
% EXAMPLE2 (Characteristic function of the Weibull distribution)
%  lambda = 1; % scale parameter of the Weibull distribution
%  k = 5;      % shape parameter of the Weibull distribution
%  a = 1; A = 1/k;
%  b = []; B = [];
%  t = linspace(-20,20,501)';
%  z = lambda*1i*t;
%  cf = FoxWrightPsi(a,A,b,B,z);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the Weibull(1,5) distribution')
%  xlabel('t')
%  ylabel('CF')

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 10:06:48

%% FUNCTION
%  [pPSIq,c] = FoxWrightPsi(a,A,b,B,z,N)

%% CHECK THE INPUT PARAMETERS
narginchk(5, 6);
if nargin < 6, N = []; end

if isempty(N)
    N = 500;
end

%% ALGORITHM
sz = size(z);
z  = z(:);
a  = a(:)';
A  = A(:)';
b  = b(:)';
B  = B(:)';
n  = (0:N)';

c  = sum(GammaLog(bsxfun(@plus,a,bsxfun(@times,A,n))),2);

if ~isempty(b) || ~isempty(B)
    c  = c - sum(GammaLog(bsxfun(@plus,b,bsxfun(@times,B,n))),2);
end

pPSIq  = sum(exp(bsxfun(@plus,c-GammaLog(n+1),n*log(z.'))));

pPSIq(z==0) = complex(exp(c(1)));
pPSIq  = reshape(pPSIq,sz);

end