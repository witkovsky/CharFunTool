function [pFq,c] = HypergeompFqSeries(a,b,z,n)
%HypergeompFqSeries  Hypergeometric function pFq(a,b,z)
%  The hypergeometric function pFq(a,b,z) is defined for complex arguments
%  ai, bj, and z. With a = [a1, a2, …, ap] and b = [b1, b2, …, bq], the
%  hypergeometric function of order p, q is defined as
%  $pFq(a;b;z) = \sum_{n=0}^\infty
%               \frac{(a_1)_n\dots(a_p)_n}{(b_1)_n\dots(b_q)_n}
%               \frac {z^n} {n!}$
%
%  This is a very simple implementation of the function pFq(a,b,z) which is
%  computed directly from the definition - as a truncated power series. The
%  sum is truncated to n terms (by default n = 500).
%
%  An important property is the convergence criteria of the hypergeometric
%  functions depending on the values of p and q (the radius of convergence
%  of a series as a variable of z).
%
%  The absolute convergence of the series depends on the values p and q, as
%  well as on the parameter values a = [a1,...,ap], b = [b1,...,bq], and
%  also on the value of argument z. For more details on computation of
%  hypergeometric functions and the convergence regions see [1,2].
%
%  Numerical stability of the present algorithm could be strongly affected
%  at z values close to the border of the convergence region.
%
% SYNTAX
%    pFq = HypergeompFqSeries(a,b,z)
%    pFq = HypergeompFqSeries(a,b,z,n)
%
% EXAMPLE1
%  a = 3;
%  b = 1.5;
%  z = 1i*(0:0.05:1)';
%  F11 = HypergeompFqSeries(a,b,z)
%
% EXAMPLE2
%  a = [3 2.5];
%  b = 1.5;
%  z = 1i*(0:0.01:0.5)';
%  F21 = HypergeompFqSeries(a,b,z)
%
% EXAMPLE3
%  t = 1i + (-3:3)';
%  a = [3*t 2.5*t];
%  b = 1.5*t;
%  z = 0.8;
%  F21 = HypergeompFqSeries(a,b,z)
%
% REFERENCES
%  [1] Luke, Y.L. (1969). The Special Functions and Their Approximations,
%      Volume I, Academic Press.
%  [2] Pearson, J. (2009). Computation of Hypergeometric Functions.
%      Worcester College. Dissertation submitted in partial fulfillment of
%      the requirements for the degree of Master of Science in Mathematical
%      Modelling and Scientific Computing University of Oxford, 4 September
%      2009.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Aug-2018 17:23:55

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, n = []; end

if isempty(n)
    n = 500;
end

%% CHECK THE COMMON SIZE of the parameters a and b
na     = size(a,1); % a is assumed to be (na x p)
nb     = size(b,1); % b is assumed to be (nb x q)
if na ~= nb
    error(message('InputSizeMismatch'));
end

%% ALGORITHM
sz = size(z);
z  = z(:);
nz = size(z,1);
k  = (0:n)';

if nz>=1 && na==1
    c  = sum(gammalog(bsxfun(@plus,a,k)),2) - sum(gammalog(a)) - ...
        (sum(gammalog(bsxfun(@plus,b,k)),2) - sum(gammalog(b)));
    pFq  = sum(exp(bsxfun(@plus,c-gammalog(k+1),bsxfun(@times,k,log(z')))));
    pFq(z==0) = complex(1);
    pFq  = reshape(pFq,sz);
elseif nz==1 && na>=1
    pFq = zeros(na,1);
    for i = 1:na
        c  = sum(gammalog(bsxfun(@plus,a(i,:),k)),2) - ...
            sum(gammalog(a(i,:))) - ...
            (sum(gammalog(bsxfun(@plus,b(i,:),k)),2) - ...
            sum(gammalog(b(i,:))));
        pFq(i)  = sum(exp(c-gammalog(k+1)+k*log(z)));
    end
    pFq(z==0) = complex(1);
else
    error(message('InputSizeMismatch'));
end



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