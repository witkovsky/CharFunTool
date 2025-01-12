function cf = cf_MedianUniform(t, a, b, N)
%  cf_MedianUniform calculates the CF of the sample median from a random
%  sample of size N from the uniform distribution U(a,b), with a < b.    
%
%  The algorithm is a realization of equation (8) in Theorem 3 of the
%  PMW manuscript.
%
% Syntax:
%  cf = cf_MedianUniform(t, a, b, N)
%
% % EXAMPLE
%  a = -3;                      % Replace with your specific values
%  b = 3;                       % Replace with your specific values
%  N = 16;                      % Replace with your specific values
%  t = linspace(-5, 5, 101);    % Range of t values
%  % Calculate CF for the entire range of t
%  cf = cf_MedianUniform(t, a, b, N);
%  % Plot the CF
%  plot(t, real(cf), t, imag(cf));
%  xlabel('t');
%  ylabel('CF \phi_M(t)');
%  title(['CF for a = ', num2str(a), ', b = ', num2str(b), ', N = ', num2str(N)]);
%  grid on;
%
% % EXAMPLE
%  a = -4;                      % Replace with your specific values
%  b = 1;                       % Replace with your specific values
%  N = 16;                      % Replace with your specific values
%  t = linspace(-5, 5, 101);    % Range of t values
%  % Calculate CF for the entire range of t
%  cf = cf_MedianUniform(t, a, b, N);
%  % Plot the CF
%  plot(t, real(cf), t, imag(cf));
%  xlabel('t');
%  ylabel('CF \phi_{Median}(t)');
%  title(['CF for a = ', num2str(a), ', b = ', num2str(b), ', N = ', num2str(N)]);
%  grid on;
%
% % EXAMPLE
%  a = -4;                      % Replace with your specific values
%  b = 1;                       % Replace with your specific values
%  N = 16;                      % Replace with your specific values
%  cf   = @(t) cf_MedianUniform(t, a, b, N);
%  x    = linspace(a,b,101)';
%  prob = [0.9, 0.95, 0.975, 0.99];
%  clear options
%  options.T = 10;
%  % options.xMin = a;
%  % options.xMax = b;
%  result = cf2DistGP(cf,x,prob,options);
%  pdf = pdf_MedianUniform(x, a, b, N);
%  figure
%  plot(x,pdf,'.-',x,result.pdf,'--')
%  xlabel('x');
%  ylabel('PDF_{Median}(x)');
%  title(['Exact PDF vs. PDF from inverted CF for a = ', num2str(a), ', b = ', num2str(b), ', N = ', num2str(N)]);
%  legend('Exact','Inverted from CF','Location','ne')
%  grid on;

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver. 20-Jan-2024 14:45:47

%% Algorithm

szt = size(t);
t   = t(:);
cf0 = exp(1i*t*(a+b)/2);
t   = t*(b-a)/2;
cf  = 0;  % Initialize the output vector

if mod(N, 2) == 0
    for k = 0:(N/2 - 1)
        binom_k = nchoosek(N/2 - 1, k);
        for j = 0:(N/2 - 1)
            binom_j = nchoosek(N/2 - 1, j);
            for l = 0:(N/2 - 1 - j)
                binom_l = nchoosek(N/2 - 1 - j, l);
                cf = cf + binom_k * binom_j *  binom_l * K(k,j,l,t) * (-2)^l/(k+j+1);
            end
        end
    end
    % Multiply by precomputed constants
    cf = nchoosek(N, N/2) * (N^2/2^N) * cf;
else
    % nchoosek(n, (n+1)/2) * ((n+1)/2^n) * HypergeompFqSeries([],1+n/2,-t.^2/4)
    cf = gamma(1+N/2) * (t.^2/4).^(-N/4) .* besselj(N/2, abs(t));
end

cf = cf0 .* cf;
cf = reshape(cf,szt);

cf(t==0) = 1;
end
%% J1 function
function f = J1(a,t)
% J1fun Auxiliary function

f = (-1)^a / (a+1);
f = f * real(HypergeompFqSeries(a/2+0.5,[0.5 a/2+1.5],-t.^2/4));
% Alternatively use hypergeom function from Symbolic Toolbox
%f = f * hypergeom(a/2+0.5,[0.5 a/2+1.5],-t.^2/4);

end
%% K function
function f = K(k,j,l,t)
% Kfun Auxiliary function

f = J1(k+j+l+1,t) - (-1)^(k+j+1) * J1(l,t);

end
%% HypergeompFqSeries function
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
%% gammalog function
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