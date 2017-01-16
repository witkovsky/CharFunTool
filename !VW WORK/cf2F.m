
%% Exact Fisher-Snedecor F-distribution by numerical inversion of the CF
% Bromwich integration

clear
%% Laplace transform
% d1 = 1;
% d2 = 5;
% M1 = @(s) (gamma((d1+d2)/2) / gamma(d2/2)) * KummerU(d1/2,1-d2/2,(d2/d1)*s);
% 
% d1 = 10;
% d2 = 7;
% M2 = @(s) (gamma((d1+d2)/2) / gamma(d2/2)) * KummerU(d1/2,1-d2/2,(d2/d1)*s);
% 
% M = @(s) M1(s) .* M2(s);

d1 = 1;
d2 = 5;
M = @(s) (gamma((d1+d2)/2) / gamma(d2/2)) * hypergeomU(d1/2,1-d2/2,(d2/d1)*s);
%M = @(s) (gamma((d1+d2)/2) / gamma(d2/2)) * KummerU(d1/2,1-d2/2,(d2/d1)*s);

%% Integrand 1 on (A,inf)
prob = 0.995;
y = finv(prob,d1,d2);
a0 = 0.1;
A = a0*y;
B = inf;

const1 = exp(a0*y)/(pi*y);

%fun1 = @(p) -const1 * (1i * (M(a0-p/y)) ./(a0-p/y)) .* exp(-p);
fun1 = @(p) const1 * imag((M(a0-p/y)) ./(a0-p/y)) .* exp(-p);

p = linspace(A,30);
plot(p,real(fun1(p)))

%% Integral 1 on (A,inf)

cdf1 =   integral(fun1,A,B)


%% Integrand on (0,1)
% transform p from (A,inf) to u from (0,1) p = A + (u./(1-u)).^2;

prob = 0.99;
y = finv(prob,d1,d2);

%fun = @(u) const1 * imag((M(a0-(A + (u./(1-u)).^2)/y)) ...
%    ./(a0-(A + (u./(1-u)).^2)/y)) .* exp(-(A + (u./(1-u)).^2)) .* (2*u./(1-u).^3);

const1 = 2/(pi);
fun = @(u) const1 * real(1i*M(-(u./(1-u)).^2/y)) .* exp(-(u./(1-u)).^2) .* (1./(u.*(1-u)));

u = linspace(0,1);
plot(u,real(fun(u)))
%% Integral

cdf = integral(fun,0,1)