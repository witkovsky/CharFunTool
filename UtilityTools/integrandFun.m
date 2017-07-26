function [intFun,fun] = integrandFun(q,x,t,k)
% Integrand function of the integrand for integration of the Weibull CF by
% method for integration of highly oscillatory functions, see Huybrechs
% (2006). 
%
% EXAMPLE (x = 0, F(x) = integral(@(q)fun(q,x,t,k),0,Inf))
% x = 0;
% k = 5;
% t = 100;
% F0 = integral(@(q)fun(q,x,t,k),0,Inf)
%
% EXAMPLE (x = 1-eps, F(x) = integral(@(q)fun(q,x,t,k),0,Inf))
% x = 1-sqrt(eps);
% k = 5;
% t = 100;
% F1 = integral(@(q)fun(q,x,t,k),0,Inf)

g   = @(u,k) (u ./ (1-u)).^(2./k);
f   = @(u) 2*u./(1-u).^3 .* exp(-(u ./ (1-u)).^2);
h   = @(u,p,k) (g(u,k)+1i*p).^(k/2) ./ (1+(g(u,k)+1i*p).^(k/2));
hD1 = @(u,p,k) (1i*k/2).*(g(u,k)+1i*p).^(-1+k/2) ./ (1+(g(u,k)+1i*p).^(k/2)).^2;

fun = @(q,x,t,k) (exp(1i*t*g(x,k))./t) .* f(h(x,q./t,k)) .* hD1(x,q./t,k) .* exp(-q);
intFun = fun(q,x,t,k);