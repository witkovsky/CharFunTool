function y = expQ(x,q)
% expQ or the Tsallis q-exponential.
%  The q-exponential is a deformation of the exponential function using the
%  real parameter q. The q-deformed exponential and logarithmic functions
%  were first introduced in Tsallis statistics in 1994. However, the
%  q-deformation is the Box–Cox transformation for q = 1-lambda, proposed
%  by George Box and David Cox in 1964. Note that the q-exponential in
%  Tsallis statistics is different from a version used elsewhere.
%
%  For x < - 1/(1-q) expQ is defined as 0. The limit for q -> 1
%  gives expQ(x,q) ~ exp(x). Otherwise expQ(x,q) = (1 + (1-q)*x)^(1/(1-q)).
%
% For more details see WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Tsallis_statistics#q-exponential
% 
% SYNTAX:
%  y = expQ(x,q)
%
% EXAMPLE:
%  x = linspace(0,1,101);
%  figure
%  q  = -1;  Em1 = expQ(x,q); plot(x,Em1,'LineWidth',2); hold on; grid on
%  q  = 0;   E0  = expQ(x,q); plot(x,E0,'LineWidth',2);
%  q  = 0.5; E05 = expQ(x,q); plot(x,E05,'LineWidth',2);
%  q  = 1;   E1  = expQ(x,q); plot(x,E1,'LineWidth',2);
%  q  = 1.5; E15 = expQ(x,q); plot(x,E15,'LineWidth',2); hold off;
%  legend({'q = -1', 'q = 0', 'q = 0.5', 'q = 1', 'q = 1.5'})
% 
% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 3-Feb-2023 10:21:12


if q == 1
    y = exp(x);
else
    y = max(0,(1 + (1-q)*x)).^(1/(1-q));
end

end