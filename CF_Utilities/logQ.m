function y = logQ(x,q)
% logQ or the Tsallis q-logarithm.
%  The q-logarithm is the inverse of the Tsallis q-exponential and a
%  deformation of the logarithm using the real parameter q. 
%  
%  The logQ is not defined for the negative values x. The limit for q -> 1
%  gives logQ(x,q) ~ log(x). Otherwise logQ(x,q) = (x^(1-q)-1)/(1-q).
%
% For more details see WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Tsallis_statistics#q-logarithm
% 
% SYNTAX:
%  y = logQ(x,q)
%
% EXAMPLE:
%  x = linspace(0,1,101);
%  figure
%  q  = -1;  Lm1 = logQ(x,q); plot(x,Lm1,'LineWidth',2); hold on; grid on
%  q  = 0;   L0  = logQ(x,q); plot(x,L0,'LineWidth',2);
%  q  = 0.5; L05 = logQ(x,q); plot(x,L05,'LineWidth',2);
%  q  = 1;   L1  = logQ(x,q); plot(x,L1,'LineWidth',2);
%  q  = 1.5; L15 = logQ(x,q); plot(x,L15,'LineWidth',2); hold off;
%  legend({'q = -1', 'q = 0', 'q = 0.5', 'q = 1', 'q = 1.5'})
%
% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 3-Feb-2023 10:21:12

id = x > 0;
y  = NaN(size(x));

if q == 1 
    y(id) = log(x(id));
else
    y(id) = (x(id).^(1-q)-1)/(1-q);
end

end