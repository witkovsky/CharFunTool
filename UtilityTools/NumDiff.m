function [diff1,diff2,diff1B,diff2B] = NumDiff(fun,x,h)
% NUMDIFF  Auxiliary function to calculate numerical (first) derivative.
%          Given below is the five point method for the first derivative
%          (five-point stencil in one dimension). See  Abramowitz & Stegun,
%          Table 25.2.
%          See also http://en.wikipedia.org/wiki/Numerical_differentiation.
%          https://en.wikipedia.org/wiki/Finite_difference_coefficient

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 10:06:48

%% FUNCTION
%  [diff1,diff2,diff1B,diff2B] = NumDiff(fun,x,h)

%% CHECK THE INPUT PARAMETERS
if nargin < 3
    h = 1e-3;
end

if nargin < 3
    x = 0;
end

%% ALGORITHM
diff1 = (-fun(x+2*h)+8*fun(x+h)-8*fun(x-h)+fun(x-2*h))/(12*h);

diff2 = (fun(x-4*h) - 16*fun(x-3*h) + 64*fun(x-2*h) + 16*fun(x-h) - 130*fun(x) + ...
    +16*fun(x+h) + 64*fun(x+2*h) - 16*fun(x+3*h) + fun(x+4*h))/(144*h^2);

% 1/280	?4/105	1/5	?4/5	0	4/5	?1/5	4/105	?1/280
diff1B = ((1/280)*fun(x-4*h)-(4/105)*fun(x-3*h)+(1/5)*fun(x-2*h)-(4/5)*fun(x-h) + ... 
        +(4/5)*fun(x+h)-(1/5)*fun(x+2*h)+(4/105)*fun(x+3*h)-(1/280)*fun(x+4*h))/h;

% second difference with order of accuracy 8
diff2B = (-(1/560)*fun(x-4*h)+(8/315)*fun(x-3*h)-(1/5)*fun(x-2*h)+(8/5)*fun(x-h) ...
    - (205/72)*fun(x)+(8/5)*fun(x+h)-(1/5)*fun(x+2*h)+(8/315)*fun(x+3*h) ...
    - (1/560)*fun(x+4*h))/h^2;
end