function g = GammaZX(z,funmode)
% GammaZX  Gamma function valid in the entire complex plane, the argument z
%  may be complex and of any size. 
%
% SYNTAX:
%  g = GammaZX(z,funmode)
%            
%     =========================================================
%     Purpose: Compute the gamma function gamma(z)or ln[gamma(z)]
%     for a complex argument
%     Input :
%     z       --- Complex argument
%     funmode --- Function mode
%                 funmode=0 for ln[gamma(z)]
%                 funmode=1 for gamma(z)
%     Output:
%     g      --- ln[gamma(z)] or gamma(z)
%     ========================================================

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Jul-2018 16:13:31

%% FUNCTION CALL
%  g = GammaZX(z,kf)

%% ALGORITHM
sz = size(z);
z  = z(:);

if nargin < 2
    funmode = 1;
end

x    = real(z);
y    = imag(z);
a    = zeros(1,10);
x1   = 0;
pi   = 3.141592653589793d0;
a(:) = [ 8.333333333333333d-02, ...
    -2.777777777777778d-03, ...
    7.936507936507937d-04, ...
    -5.952380952380952d-04, ...
    8.417508417508418d-04, ...
    -1.917526917526918d-03, ...
    6.410256410256410d-03, ...
    -2.955065359477124d-02, ...
    1.796443723688307d-01, ...
    -1.39243221690590d+00 ];

if (and(and(y == 0,x == fix(x)),x <= 0))
    gr = 1e300;
    gi = 0;
    g  = gr + 1i*gi;
    return
elseif x < 0
    x1 = x;
    x  = -x;
    y  = -y;
end

x0 = x;
if x <= 7
    na = fix(7-x);
    x0 = x + na;
end

z1  = sqrt(x0 .* x0 + y .* y);
th  = atan(y ./ x0);
gr  = (x0 - 0.5) .* log(z1) - th .* y - x0 + 0.5*log(2*pi);
gi  = th .* (x0 - 0.5) + y .* log(z1) - y;
for k = 1:10
    t  = z1.^(1 - 2*k);
    gr = gr + a(k) .* t .* cos((2*k-1).*th);
    gi = gi - a(k) .* t .* sin((2*k-1).*th);
end

if x <= 7
    gr1 = 0;
    gi1 = 0;
    for j = 0:na-1
        gr1 = gr1 + 0.5*log((x+j).^2 + y .* y);
        gi1 = gi1 + atan(y ./ (x+j));
    end
    gr = gr - gr1;
    gi = gi - gi1;
end

if x1 < 0
    z1  = sqrt(x .* x + y .* y);
    th1 = atan(y ./ x);
    sr  = -sin(pi .* x) .* cosh(pi .* y);
    si  = -cos(pi .* x) .* sinh(pi .* y);
    z2  = sqrt(sr .* sr + si .* si);
    th2 = atan(si ./ sr);
    if any(sr < 0)
        th2 = pi + th2;
    end
    gr = log(pi ./(z1.*z2)) - gr;
    gi = -th1 - th2 - gi;
end

if funmode == 1
    g0 = exp(gr);
    gr = g0 .* cos(gi);
    gi = g0 .* sin(gi);
end
g = gr + 1i*gi;

g = reshape(g,sz);
end