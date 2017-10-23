function g = GammaZX(z,kf)
% GammaZX  Gamma function valid in the entire complex plane.
%
% SYNTAX:
%  f = GammaZX(z,kf)
%             z may be complex and of any size.
%             Also  n! = prod(1:n) = exp(gammalog(n+1))
%     =========================================================
%     Purpose: Compute the gamma function gamma(z)or ln[gamma(z)]
%     for a complex argument
%     Input :
%     z  --- Complex argument
%     KF --- Function code
%            KF=0 for ln[gamma(z)]
%            KF=1 for gamma(z)
%     Output:
%     G --- ln[gamma(z)] or gamma(z)
%     ========================================================

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 9-Oct-2017 12:33:54

%% FUNCTION CALL
%  g = GammaZX(z,kf)

%% ALGORITHM

if nargin < 2
    kf = 1;
end

x    = real(z);
y    = imag(z);
a    = zeros(1,10);
x1   = 0.0;
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

if (y == 0.0d0 && x == fix(x) && x <= 0.0d0)
    gr = 1.0d+300;
    gi = 0.0d0;
    g = gr +1i*gi;
    return
elseif x < 0.0d0
    x1 = x;
    x  = -x;
    y  = -y;
end

x0 = x;
if x <= 7.0
    na = fix(7-x);
    x0 = x+na;
end

z1  = sqrt(x0.*x0+y.*y);
th  = atan(y./x0);
gr  = (x0-.5d0).*log(z1)-th.*y-x0+0.5d0.*log(2.0d0.*pi);
gi  = th.*(x0-0.5d0)+y.*log(z1)-y;
for k = 1:10
    t  = z1.^(1-2.*k);
    gr = gr+a(k).*t.*cos((2.0d0.*k-1.0d0).*th);
    gi = gi-a(k).*t.*sin((2.0d0.*k-1.0d0).*th);
end

if x <= 7.0
    gr1 = 0.0d0;
    gi1 = 0.0d0;
    for j = 0:na-1
        gr1 = gr1+.5d0.*log((x+j).^2+y.*y);
        gi1 = gi1+atan(y./(x+j));
    end
    gr = gr-gr1;
    gi = gi-gi1;
end

if x1 < 0.0d0
    z1  = sqrt(x.*x+y.*y);
    th1 = atan(y./x);
    sr  = -sin(pi.*x).*cosh(pi.*y);
    si  = -cos(pi.*x).*sinh(pi.*y);
    z2  = sqrt(sr.*sr+si.*si);
    th2 = atan(si./sr);
    if any(sr < 0.0d0)
        th2 = pi+th2;
    end
    gr = log(pi./(z1.*z2))-gr;
    gi = -th1-th2-gi;
end

if kf == 1
    g0 = exp(gr);
    gr = g0.*cos(gi);
    gi = g0.*sin(gi);
end
g = gr + 1i*gi;

return
end