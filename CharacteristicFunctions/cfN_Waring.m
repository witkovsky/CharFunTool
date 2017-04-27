function cf = cfN_Waring(t,a,b,r,cfX)
%cfN_Waring(t,a,b,r) evaluates the characteristic function cf(t) of the
% Waring distribution, with the parameters a (a > 0), b (b > 0), and r
% (r > 0), i.e.
%   cf(t) = cfN_Waring(t,a,b,r)
%         = ((gamma(a+r)*gamma(a+b)) / (gamma(a)*gamma(a+b+r)))
%           * 2F1(r,b,a+b+r,e^(1i*t));
% where 2F1 denotes the Gauss hypergeometric function. The Waring
% distribution is also known as beta negative binomial distribution. For
% more details see [4], p. 643, and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Beta_negative_binomial_distribution
%
% SYNTAX
%  cf = cfN_Waring(t,a,b,r)
%  cf = cfN_Waring(t,a,b,m,cfX)
%
% EXAMPLE1 (CF of the Waring distribution with a = 2.2, b = 3.3, r = 4)
%  % The CF is not computed correctly!! Because the hypergeomtric function
%  % hypergeom2F1(r,b,a+b+r,z) does not converege abs(z)>=1. Here z =
%  % exp(1i*t), ans abs(exp(1i*t)) = 1.
%  a = 2.2;
%  b = 3.3;
%  r = 4;
%  t = linspace(-5,5,1001);
%  cf = cfN_Waring(t,a,b,r);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('(CF of the Waring distribution with a = 2.2, b = 3.3, r = 4')
%
% EXAMPLE2 (CF of the compound Waring-Exponential distribution)
%  a = 2.2;
%  b = 3.3;
%  r = 4;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  t = linspace(-10,10,501);
%  cf = cfN_Waring(t,a,b,r,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the compound Waring-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Waring-Exponential distribution)
%  a = 2.2;
%  b = 3.3;
%  r = 4;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Waring(t,a,b,r,cfX);
%  x = linspace(0,35,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Preprint submitted
%     to Insurance: Mathematics and Economics.
% [2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models. Preprint submitted
%     to Journal of Statistical Software.
% [3] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
% [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
%     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN
%     3-87773-025-6.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfN_Waring(t,a,b,r,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(4, 5);
if nargin < 5, cfX = []; end

%% Characteristic function of the (compound) Waring distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf    = exp(gammaln(a+r) + gammaln(a+b) - gammaln(a) - gammaln(a+b+r));
cf    = cf * hypergeom2F1(r,b,a+b+r,expit);
cf    = reshape(cf,szt);
cf(t==0) = 1;

end
%% HYPERGEOM2F1
function f = hypergeom2F1(a,b,c,z)
%HYPERGEOM2F1 Computes the hypergeometric function 2F1(a,b,c,z),
% for the real parameters a, b and c (here assumed to be scalars), and the
% complex argument z (could be scalar, vector or array).
%
% The algorithm is based on a Fortran program in S. Zhang & J. Jin
% "Computation of Special Functions" (Wiley, 1996).
% http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
% Converted by using f2matlab: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%
% SYNTAX
%   f = hypergeom2F1(a,b,c,z)
%
%
% EXAMPLE1
% a = 3;
% b = 2.5;
% c = 1.5;
% z = 1i*(0:0.05:1)';
% f =  hypergeom2F1(a,b,c,z)
%
% EXAMPLE2 (CF of Beta(1/2,1/2) distribution)
% a = 1/2;
% b = 1;
% t = linspace(-50,50,2^11)';
% z = 1i*t;
% cf =  hypergeom1F1(a,b,z);
% figure; plot(t,real(cf),t,imag(cf))
% title('Characteristic function of Beta(1/2,1/2) distribution')
% xlabel('t')
% ylabel('CF')
%
% See also: hypergeom1F1, hypergeom2F1, hypergeomU, kummerM, kummerU,
%        whittakerM, whittakerW, gammaincc,
%
% NOTE:  the functions based on hypergeom2F1 is badly conditioned for when
%        c is negative integer. In such situations, approximate the the
%        function value by using the noninteger parameter c, say c = c +
%        eps, for some small eps.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Aug-2015 10:19:29

%% CHECK THE PARAMETERS

sz = size(z);

if max(sz) > 1
    z = z(:);
    n = length(z);
    f = zeros(n,1);
    for i=1:n
        f(i) = hygfz(a,b,c,z(i));
    end
    f = reshape(f,sz);
else
    f = hygfz(a,b,c,z);
end
end
%% FUNCTION HYGFZ
function [zhf]=hygfz(a,b,c,z)

%     ======================================================
%     Purpose: Compute the hypergeometric function for a
%     complex argument, F(a,b,c,z)
%     Input :
%     a --- Parameter
%     b --- Parameter
%     c --- Parameter,  c <> 0,-1,-2,...
%     z --- Complex argument
%     Output:
%     ZHF --- F(a,b,c,z)
%     Routines called:
%     (1)GAMMA for computing gamma function
%     (2)PSI for computing psi function
%     ======================================================

%% ALGORITHM
zw  = 0;
x   = real(z);
y   = imag(z);
eps = 1e-15;
l0  = (c == fix(c) && c < 0);
l1  = (abs(1-x) < eps && y == 0 && c-a-b <= 0);
l2  = (abs(z+1) < eps && abs(c-a+b-1) < eps);
l3  = (a == fix(a) && a < 0);
l4  = (b == fix(b)&& b < 0);
l5  = (c-a == fix(c-a) && c-a <= 0);
l6  = (c-b == fix(c-b) && c-b <= 0);
aa  = a;
bb  = b;
a0  = abs(z);
if (a0 > 0.95)
    eps = 1e-8;
end
pi = 3.141592653589793;
el = 0.5772156649015329;
if (l0 || l1)
    fprintf(1,'%s \n','the hypergeometric series is divergent');
    return
end
if (a0 == 0 || a == 0 || b == 0)
    zhf = 1;
elseif (z == 1 && c-a-b > 0)
    gc   = gamma(c);
    gcab = gamma(c-a-b);
    gca  = gamma(c-a);
    gcb  = gamma(c-b);
    zhf  = gc * gcab / (gca * gcb);
elseif (l2)
    g0   = sqrt(pi) * 2^(-a);
    g1   = gamma(c);
    g2   = gamma(1 + a/2 - b);
    g3   = gamma(0.5 + 0.5 * a);
    zhf  = g0 * g1 / (g2 * g3);
elseif (l3 || l4)
    if (l3)
        nm = fix(abs(a));
    end
    if (l4)
        nm = fix(abs(b));
    end
    zhf  = 1;
    zr   = 1;
    for  k = 1:nm
        zr  = zr * (a+k-1) * (b+k-1) / (k * (c+k-1)) * z;
        zhf = zhf + zr;
    end
elseif (l5 || l6)
    if (l5)
        nm = fix(abs(c-a));
    end
    if (l6)
        nm = fix(abs(c-b));
    end
    zhf = 1;
    zr  = 1;
    for  k = 1:nm
        zr  = zr * (c-a+k-1) * (c-b+k-1) / (k * (c+k-1)) * z;
        zhf = zhf + zr;
    end
    zhf  = (1-z)^(c-a-b) * zhf;
elseif (a0 <= 1)
    if (x < 0)
        z1  = z / (z-1);
        if (c > a && b < a && b > 0)
            a = bb;
            b = aa;
        end
        zc0 = 1 / ((1-z)^a);
        zhf = 1;
        zr0 = 1;
        for  k = 1:500
            zr0 = zr0 * (a+k-1) * (c-b+k-1) / (k * (c+k-1)) * z1;
            zhf = zhf + zr0;
            if (abs(zhf-zw) < abs(zhf) *eps)
                break
            end
            zw = zhf;
        end
        zhf = zc0 * zhf;
    elseif (a0 >= 0.9)
        gm   = 0;
        mcab = fix(c-a-b+eps * sign(c-a-b));
        if (abs(c-a-b-mcab) < eps)
            m = fix(c-a-b);
            ga  = gamma(a);
            gb  = gamma(b);
            gc  = gamma(c);
            gam = gamma(a + m);
            gbm = gamma(b + m);
            pa  = psi(a);
            pb  = psi(b);
            if (m ~= 0)
                gm = 1;
            end
            for  j = 1:abs(m)-1
                gm = gm * j;
            end
            rm  = 1;
            for  j = 1:abs(m)
                rm = rm * j;
            end
            zf0  = 1;
            zr0  = 1;
            zr1  = 1;
            sp0  = 0;
            sp   = 0;
            if (m >= 0)
                zc0 = gm * gc / (gam * gbm);
                zc1 = -gc * (z-1)^m / (ga * gb * rm);
                for  k = 1:m-1
                    zr0 = zr0 *(a+k-1) * (b+k-1) / (k * (k-m)) * (1-z);
                    zf0 = zf0 + zr0;
                end
                for  k = 1:m
                    sp0 = sp0 + 1 / (a+k-1) + 1 / (b+k-1) - 1 / k;
                end
                zf1 = pa + pb + sp0 + 2 * el + log(1-z);
                for  k = 1:500
                    sp = sp + (1-a) / (k * (a+k-1)) + (1-b) / (k * (b+k-1));
                    sm = 0;
                    for  j = 1:m
                        sm = sm + (1-a) / ((j+k) * (a+j+k-1)) + 1 / (b+j+k-1);
                    end
                    zp  = pa + pb + 2 * el + sp + sm + log(1-z);
                    zr1 = zr1 * (a+m+k-1) * (b+m+k-1) / (k * (m+k)) * (1-z);
                    zf1 = zf1 + zr1 * zp;
                    if(abs(zf1-zw) < abs(zf1) * eps)
                        break
                    end
                    zw = zf1;
                end
                zhf = zf0 * zc0 + zf1 * zc1;
            elseif (m < 0)
                m   = -m;
                zc0 = gm * gc / (ga * gb * (1-z)^m);
                zc1 =  -(-1)^m * gc / (gam * gbm * rm);
                for  k = 1:m-1
                    zr0 = zr0 * (a-m+k-1) * (b-m+k-1) / (k * (k-m)) * (1-z);
                    zf0 = zf0 + zr0;
                end
                for  k = 1:m
                    sp0 = sp0 + 1 / k;
                end
                zf1 = pa + pb - sp0 + 2 * el + log(1-z);
                for  k = 1:500
                    sp = sp + (1-a) / (k * (a+k-1)) + (1-b) / (k * (b+k-1));
                    sm = 0;
                    for  j = 1:m
                        sm = sm + 1 / (j+k);
                    end
                    zp  = pa + pb + 2 * el + sp - sm + log(1-z);
                    zr1 = zr1 * (a+k-1) * (b+k-1) / (k * (m+k)) * (1-z);
                    zf1 = zf1 + zr1 * zp;
                    if (abs(zf1-zw) < abs(zf1) * eps)
                        break
                    end
                    zw = zf1;
                end
                zhf = zf0 * zc0 + zf1 * zc1;
            end
        else
            ga   = gamma(a);
            gb   = gamma(b);
            gc   = gamma(c);
            gca  = gamma(c - a);
            gcb  = gamma(c - b);
            gcab = gamma(c - a - b);
            gabc = gamma(a + b - c);
            zc0  = gc * gcab / (gca * gcb);
            zc1  = gc * gabc / (ga * gb) * (1-z)^(c-a-b);
            zhf  = 0;
            zr0  = zc0;
            zr1  = zc1;
            for  k = 1:500
                zr0 = zr0 * (a+k-1) * (b+k-1) / (k * (a+b-c+k)) * (1-z);
                zr1 = zr1 * (c-a+k-1) * (c-b+k-1) / (k * (c-a-b+k)) * (1-z);
                zhf = zhf + zr0 + zr1;
                if (abs(zhf-zw) < abs(zhf) * eps)
                    break
                end
                zw = zhf;
            end
            zhf = zhf + zc0 + zc1;
        end
    else
        z00 = 1;
        if (c-a < a && c-b < b)
            z00 = (1-z)^(c-a-b);
            a   = c - a;
            b   = c - b;
        end
        zhf = 1;
        zr  = 1;
        for  k = 1:1500
            zr  = zr * (a+k-1) * (b+k-1) / (k * (c+k-1)) * z;
            zhf = zhf + zr;
            if (abs(zhf-zw) <= abs(zhf) * eps)
                break
            end
            zw = zhf;
        end
        zhf = z00 * zhf;
    end
elseif (a0 > 1.0d0)
    mab = fix(a-b+eps * sign(a-b));
    if(abs(a-b-mab) < eps && a0 <= 1.1)
        b = b + eps;
    end
    if (abs(a-b-mab) > eps)
        ga  = gamma(a);
        gb  = gamma(b);
        gc  = gamma(c);
        gab = gamma(a - b);
        gba = gamma(b - a);
        gca = gamma(c - a);
        gcb = gamma(c - b);
        zc0 = gc * gba / (gca * gb * (-z) ^a);
        zc1 = gc * gab / (gcb * ga *(-z) ^b);
        zr0 = zc0;
        zr1 = zc1;
        zhf = 0;
        for  k = 1:500
            zr0 = zr0 * (a+k-1) * (a-c+k) / ((a-b+k) * k *z);
            zr1 = zr1 * (b+k-1) * (b-c+k) / ((b-a+k) * k *z);
            zhf = zhf + zr0 + zr1;
            if (abs((zhf-zw) / zhf)<= eps)
                break
            end
            zw = zhf;
        end
        zhf = zhf + zc0 + zc1;
    else
        if (a-b < 0)
            a = bb;
            b = aa;
        end
        ca  = c - a;
        cb  = c - b;
        nca = fix(ca + eps * sign(ca));
        ncb = fix(cb + eps * sign(cb));
        if(abs(ca-nca) < eps || abs(cb-ncb) < eps)
            c = c + eps;
        end
        ga  = gamma(a);
        gc  = gamma(c);
        gcb = gamma(c - b);
        pa  = psi(a);
        pca = psi(c - a);
        pac = psi(a - c);
        mab = fix(a - b + eps);
        zc0 = gc / (ga * (-z)^b);
        gm  = gamma(a - b);
        zf0 = gm / gcb * zc0;
        zr  = zc0;
        for  k = 1:mab-1
            zr   = zr * (b+k-1) / (k * z);
            t0   = a - b - k;
            g0   = gamma(t0);
            gcbk = gamma(c - b - k);
            zf0  = zf0 + zr * g0 / gcbk;
        end
        if (mab == 0)
            zf0 = 0;
        end
        zc1 = gc / (ga * gcb * (-z) ^a);
        sp  = -2 * el - pa - pca;
        for  j = 1:mab
            sp = sp + 1 / j;
        end
        zp0 = sp + log(-z);
        sq  = 1;
        for  j = 1:mab
            sq = sq * (b + j - 1) * (b - c + j) / j;
        end
        zf1 = (sq * zp0) * zc1;
        zr  = zc1;
        rk1 = 1;
        sj1 = 0;
        %           VW added next line: w0 = 0;
        w0  = 0;
        for  k = 1:10000
            zr  = zr / z;
            rk1 = rk1 * (b+k-1) * (b-c+k) / (k*k);
            rk2 = rk1;
            for  j = k+1:k+mab
                rk2 = rk2 * (b + j - 1) * (b - c + j) / j;
            end
            sj1 = sj1 + (a-1) / (k * (a+k-1)) + (a-c-1) / (k * (a-c+k-1));
            sj2 = sj1;
            for  j = k+1:k+mab
                sj2 = sj2 + 1 / j;
            end
            zp  = -2 * el - pa - pac + sj2 - ...
                1 / (k+a-c) - pi / tan(pi * (k+a-c)) + log(-z);
            zf1 = zf1 + rk2 * zr * zp;
            ws  = abs(zf1);
            if (abs((ws-w0) / ws) < eps)
                break
            end
            w0 = ws;
        end
        zhf = zf0 + zf1;
    end
end
% if (k > 150)
%     fprintf(1,[repmat(' ',1,1),'warning% you should check the accuracy' ' \n']);
% end
return
end
%% FUNCTION GAMMA
function ga = gamma(x)
%     ==================================================
%     Purpose: Compute gamma function gamma(x)
%     Input :
%     x  --- Argument of gamma(x)
%     (x is not equal to 0,-1,-2,...)
%     Output:
%     GA --- â(x)
%     ==================================================
%

g  = zeros(1,26);
pi = 3.141592653589793;
if (x == fix(x))
    if (x > 0)
        ga = 1;
        m1 = x - 1;
        for  k = 2:m1
            ga = ga * k;
        end
    else
        ga = 1e+300;
    end
else
    if (abs(x)> 1)
        z = abs(x);
        m = fix(z);
        r = 1;
        for  k = 1:m
            r = r * (z-k);
        end
        z = z - m;
    else
        z = x;
    end
    g(:) = [1, 0.5772156649015329, -0.6558780715202538, ...
        -0.420026350340952d-1, 0.1665386113822915, ...
        -0.421977345555443d-1, -0.96219715278770d-2, ...
        0.72189432466630d-2, -0.11651675918591d-2, ...
        -0.2152416741149d-3, 0.1280502823882d-3, -0.201348547807d-4, ...
        -0.12504934821d-5, 0.11330272320d-5, -0.2056338417d-6, ...
        0.61160950d-8, 0.50020075d-8, -0.11812746d-8, 0.1043427d-9, ...
        0.77823d-11, -0.36968d-11, 0.51d-12, -0.206d-13, -0.54d-14, ...
        0.14d-14, 0.1d-15];
    gr=g(26);
    for  k = 25:-1:1
        gr = gr * z + g(k);
    end
    ga = 1 / (gr * z);
    if (abs(x) > 1)
        ga = ga * r;
        if(x < 0)
            ga = -pi / (x * ga * sin(pi*x));
        end
    end
end
return
end
%% FUNCTION PSI
function [ps]=psi(x)
%     ======================================
%     Purpose: Compute Psi function
%     Input :
%     x  --- Argument of psi(x)
%     Output:
%     PS --- psi(x)
%     ======================================
%

xa = abs(x);
pi = 3.141592653589793;
el = 0.5772156649015329;
s  = 0;
if (x == fix(x) && x <= 0)
    ps = 1e+300;
    return
elseif (xa == fix(xa))
    n = xa;
    for  k = 1:n-1
        s = s + 1 / k;
    end
    ps = -el + s;
elseif (xa + 0.5 == fix(xa + 0.5))
    n = xa - 0.5;
    for  k = 1:n
        s = s + 1 / (2 * k - 1);
    end
    ps = -el + 2 * s - 1.386294361119891;
else
    if (xa < 10)
        n = 10 - fix(xa);
        for  k = 0:n-1
            s = s + 1 / (xa + k);
        end
        xa = xa+n;
    end
    x2 = 1 / (xa * xa);
    a1 = -0.8333333333333d-01;
    a2 = 0.83333333333333333d-02;
    a3 = -0.39682539682539683d-02;
    a4 = 0.41666666666666667d-02;
    a5 = -0.75757575757575758d-02;
    a6 = 0.21092796092796093d-01;
    a7 = -0.83333333333333333d-01;
    a8 = 0.4432598039215686d0;
    ps = log(xa) - 0.5 / xa + x2 * ...
        (((((((a8 * x2 + a7) * x2 +a6) * x2 + a5) * x2 + a4) * x2 + a3) ...
        * x2 + a2) * x2 + a1);
    ps = ps - s;
end
if (x < 0)
    ps = ps - pi * cos(pi * x) / sin(pi * x) -1 / x;
end
return
end