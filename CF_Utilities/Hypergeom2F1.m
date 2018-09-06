function f = Hypergeom2F1(a,b,c,z)
%HYPERGEOM2F1 Computes the Gauss hypergeometric function 2F1(a,b,c,z),
% for the real parameters a, b and c (here assumed to be scalars), and the
% complex argument z (could be scalar, vector or array).
%
% SYNTAX
%   f = Hypergeom2F1(a,b,c,z)
%
% EXAMPLE 1
%  a = 3;
%  b = 2.5;
%  c = 1.5;
%  z = 1i*(0:0.05:1)';
%  f =  Hypergeom2F1(a,b,c,z)
%
% EXAMPLE 2
%  t = 1i*linspace(-5,5,11)';
%  a = 3*t;
%  b = 2.5*t;
%  c = 1.5*t;
%  z = 0.75;
%  f =  Hypergeom2F1(a,b,c,z)
%
% NOTE:  
%  The functions based on Hypergeom2F1 is badly conditioned for when c is
%  negative integer. In such situations, approximate the function value by
%  using the noninteger parameter c, say c = c + eps, for some small eps.
%
% CREDENTIALS:
%  The algorithm is based on a Fortran program in S. Zhang & J. Jin
%  "Computation of Special Functions" (Wiley, 1996). Converted by Ben
%  Barrowes (barrowes@alum.mit.edu) 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Aug-2018 18:32:27

%% FUNCTION CALL
%  f = Hypergeom2F1(a,b,c,z)

%% CHECK THE INPUT PARAMETERS
narginchk(4, 4);

[errorcode,a,b,c] = distchck(3,a(:),b(:),c(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end
na = length(a);

sz = size(z);
z  = z(:);
nz = length(z);

if nz>=1 && na==1
    f = zeros(nz,1);
    for i=1:nz
        f(i) = hygfz(a,b,c,z(i));
    end
    f = reshape(f,sz);
elseif nz==1 && na>1
    f = zeros(na,1);
    for i=1:na
        f(i) = hygfz(a(i),b(i),c(i),z);
    end
else
    error(message('InputSizeMismatch'));
end
end
%% FUNCTION HYGFZ
function [zhf]=hygfz(a,b,c,z)
%  HYGFZ Compute the hypergeometric function for a complex argument,
%  F(a,b,c,z) 
%  Input :
%     a --- Parameter
%     b --- Parameter
%     c --- Parameter,  c <> 0,-1,-2,...
%     z --- Complex argument
%     Output:
%     ZHF --- F(a,b,c,z)

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

return
end