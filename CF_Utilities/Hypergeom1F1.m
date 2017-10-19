function f = Hypergeom1F1(a,b,z)
%HYPERGEOM1F1 Computes the confluent hypergeometric function 1F1(a,b,z),
% also known as the Kummer function M(a,b,z), for the parameters a
% and b (here assumed to be scalars), and the complex argument z
% (could be scalar, vector or array).
%
% The algorithm is based on a Fortran program in S. Zhang & J. Jin
% "Computation of Special Functions" (Wiley, 1996).
% http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
% Converted by using f2matlab: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%
% SYNTAX
%   f = Hypergeom1F1(a,b,z)
%
% EXAMPLE1 (CF of Beta(1/2,1/2) distribution)
% a  = 1/2;
% b  = 1/2;
% t  = linspace(-50,50,2^11)';
% cf =  Hypergeom1F1(a,b+b,1i*t);
% figure; plot(t,real(cf),t,imag(cf))
% title('Characteristic function of Beta(1/2,1/2) distribution')
% xlabel('t')
% ylabel('CF')
%
% NOTE1: the functions Hypergeom1F1 and KummerM are equivalent

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Oct-2017 12:16:34

%% FUNCTION
%  f = Hypergeom1F1(a,b,z)

%% CHECK THE INPUT PARAMETERS
sz = size(z);

if max(sz) > 1
    z = z(:);
    n = length(z);
    f = zeros(n,1);
    for i=1:n
        f(i) = cchg(a,b,z(i));
    end
    f = reshape(f,sz);
else
    f = cchg(a,b,z);
end
end
%% FUNCTION CCHG
function [chg]=cchg(a,b,z)
%     ===================================================
%     Purpose: Compute confluent hypergeometric function
%     M(a,b,z) with parameters a, b and complex argument z
%     Input :
%     a --- Parameter
%     b --- Parameter
%     z --- Complex argument
%     Output:
%     CHG --- M(a,b,z)
%     Routine called: GammaZX for computing Gamma function
%     ===================================================

%% ALGORITHM
chw = 0;
pi  = 3.141592653589793;
ci  = 1i;
a0  = a;

if(b == 0 || b == -fix(abs(b)))
    chg = 1e+300;
elseif(a == 0 || z == 0)
    chg = 1;
elseif(a == -1)
    chg = 1 - z./b;
elseif (a == b)
    chg = exp(z);
elseif(a-b == 1)
    chg = (1 + z / b) * exp(z);
elseif(a == 1 && b == 2)
    chg = (exp(z) - 1) / z;
elseif(a == fix(a) && a < 0)
    m   = fix(-a);
    cr  = 1;
    chg = 1;
    for  k = 1:m
        cr  = cr * (a + k - 1) / k / (b + k - 1) * z;
        chg = chg + cr;
    end
else
    x0 = real(z);
    if(x0 < 0)
        a  = b-a;
        a0 = a;
        z  = -z;
    end
    if (real(a) < 2)
        nl = 0;
    end
    if(real(a) >= 2)
        nl = 1;
        la = fix(real(a));
        a  = a - la - 1;
    end
    for  n = 0:nl
        if(real(a0) >= 2)
            a = a + 1;
        end
        if (abs(z)< 20 + abs(b) || real(a) < 0)
            chg = 1;
            crg = 1;
            for  j = 1:500
                crg = crg * (a+j-1) / (j * (b + j - 1)) * z;
                chg = chg + crg;
                if(abs((chg - chw) / chg) < 1e-15)
                    break
                end
                chw = chg;
            end
        else
            g1 = GammaZX(a);
            g2 = GammaZX(b);
            ba = b - a;
            g3 = GammaZX(ba);
            cs1 = 1;
            cs2 = 1;
            cr1 = 1;
            cr2 = 1;
            for  i = 1:8
                cr1 = -cr1 * (a + i - 1) * (a - b + i) / (z * i);
                cr2 =  cr2 * (b - a + i - 1) * (i - a) / (z * i);
                cs1 =  cs1 + cr1;
                cs2 =  cs2 + cr2;
            end
            x = real(z);
            y = imag(z);
            if (x == 0 && y >= 0)
                phi = 0.5 * pi;
            elseif (x == 0 && y <= 0)
                phi = -0.5 * pi;
            else
                phi = atan(y / x);
            end
            if (phi > -0.5 * pi && phi < 1.5 * pi)
                ns = 1;
            end
            if (phi > -1.5 * pi && phi <= -0.5 * pi)
                ns = -1;
            end
            cfac = exp(ns * ci * pi * a);
            if (y == 0)
                cfac = cos(pi * a);
            end
            chg1 = g2 / g3 * z^(-a) * cfac * cs1;
            chg2 = g2 / g1 * exp(z) * z^(a-b) * cs2;
            chg = chg1 + chg2;
        end
        if(n == 0)
            cy0 = chg;
        end
        if(n == 1)
            cy1 = chg;
        end
    end
    if(a0 >= 2)
        for  i = 1:la-1
            chg = ((2 * a - b + z) * cy1 + (b-a) * cy0) / a;
            cy0 = cy1;
            cy1 = chg;
            a   = a + 1;
        end
    end
    if(x0 < 0)
        chg = chg * exp(-z);
    end
end
return
end