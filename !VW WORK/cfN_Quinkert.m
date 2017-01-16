function cf = cfN_Quinkert(t,a,b,cfX)
%cfN_Quinkert(t,a,b) evaluates the characteristic function cf(t) of the
% Quinkert distribution, with the parameters a (a > 0) and b (b > 0), i.e. 
%   cf(t) = cfN_Quinkert(t,a,b) = 1F1(a,a+b,e^(1i*t)-1);
% where 1F1 denotes the confluent hypergeometric (Kummer's) function. For
% more details see [4], p. 564.
%
% SYNTAX
%  cf = cfN_Quinkert(t,a,b)
%  cf = cfN_Quinkert(t,a,b,cfX)
%
% EXAMPLE1 (CF of the Quinkert distribution with the parameter a=3, b=5)
%  a = 3;  
%  b = 5;
%  t = linspace(-15,15,501);
%  cf = cfN_Quinkert(t,a,b);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Quinkert distribution with a = 3, b = 5')
%
% EXAMPLE2 (CF of the compound Quinkert-Exponential distribution)
%  a = 3;  
%  b = 5;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Quinkert(t,a,b,cfX)
%  t = linspace(-15,15,501);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('CF of the compound Quinkert-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Quinkert-Exponential distribution)
%  a = 3;  
%  b = 5;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Quinkert(t,a,b,cfX);
%  x = linspace(0,1.5,101);
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
%cf = cfN_Quinkert(t,a,b,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, a = []; end
if nargin < 3, b = []; end
if nargin < 4, cfX = []; end

%%
if isempty(a), a = 1; end
if isempty(b), b = 1; end

%% Characteristic function of the Quinkert distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf = hypergeom1F1(a,a+b,expit-1);
cf = reshape(cf,szt);
cf(t==0) = 1;

end
%% HYPERGEOM1F1
function f = hypergeom1F1(a,b,z)
%HYPERGEOM1F1 Computes the confluent hypergeometric function 1F1(a,b,z), 
% also known as the Kummer function M(a,b,z), for the real parameters a 
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
%   f = hypergeom1F1(a,b,z)
%
% EXAMPLE1 (CF of Beta(1/2,1/2) distribution)
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
% NOTE1: the functions hypergeom1F1 and kummerM are equivalent
%
% NOTE2: the functions based on hypergeomU resp. kummerU are badly
%        conditioned for when b is negative integer. In such situations,
%        approximate the the function value by using the noninteger 
%        parameter b, say b = b + eps, for some small eps.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Aug-2015 10:19:29

%% CHECK THE PARAMETERS

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
%     M(a,b,z) with real parameters a, b and a
%     complex argument z
%     Input :  
%     a --- Parameter
%     b --- Parameter
%     z --- Complex argument
%     Output:  
%     CHG --- M(a,b,z)
%     Routine called: GAMMA for computing gamma function
%     ===================================================

%% ALGORITHM
chw = 0;
pi  = 3.141592653589793;
ci  = 1i;
a0  = a;

if(b == 0 || b == -fix(abs(b)))
    chg = 1e+300;
elseif(a == 0 || z == 0);
    chg = 1;
elseif(a == -1),
    chg = 1 - z./b;
elseif (a == b), 
    chg = exp(z);
elseif(a-b == 1),
    chg = (1 + z / b) * exp(z);
elseif(a == 1 && b == 2),
    chg = (exp(z) - 1) / z;
elseif(a == fix(a) && a < 0);
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
    if (a < 2)
        nl = 0;
    end
    if(a >= 2)
        nl = 1;
        la = fix(a);
        a  = a - la - 1;
    end
    for  n = 0:nl;
        if(a0 >= 2)
            a = a + 1; 
        end
        if (abs(z)< 20 + abs(b) || a < 0)
            chg = 1;
            crg = 1;
            for  j = 1:500,
                crg = crg * (a+j-1) / (j * (b + j - 1)) * z;
                chg = chg + crg;
                if(abs((chg - chw) / chg) < 1e-15)
                    break
                end
                chw = chg;
            end
        else
            g1 = gamma(a);
            g2 = gamma(b);
            ba = b - a;
            g3 = gamma(ba);
            cs1 = 1;
            cs2 = 1;
            cr1 = 1;
            cr2 = 1;
            for  i = 1:8,
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
        for  i = 1:la-1,
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
%% FUNCTION GAMMA
function ga=gamma(x)
%     ==================================================
%     Purpose: Compute gamma function gamma(x)
%     Input :  x  --- Argument of gamma(x)
%     (x is not equal to 0,-1,-2,...)
%     Output:  GA --- gamma(x)
%     ==================================================
g  = zeros(1,26);
pi = 3.141592653589793;
if(x == fix(x));
    if (x > 0)
        ga = 1;
        m1 = x - 1;
        for  k = 2:m1,
            ga = ga * k;
        end
    else
        ga = 1e+300;
    end
else
    if(abs(x)> 1)
        z = abs(x);
        m = fix(z);
        r = 1;
        for  k = 1:m,
            r = r .* (z - k);
        end
        z = z - m;
    else
        z = x;
    end
    g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0, ...
        -0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,...
        -.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,...
        -.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,...
        -.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,...
        .50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,...
        .51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
    gr=g(26);
    for k = 25:-1:1,
        gr = gr * z + g(k);
    end
    ga = 1 / (gr * z);
    if(abs(x)> 1)
        ga = ga * r;
        if(x < 0)
            ga = -pi / (x * ga * sin(pi * x)); 
        end
    end
end
return
end
