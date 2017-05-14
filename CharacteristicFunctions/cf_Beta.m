function cf = cf_Beta(t,alpha,beta,coef,n)
%%cf_Beta Characteristic function of a linear combination (resp.
%  convolution) of independent BETA random variables, X ~ Beta(alpha,beta),
%  with the shape parameters alpha > 0 and beta >0, defined on the interval
%  (0,1) with the mean = alpha / (alpha + beta) and the variance =
%  (alpha*beta) / ((alpha+beta)^2*(alpha+beta+1)).   
%
%  The characteristic function of X ~ Beta(alpha,beta) is 
%    cf(t) = cf_Beta(t,alpha,beta) = 1F1(alpha; alpha +beta; i*t),
%  where 1F1(.;.;.) is the Confluent hypergeometric function. Hence,the
%  characteristic function of Y  = coef(1)*X_1 + ... + coef(N)*X_N 
%  is  
%    cf(t) =  cf_X_1(coef(1)*t) * ... * cf_X_N(coef(N)*t), 
%  where X_i ~ Beta(alpha(i),beta(i)) with cf_X_i(t).
%
% SYNTAX
%  cf = cf_Beta(t,alpha,beta,coef,n)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'shape' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  n     - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of n iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically
%          distributed random variable. If empty, default value is n = 1.  
%
% EXAMPLE 1:
% % CF of a Beta RV
%   alpha = 1;
%   beta  = 3;
%   t = linspace(-50,50,501);
%   cf = cf_Beta(t,alpha,beta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of a Beta RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a Beta RV
%   alpha = 1;
%   beta  = 3;
%   cf = @(t)cf_Beta(t,alpha,beta);
%   clear options
%   options.N = 2^12;
%   options.xMin = 0;
%   options.xMax = 1;
%   x = linspace(0,1,201);
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3:
% % CF of a linear combination of independent Beta RVs
%   alpha = 1;
%   beta  = 3;
%   coef  = 1./(((1:50) - 0.5)*pi).^2;
%   weights = coef/sum(coef);
%   t = linspace(-100,100,501);
%   cf = cf_Beta(t,alpha,beta,weights);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of a weighted linear combination of independent Beta RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a weighted linear combination of independent Beta RVs
%   alpha = 1;
%   beta  = 3;
%   coef  = 1./(((1:50) - 0.5)*pi).^2;
%   weights = coef/sum(coef);
%   cf = @(t)cf_Beta(t,alpha,beta,weights);
%   clear options
%   options.xMin = 0;
%   options.xMax = 1;
%   x = linspace(0,1,201);
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,x,prob,options)
%
% WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Beta_distribution

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 14-May-2017 12:08:24

%% ALGORITHM
% cf = cf_Beta(t,alpha,beta,coef,n)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, n = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end

%%
if isempty(beta) && ~isempty(alpha)
    beta = 1;
elseif isempty(beta) && ~isempty(coef)
    beta = 1;
elseif ~any(beta)
    beta = 1;
end

if isempty(alpha) && ~isempty(coef)
    alpha = 1;
elseif isempty(alpha) && ~isempty(beta)
    alpha = 1;
end

if isempty(coef) && ~isempty(beta)
    coef = 1;
elseif isempty(coef) && ~isempty(alpha)
    coef = 1;
end

if isempty(n)
    n = 1;
end

%% Equal size of the parameters   
if ~isempty(coef) && isscalar(alpha) && isscalar(beta) && isempty(n)
    coef = sort(coef);
    m    = length(coef);
    [coef,idx] = unique(coef);
    alpha = alpha * diff([idx;m+1]);
end

[errorcode,coef,alpha,beta] = distchck(3,coef(:)',alpha(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szc = length(coef);
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:szc
    cf = cf .* hypergeom1F1(alpha(i),alpha(i)+beta(i),1i*coef(i)*t);
end
cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(n)
    if isscalar(n)
        cf = cf .^ n;
    else
        error('n should be a scalar (positive integer) value');
    end
end

end
%% Function hypergeom1F1
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
% a  = 1/2;
% b  = 1/2;
% t  = linspace(-50,50,2^11)';
% cf =  hypergeom1F1(a,b+b,1i*t);
% figure; plot(t,real(cf),t,imag(cf))
% title('Characteristic function of Beta(1/2,1/2) distribution')
% xlabel('t')
% ylabel('CF')
%
% NOTE1: the functions hypergeom1F1 and kummerM are equivalent

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
    if (a < 2)
        nl = 0;
    end
    if(a >= 2)
        nl = 1;
        la = fix(a);
        a  = a - la - 1;
    end
    for  n = 0:nl
        if(a0 >= 2)
            a = a + 1;
        end
        if (abs(z)< 20 + abs(b) || a < 0)
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
            g1 = gamma(a);
            g2 = gamma(b);
            ba = b - a;
            g3 = gamma(ba);
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
if(x == fix(x))
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
    if(abs(x)> 1)
        z = abs(x);
        m = fix(z);
        r = 1;
        for  k = 1:m
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
    for k = 25:-1:1
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