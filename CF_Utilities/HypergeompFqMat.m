function [s,ss] = HypergeompFqMat(a,b,x,y,alpha,MAX,lam)
% HypergeompFqMat - The hypergeometric function of a matrix argument.
%  Computes truncated hypergeometric function pFq^alpha(a;b;x;y) with
%  the parameters a = [a1,...,ap] and b = [b1,...,bq] and of two matrix
%  arguments, x and y. 
%
%  Here, hypergeometric function pFq^alpha(a;b;x;y) is defined as
%  pFq^alpha(a;b;x;y) = pFq^alpha(a1,...,ap;b1,...,bp;x;y) = 
%    = sum_k sum_kappa [((a1)_kappa^(alpha) ... (a1)_kappa^(alpha)) ...
%       / (k!(b1)_kappa^(alpha) ... (b1)_kappa^(alpha)] C_kappa^alpha(X),
%  where kappa = (kappa1,kappa2,...) is a partition of k and
%  (ak)_kappa^(alpha), (bk)_kappa^(alpha) denote the generalized 
%  Pochhammer symbols, and C_kappa^alpha(X) is the Jack function.
%
%  For statistical applications considered here we need to use the
%  confluent hypergeometric function 1F1(a;b;X) and the Gauss
%  hypergeometric function 2F1(a,b;c;X) in matrix argument. These can be
%  expressed as special cases of the generalized hypergeometric function
%  pFq^alpha(a;b;X), with the parameter alpha = 2 (the case with zonal
%  polynomials). 
%
%  For more details and definition of the hypergeometric function
%  pFq^alpha(a;b;x;y) with matrix argument see, e.g., Koev and Edelman
%  (2006) or Muirhead (2009).
%
% SYNTAX:
%  [s,ss] = HypergeompFqMat(a,b,x,y,alpha,MAX,lam)
%
% INPUTS:
%  a      - possibly complex (1xp) vector or (nxp) matrix parameter of the
%           hypergeometric function pFq^alpha(a;b;x;y), a = [a1,...,ap],
%           p columns of n-dimensional vector parameters  
%  b      - possibly complex (1xq) vector or (nxq) matrix parameter of the
%           hypergeometric function pFq^alpha(a;b;x;y), b = [b1,...,bq], 
%           q columns of n-dimensional vector parameters   
%  x      - matrix argument (specified as vector of its eigenvalues),
%  y      - optional second matrix argument (pecified as vector of its
%           eigenvalues), 
%  alpha  - parameter of the hypergeometric function pFq^alpha(a;b;x;y),
%           default value is alpha = 2,
%  MAX    - maximum number of partitions, |kappa|<=MAX, default value is
%           MAX = 20,
%  lam    - optional parameter, kappa<=lam.
%
% OUTPUTS:
%  s     - hypergeometric sum, pFq^alpha(a;b;x;y).
%  ss    - partial sums.
%
% EXAMPLE 1:
%  % Evaluate 2F3^9([3 4];[5 6 7];[1 2];[8,9]) & kappa<=(4,3), |kappa|<=6
%  a      = [3 4];
%  b      = [5 6 7];
%  x      = [1,2];
%  y      = [8 9];
%  alpha  = 9;
%  MAX    = 6;
%  lam    = [4,3];
%  [s,ss] = HypergeompFqMat(a,b,x,y,alpha,MAX,lam)
%
% EXAMPLE 2: 
% % CF of minus log of noncentral Wilks Lambda RV distribution
% % cf_LogRV_WilksLambdaNC is using HypergeompFqMat
%   p     = 10;
%   m     = 30; % elsewhere it is denoted as n (d.f. of within SS&P)
%   n     = 5;  % elsewhere it is denoted as q (d.f. of between SS&P)
%   X     = [0 0 0.1 1 20];
%   coef  = -1;
%   MAX   = 25;
%   t     = linspace(-10,10,201);
%   cf    = cf_LogRV_WilksLambdaNC(t,p,m,n,X,coef,MAX);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of log of noncentral Wilks Lambda RV')
%
% EXAMPLE 3:
% % PDF/CDF of minus log Wilks Lambda RV (p=10, m=30, n=5) from its CF
%   p     = 10;
%   m     = 30;
%   n     = 5;
%   X     = [0 0 0.1 1 20];
%   MAX   = 30;
%   coef  = -1;
%   cf0   = @(t) cf_LogRV_WilksLambdaNC(t,p,m,n,[],coef,MAX);
%   cf    = @(t) cf_LogRV_WilksLambdaNC(t,p,m,n,X,coef,MAX);
%   prob  = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result0 = cf2DistGP(cf0,[],prob,options);
%   figure
%   result  = cf2DistGP(cf,[],prob,options);
%   disp(result)
%   figure
%   plot(result0.x,result0.cdf,result.x,result.cdf);grid on
%   xlabel('x')
%   ylabel('CDF')
%   title('CDFs of -log(\Lambda) under null and alternative hypothesis')
%   figure
%   plot(result0.x,result0.pdf,result.x,result.pdf);grid on
%   xlabel('x')
%   ylabel('PDF')
%   title('PDFs of -log(\Lambda) under null and alternative hypothesis')
%
% EXAMPLE 4: (!!!LONG CALCULATION of Hypergeom2F1Mat!!!)
% % Non-null distribution of log-transformed T statistic, W = -log(T)
% % Test statistic for testing equality of 2 covariance matrices
%   p   = 5;         % p   - length of the sample vectors (dimensionality)
%   n1  = 15;        % n1  - sample size from the first normal poulation
%   n2  = 20;        % n1  - sample size from the second normal poulation
%   nu1 = n1 - 1;    % nu1 - degrees of freedom
%   nu2 = n2 - 1;    % nu1 - degrees of freedom
%   nu  = nu1 + nu2; % nu1 - degrees of freedom
%   % delta - eigenvalues of the non-centrality matrix Delta
%   delta = [1.1, 1.325, 1.55, 1.775, 2.0]';
%   % Create the characteristic function of null distribution
%   c    = GammaMultiLog(nu/2,p)-GammaMultiLog(nu1/2,p)-GammaMultiLog(nu2/2,p);
%   cfH0 = @(t) exp(c + GammaMultiLog((1-1i*t)*nu1/2,p) + ...
%       GammaMultiLog((1-1i*t)*nu2/2,p) - GammaMultiLog((1-1i*t)*nu/2,p));
%   % Create the characteristic function of non-null distribution
%   MAX  = 50;
%   cfHA = @(t) cfH0(t(:)) .* prod(delta).^(-nu1/2) .* ...
%       Hypergeom2F1Mat(nu/2,(1-1i*t)*nu1/2,(1-1i*t)*nu/2,1-1./delta,MAX);
%   % Evaluate PDF/CDF and selected quantiles of W = -log(T)
%   x   = linspace(50,90);
%   prob = [0.9 0.95 0.975 0.99 0.999];
%   result = cf2DistGP(cfHA,x,prob);
%
% AUTHOR/CREDITS:
%  Copyright (c) 2004 Plamen Koev.
%  MHG is a collection of MATLAB functions written in C for computing the
%  hypergeometric function of a matrix argument.
%
% MHG LICENCE:
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by the
%  Free Software Foundation; either version 2 of the License, or (at your
%  option) any later version.
%
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
%  Public License for more details.
%
%  You should have received a copy of the GNU General Public License along
%  with this program; if not, write to the Free Software Foundation, Inc.,
%  59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
%
%  If you use MHG in any program or publication, please acknowledge
%  its author by adding a reference.
%
% REFERENCES:
% [1] Koev, P. and Edelman, A., 2006. The efficient evaluation of the
%     hypergeometric function of a matrix argument. Mathematics of
%     Computation, 75(254), 833-846.
% [2] Muirhead RJ. Aspects of multivariate statistical theory. John Wiley &
%     Sons; 2009 Sep 25. 

% This is modified version of the original MATLAB code hg.m by Plamen Koev
% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 25-Sep-2019 15:45:11

%% FUNCTION CALL
% [s,ss] = HypergeompFqMat(a,b,x,y,alpha,MAX,lam)

%% MEX FUNCTION CALL (all parametrs should be specified)
% [s,ss] = HypergeompFqMat(a,b,x,y,alpha,MAX,lam)
%
% % e.g. with specific parameters y = [], alpha = 2, MAX = 50, lam = []
% [s,ss] = HypergeompFqMat_mex(a,b,x,[],2,50,[])

%% CHECK THE INPUT PARAMETERS
narginchk(3, 7);
if nargin < 7, lam   = []; end
if nargin < 6, MAX   = []; end
if nargin < 5, alpha = []; end
if nargin < 4, y     = []; end

if isempty(alpha), alpha = 2; end
if isempty(MAX), MAX  = 20; end

%% CHECK THE COMMON SIZE of the parameters a and b
na     = size(a,1);
nb     = size(b,1);
if na ~= nb
    error('dimension');
end

%% ALGORITHM
n      = length(x);
lambda = floor(MAX./(1:n));

if ~isempty(lam)
    lambda = min(lam,lambda(1:length(lam)));
    MAX    = min(sum(lambda),MAX);
end

Lp = length(lambda);
while lambda(Lp) == 0
    Lp = Lp-1;
end
lambda = lambda(1:Lp);

f = 1:MAX+1;
for i = 2:Lp-1
    for j = i + 1:MAX+1
        f(j) = f(j) + f(j-i);
    end
end
f = f(end);

D  = zeros(f,1);
Sx = zeros(f,n);
Sx(1,:) = 1;

xn = ones(n,MAX+1);
if size(x,2) > 1
    x = x.';
end

for i = 2:MAX+1
    xn(:,i) = xn(:,i-1) .* x;
end
prodx = cumprod(x);

XY = false;
if nargin>5 && ~isempty(y)
    XY = true;
end

if XY
    Sy = Sx;
    yn = ones(n,MAX+1);
    if size(y,2) > 1
        y = y.';
    end
    for i = 2:MAX+1
        yn(:,i) = yn(:,i-1) .* y;
    end
    prody = cumprod(y);
else
    prody = ones(size(prodx));
    Sy = ones(size(Sx));
    yn = ones(n,MAX+1);
end

l     = zeros(1,Lp);
z     = complex(ones(na,Lp));
kt    = -(1:Lp);
cc1   = 0;
ss    = complex(zeros(na,MAX+1));
ss(:,1) = 1;
sl    = 1;
h     = 1;
ww    = ones(1,Lp);
heap  = lambda(1) + 2;
d     = zeros(1,Lp);
while h > 0
    if (l(h)<lambda(h)) && (h==1||l(h)<l(h-1)) && (MAX>=sl) && (z(h)~=0)
        l(h) = l(h) + 1;
        if l(h)==1 && h>1 && h<n
            D(ww(h)) = heap;
            ww(h)    = heap;
            m        = min(lambda(h),MAX-sl+l(h));
            heap     = heap + min(m,l(h-1));
        else
            ww(h) = ww(h) + 1;
        end
        w  = ww(h);
        c  = -(h-1)/alpha + l(h) - 1;       
        zn = prod(a+c,2) * alpha;
        dn = prod(b+c,2) * (kt(h)+h+1);
        if XY
            zn = zn * alpha * l(h);
            dn = dn * (n+alpha*c);
            for j = 1:h-1
                delta = kt(j) - kt(h);
                zn    = zn * delta;
                dn    = dn * (delta-1);
            end
        end
        kt(h) = kt(h) + alpha;
        for j = 1:h-1
            delta = kt(j) - kt(h);
            zn    = zn * delta;
            dn    = dn * (delta+1);
        end
        z(:,h) = z(:,h) .* zn./dn;
        sl   = sl + 1;
        if h < n
            if h > 1
                d(h-1) = d(h-1) - 1;
            end
            d(h) = l(h);
            cc   = prod(h+1-alpha+kt(1:h)) / prod(h+kt(1:h));
            pp   = l(1);
            k    = 2;
            while k<=h && l(k)>1
                pp = D(pp) + l(k) - 2;
                k  = k + 1;
            end
            Sx(w,h) = cc * prodx(h) * Sx(pp,h);
            if XY
                Sy(w,h) = cc * prody(h) * Sy(pp,h);
            end
            g       = find(d(1:h)>0);
            lg      = length(g);
            slm     = 1;
            nhstrip = prod(d(g)+1) - 1;
            mu      = l;
            mt      = kt;
            blm     = ones(1,lg);
            lmd     = l(g) - d(g);
            for i  = 1:nhstrip
                j  = lg;
                gz = g(lg);
                while mu(gz) == lmd(j)
                    mu(gz) = l(gz);
                    mt(gz) = kt(gz);
                    slm    = slm - d(gz);
                    j      = j - 1;
                    gz     = g(j);
                end
                t      = kt(gz)-mt(gz);
                blm(j) = blm(j) * (1+t);
                dn     = t + alpha;
                for r = 1:gz-1
                    q1     = mt(r) - mt(gz);
                    q2     = kt(r) - mt(gz);
                    blm(j) = blm(j) * (q1+alpha-1) * (1+q2);
                    dn     = dn * q1 * (alpha+q2);
                end
                blm(j) = blm(j)/dn;
                mu(gz) = mu(gz)-1;
                mt(gz) = mt(gz)-alpha;
                slm    = slm + 1;
                if j < lg
                    blm(j+1:end) = blm(j);
                end
                nmu = mu(1) + 1;
                for k = 2:h-(mu(h)==0)
                    nmu = D(nmu) + mu(k) - 1;
                end
                for k = h+1:n
                    Sx(w,k) = Sx(w,k) + blm(j) * Sx(nmu,k-1) * xn(k,slm);
                end
                if XY
                    for k = h+1:n
                        Sy(w,k) = Sy(w,k) + blm(j) * Sy(nmu,k-1) * yn(k,slm);
                    end
                end
            end
            for k = h+1:n
                Sx(w,k) = Sx(w,k) + Sx(w,k-1);
            end
            if XY
                for k = h+1:n
                    Sy(w,k) = Sy(w,k) + Sy(w,k-1);
                end
                ss(:,sl) = ss(:,sl) + z(:,h) * Sx(w,n) * Sy(w,n);
            else
                ss(:,sl) = ss(:,sl) + z(:,h) * Sx(w,n);
            end
        else
            pp = l(1) + 1 - l(n);
            k  = 2;
            while k<=n-1 && l(k)>l(n)
                pp = D(pp) + l(k) - 1 - l(n);
                k  = k + 1;
            end
            k   = (l(n)==1);
            cc1 = k + (1-k) * cc1;
            if XY
                cc1 = cc1*(prod(1+kt(1:n-1)-kt(n)) ...
                    * (1/alpha+l(n)-1) ...
                    / (prod(alpha+kt(1:n-1)-kt(n))*l(n)))^2 ...
                    * prodx(n) * prody(n);
                ss(:,sl)= ss(:,sl) + z(:,n) * cc1 * Sx(pp,n) * Sy(pp,n);
            else
                cc1 = cc1 * prod(1+kt(1:n-1)-kt(n)) ...
                    * prodx(n) * (1/alpha+l(n)-1) ...
                    / (prod(alpha+kt(1:n-1)-kt(n))*l(n));
                ss(:,sl) = ss(:,sl) + z(:,n) * cc1 * Sx(pp,n);
            end
        end
        if h < Lp
            z(:,h+1)  = z(:,h);
            ww(h+1) = w;
            h       = h + 1;
        end
    else
        sl    = sl - l(h);
        l(h)  = 0;
        kt(h) = -h;
        h     = h - 1;
    end
end
s = sum(ss,2);
end