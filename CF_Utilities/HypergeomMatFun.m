function [s,ss] = HypergeomMatFun(MAX,alpha,p,q,x,y,lam)
% HypergeomMatFun - the hypergeometric function of a matrix argument.
% Computes truncated hypergeometric function pFq^alpha(p;q;x;y) of two
% matrix arguments, x and y.
%
% SYNTAX:
%  [s,ss] = HypergeomMatFun(MAX,alpha,p,q,x,y,lam)
%
% INPUTS:
%  MAX    - maximum number of partitions, |kappa|<=MAX.
%  alpha  - parameter of the hypergeometric function pFq^alpha(p;q;x;y),
%  p      - parameter of the hypergeometric function pFq^alpha(p;q;x;y),
%  q      - parameter of the hypergeometric function pFq^alpha(p;q;x;y),
%  x      - matrix argument (vector of eigenvalues),
%  y      - optional second matrix argument (vector of eigenvalues),
%  lam    - optional parameter, kappa<=lam
%
% OUTPUTS:
%  s     - hypergeometric sum, pFq^alpha(p;q;x;y).
%  ss    - partial sums.
%
% EXAMPLE (2F3^9([3 4];[5 6 7];[1 2];[8,9]) & kappa<=(4,3) & |kappa|<=6)
% [s,ss] = HypergeomMatFun(6,9,[3 4],[5 6 7],[1,2],[8 9],[4 3])
%
% AUTHOR/CREDITS:
% Plamen Koev, February 2004, May 22, 2004. See also MHG.
% MHG is a collection of MATLAB functions written in C for computing the
% hypergeometric function of a matrix argument.
%
% MHG is Copyright (c) 2004 Plamen Koev
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
%
% If you use MHG in any program or publication, please acknowledge
% its author by adding a reference to:
%
% Plamen Koev and Alan Edelman, The Efficient Evaluation of the
% Hypergeometric Function of a Matrix Argument, Mathematics of Computation
% 75 (2006), 833-846.

% This is modified version of the original MATLAB code hg.m by Plamen Koev
% 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Oct-2017 15:00:12

%% ALGORITHM

n      = length(x);
lambda = floor(MAX./(1:n));

if nargin == 7
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
end

l     = zeros(1,Lp);
z     = ones(1,Lp);
kt    = -(1:Lp);
cc1   = 0;
ss    = zeros(1,MAX+1);
ss(1) = 1;
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
        zn = prod(p+c) * alpha;
        dn = prod(q+c) * (kt(h)+h+1);
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
        z(h) = z(h) * zn/dn;
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
                ss(sl) = ss(sl) + z(h) * Sx(w,n) * Sy(w,n);
            else
                ss(sl) = ss(sl) + z(h) * Sx(w,n);
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
                ss(sl)= ss(sl) + z(n) * cc1 * Sx(pp,n) * Sy(pp,n);
            else
                cc1 = cc1 * prod(1+kt(1:n-1)-kt(n)) ...
                      * prodx(n) * (1/alpha+l(n)-1) ...
                      / (prod(alpha+kt(1:n-1)-kt(n))*l(n));
                ss(sl) = ss(sl) + z(n) * cc1 * Sx(pp,n);
            end
        end
        if h < Lp
            z(h+1)  = z(h);
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
s = sum(ss);
end
