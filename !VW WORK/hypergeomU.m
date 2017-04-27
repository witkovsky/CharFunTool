function [y,md] = hypergeomU(a,b,z)
%HYPERGEOMU Computes the Kummer confluent hypergeometric function 
%  of the second kind, U(a,b,z). 
%
% SYNTAX
%   f = hypergeomU(a,b,z)
%
% The function U(a,b,z) is computed from the confluent hyperhgeometric 
% function chgu.m.
% chgu.m program is a direct conversion of the corresponding Fortran
% program in S. Zhang & J. Jin "Computation of Special Functions" (Wiley,
% 1996). online:
% http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html. This version
% was corrected by Tomy Duby (December 20, 2015).
%
% EXAMPLE
% a = 3;
% b = 2.5;
% z = 1i*(0:0.05:1)';
% f =  hypergeomU(a,b,z)
%
% EXAMPLE2 (CF of the Fisher-Snedecor F(d1,d2) distribution)
% % See e.g. https://en.wikipedia.org/wiki/F-distribution
% d1 = 5;
% d2 = 3;
% t = linspace(-50,50,2^11)';
% cf = gamma((d1+d2)/2) *  hypergeomU(d1/2,1-d2/2,-(d2/d1)*1i*t) / gamma(d2/2);
% figure; plot(t,real(cf),t,imag(cf))
% title('CF of the Fisher-Snedecor F(5,3) distribution')
% xlabel('t')
% ylabel('CF')
%

% See also: hypergeom1F1, hypergeom2F1, hypergeomU, kummerM, kummerU, 
%        whittakerM, whittakerW, gammaincc,
%
% NOTE1: the functions hypergeomU and kummerU are equivalent

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Dec-2015 10:19:29

%% EVALUATE KUMMER U(a,b,z) function 
%% CHECK THE PARAMETERS

sz = size(z);

if max(sz) > 1
    z = z(:);
    n = length(z);
    y = zeros(n,1);
    for i=1:n
        y(i) = chgu(a,b,z(i));
    end
    y = reshape(y,sz);
else
    [y,md] = chgu(a,b,z);
end
                
end

%% FUNCTION CHGU
function [hu,md] = chgu(a,b,x)
%CHGU: Confluent Hypergeometric function of Kummer of second kind.
%
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%     =======================================================
%     Purpose: Compute the confluent hypergeometric function
%     U(a,b,x)
%     Input  : 
%     a  --- Parameter
%     b  --- Parameter
%     x  --- Argument(x > 0)
%     Output:  
%     HU --- U(a,b,x)
%     MD --- Method code
%     Routines called:
%     (1)CHGUS for small x(MD=1)
%     (2)CHGUL for large x(MD=2)
%     (3)CHGUBI for integer b(MD=3)
%     (4)CHGUIT for numerical integration(MD=4)
%
% History:
%   Modified Tomy Duby OAA Computing Ltd, 20-Dec-15: Made sure 
%       that the main function works fine when x is complex.
%
%     =======================================================

%% Check the arguments
%     if (nargin() == 3)
% %        hu = [];
% %        md = [];
%     end
    
    id  = [];
    aa  = a - b + 1;
    il1 = a  == fix(a)  & a  <= 0.0;
    il2 = aa == fix(aa) & aa <= 0.0;

    % il3 = abs(a.*(a-b+1.0))./x <= 2.0; commented out and replaced TD,
    % 20-Dec-15
    il3 = abs(a.*(a-b+1.0)./x) <= 2.0;
    
%    bl1 = x <= 5.0 |(x <= 10.0 & a <= 2.0);
%    bl2 = (x > 5.0 & x <= 12.5) & (a >= 1.0&b >= a+4.0);
%    bl3 = x > 12.5 &a >= 5.0&b >= a+5.0;
%   Replaced the above three lines with:
    bl1 = abs( x ) <= 5.0 |(abs( x ) <= 10.0 & a <= 2.0);
    bl2 = (abs( x ) > 5.0 & abs( x ) <= 12.5) & (a >= 1.0&b >= a+4.0);
    bl3 = abs( x ) > 12.5 & a >= 5.0 & b >= a+5.0;
    
    bn = b == fix(b) & b ~= 0.0;
    id1 = -100;
    if(b ~= fix(b));
        [hu,id1] = chgus(a,b,x);
        md = 1;
        if(id1 >= 6)
            return 
        end
    
        hu1 = hu;
    end
    if(il1 || il2 || il3)
        [hu,id] = chgul(a,b,x);
        md = 2;
        if(id >= 6)
            return
        end
        
        if(id1 > id);
            md = 1;
            id = id1;
            hu = hu1;
        end
    end
    
    if (a >= 0.0)
        if (bn && (bl1 || bl2 || bl3))
            [hu,id] = chgubi(a,b,x);
            md = 3;
        else
            [hu,id] = chguit(a,b,x);
            md = 4;
        end
    else
        if (b < a) % VW (b <= a)
            b00 = b;
            a = a-b+1.0d0;
            b = 2.0d0-b;
            [hu,id] = chguit(a,b,x);
            hu = x.^(1.0d0-b00).*hu;
            md = 4;
        elseif (bn && (~il1)) % (bn && (~il1))
            [hu,id] = chgubi(a,b,x);
            md = 3;
        end
    end
    
    if(id < 6)
        fprintf(1,'%s \n','no accurate result obtained');
    end
    return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CHGUS
function [hu,id]=chgus(a,b,x)
%     ======================================================
%     Purpose: Compute confluent hypergeometric function
%     U(a,b,x) for small argument x
%     Input  : 
%     a  --- Parameter
%     b  --- Parameter(b <> 0,-1,-2,...)
%     x  --- Argument
%     Output:  HU --- U(a,b,x)
%     ID --- Estimated number of significant digits
%     Routine called: GAMMA for computing gamma function
%     ======================================================

    h0=0.0;
    ga=gamma(a);
    gb=gamma(b);
    xg1=1.0d0+a-b;
    gab=gamma(xg1);
    xg2=2.0d0-b;
    gb2=gamma(xg2);
    hu0=pi./sin(pi.*b);
    r1=hu0./(gab.*gb);
    r2=hu0.*x.^(1.0d0-b)./(ga.*gb2);
    hu=r1-r2;
    hmax=0.0d0;
    hmin=1.0d+300;
    
    for  j=1:150;
        r1=r1.*(a+j-1.0d0)./(j.*(b+j-1.0d0)).*x;
        r2=r2.*(a-b+j)./(j.*(1.0d0-b+j)).*x;
        hu=hu+r1-r2;
        hua=abs(hu);

        if(hua > hmax)
            hmax=hua; 
        end

        if(hua < hmin)
            hmin=hua; 
        end

        if(abs(hu-h0)< abs(hu) * 1.0d-15)
            break 
        end
        
        h0=hu;
    end

    d1=log10(hmax);
    if(hmin ~= 0.0)
        d2=log10(hmin); 
    end
    
    id=fix(15-abs(d1-d2));
return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CHGUL
function [hu,id]=chgul(a,b,x)
%     =======================================================
%     Purpose: Compute the confluent hypergeometric function
%     U(a,b,x)for large argument x
%     by its asymptotic expansion, DLMF dquation 13.7.3

%     Input  : 
%     a  --- Parameter
%     b  --- Parameter
%     x  --- Argument
%     Output:  HU --- U(a,b,x)
%     ID --- Estimated number of significant digits
%     =======================================================

    aa=a-b+1.0d0;
    il1 = a == fix(a) & a <= 0.0;
    il2 = aa == fix(aa) & aa <= 0.0;
    if(il1)
        nm=abs(a); 
    end
    if(il2)
        nm=abs(aa); 
    end
    
    if(il1 || il2);
        % if a or a - b + 1 is a negative integer. The asymptotic expansion
        % has a finite number of terms
        hu=1.0d0;
        r=1.0d0;
        for  k=1:nm;
            r=-r.*(a+k-1.0d0).*(a-b+k)./(k.*x);
            hu=hu+r;
        end  
        
        hu = x .^(-a) .* hu;
        id=10;
        
    else
        % Neither a nor a - b + 1 is a negative integer.
        hu=1.0d0;
        r=1.0d0;

        % Consider 25 terms in the asymptotic expansion
        for  k=1:25;
            r=-r.*(a+k-1.0d0).*(a-b+k)./(k.*x);
            ra=abs(r);
            if(k > 5 && ra >= r0 || ra < 1.0d-15)
                break 
            end
            r0=ra;
            hu=hu+r;
        end

        id = fix( abs( log10( ra ) ) );
        hu = x.^(-a) .* hu;
    end
    
    return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CHGUBI
function [hu,id]=chgubi(a,b,x)
%     ======================================================
%     Purpose: Compute confluent hypergeometric function
%     U(a,b,x)with integer b(b = ס1,ס2,...)
%     Input  :
%     a  --- Parameter
%     b  --- Parameter
%     x  --- Argument
%     Output:  HU --- U(a,b,x)
%              ID --- Estimated number of significant digits
%     Routines called:
%    (1)GAMMA for computing gamma function ג(x)
%    (2)PSI for computing psi function
%     ======================================================
%ps=[];
%ga=[];
%a1=[];
%ga1=[];
%id=-100;
el  = 0.5772156649015329d0;
n   = abs(b-1);
rn1 = 1.0d0;
rn  = 1.0d0;

for j  = 1:n,
    rn = rn .* j;
    if (j == n-1)
        rn1 = rn;
    end
end

ps = psi(a);
ga = gamma(a);
if (b > 0.0)
    a0  = a;
    a1  = a - n;
    a2  = a1;
    ga1 = gamma(a1);
    ua  = (-1).^(n-1) ./ (rn .* ga1);
    ub  = rn1 ./ ga .* x.^(-n);
else
    a0  = a + n;
    a1  = a0;
    a2  = a;
    ga1 = gamma(a1);
    ua  = (-1).^(n-1) ./ (rn .* ga) .* x.^n;
    ub  = rn1 ./ ga1;
end

hm1  = 1.0d0;
h0   = hm1;       % Added here TD, 27-feb-13
r    = 1.0d0;
hmax = 0.0d0;
hmin = 1.0d+300;

for k = 1:150;
    r = r .* (a0 + k - 1) .* x ./ ((n+k) .* k);
    hm1 = hm1 + r;
    hu1 = abs(hm1);
    if(hu1 > hmax)
        hmax = hu1;
    end
    if (hu1 < hmin)
        hmin = hu1;
    end
    
    if (abs(hm1-h0) < abs(hm1) * 1.0d-15)
        break
    end
    
    h0 = hm1;
end
da1 = log10(hmax);

if (hmin ~= 0.0)
    da2 = log10(hmin);
end

id  = fix(15 - abs(da1-da2));
hm1 = hm1 .* log(x);
s0  = 0.0d0;

for m = 1:n;
    if (b >= 0.0)
        s0 = s0 - 1/m;
    end
    if (b < 0.0)
        s0 = s0 + (1-a) ./ (m .* (a+m-1));
    end
end

hm2  = ps + 2.0d0 .* el + s0;
r    = 1.0d0;
hmax = 0.0d0;
hmin = 1.0d+300;
for k = 1:150;
    s1 = 0.0d0;
    s2 = 0.0d0;
    
    if(b > 0.0);
        for m = 1:k;
            s1 = s1 - (m + 2.0d0 .* a - 2.0d0) ./ (m .* (m+a-1));
        end
        
        for m = 1:n;
            s2 = s2 + 1.0d0 ./ (k+m);
        end
    else
        for m = 1:k+n;
            s1 = s1 + (1.0d0-a) ./ (m .* (m+a-1));
        end
        
        for m = 1:k;
            s2 = s2 + 1.0d0 ./ m;
        end
    end
    
    hw  = 2.0d0 .* el + ps + s1 - s2;
    r   = r .* (a0+k-1) .* x ./ ((n+k).*k);
    hm2 = hm2 + r .* hw;
    hu2 = abs(hm2);
    if (hu2 > hmax)
        hmax = hu2;
    end
    
    if (hu2 < hmin)
        hmin = hu2;
    end
    
    if (abs((hm2-h0) ./ hm2) <  1.0d-15)
        break
    end
    
    h0 = hm2;
end
db1 = log10(hmax);

if (hmin ~= 0.0)
    db2 = log10(hmin);
end

id1 = 15 - abs(db1-db2);
if (id1 < id)
    id = id1;
end

hm3 = 1.0d0;
if (n == 0)
    hm3 = 0.0d0;
end

r = 1.0d0;

for k = 1:n-1;
    r = r.*(a2+k-1.0d0)./((k-n).*k).*x;
    hm3 = hm3+r;
end

sa = ua * (hm1+hm2);
sb = ub .* hm3;
hu = sa + sb;
if(sa ~= 0.0)
    id1 = fix(log10(abs(sa)));
end

if (hu ~= 0.0)
    id2 = fix(log10(abs(hu)));
end

if (sa * sb < 0.0)
    id = id - abs(id1-id2);
end

return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION CHGUIT
function [hu,id]=chguit(a,b,x)
%     ======================================================
%     Purpose: Compute hypergeometric function U(a,b,x) by
%     using Gaussian-Legendre integration (n=60)
%     Input  : 
%     a  --- Parameter(a > 0)
%     b  --- Parameter
%     x  --- Argument(x > 0)
%     Output:  
%     HU --- U(a,b,z)
%     ID --- Estimated number of significant digits
%     Routine called: GAMMA for computing gamma(x)
%     ======================================================

t = zeros(1,30);
w = zeros(1,30);
t(:)=[.259597723012478d-01,.778093339495366d-01,...
    .129449135396945d+00,.180739964873425d+00,...
    .231543551376029d+00,.281722937423262d+00,...
    .331142848268448d+00,.379670056576798d+00,...
    .427173741583078d+00,.473525841761707d+00,...
    .518601400058570d+00,.562278900753945d+00,...
    .604440597048510d+00,.644972828489477d+00,.683766327381356d+00,...
    .720716513355730d+00,.755723775306586d+00,.788693739932264d+00,...
    .819537526162146d+00,.848171984785930d+00,.874519922646898d+00,...
    .898510310810046d+00,.920078476177628d+00,.939166276116423d+00,...
    .955722255839996d+00,.969701788765053d+00,.981067201752598d+00,...
    .989787895222222d+00,.995840525118838d+00,.999210123227436d+00];

w(:)=[.519078776312206d-01,.517679431749102d-01,.514884515009810d-01,...
    .510701560698557d-01,.505141845325094d-01,.498220356905502d-01,...
    .489955754557568d-01,.480370318199712d-01,.469489888489122d-01,...
    .457343797161145d-01,.443964787957872d-01,.429388928359356d-01,...
    .413655512355848d-01,.396806954523808d-01,.378888675692434d-01,...
    .359948980510845d-01,.340038927249464d-01,.319212190192963d-01,...
    .297524915007890d-01,.275035567499248d-01,.251804776215213d-01,...
    .227895169439978d-01,.203371207294572d-01,.178299010142074d-01,...
    .152746185967848d-01,.126781664768159d-01,.100475571822880d-01,...
    .738993116334531d-02,.471272992695363d-02,.202681196887362d-02];
id = 7;
a1 = a - 1;
b1 = b - a -1;
c  = 12 ./ x;
for m = 10:5:100;
    hu1 = 0;
    hu0 = 0;    % Added here TD, 13-feb-11
    g   = 0.5 .* c./m;
    d   = g;
    for j = 1:m;
        s = 0;
        for k = 1:30;
            t1 = d + g .* t(k);
            t2 = d - g .* t(k);
            f1 = exp(-x .* t1) .* t1.^a1 .* (1+t1).^b1;
            f2 = exp(-x .* t2) .* t2.^a1 .* (1+t2).^b1;
            s  = s + w(k) .* (f1+f2);
        end      
        hu1 = hu1 + s .* g;
        d   = d + 2 .* g;
    end
    if(abs(1-hu0./hu1) < 1.0d-7)
        break
    end
    hu0 = hu1;
end
ga = gamma(a);
hu1= hu1 ./ ga;
for m = 2:2:10;
    hu2 = 0;
    g   = 0.5 ./ m;
    d   = g;
    
    for j = 1:m
        s = 0;
        
        for k = 1:30
            t1 = d + g .* t(k);
            t2 = d - g .* t(k);
            t3 = c ./ (1-t1);
            t4 = c ./ (1-t2);
            f1 = t3 .* t3./c .* exp(-x.*t3) .* t3.^a1 .* (1+t3).^b1;
            f2 = t4 .* t4./c .* exp(-x.*t4) .* t4.^a1 .* (1+t4).^b1;
            s  = s + w(k) .* (f1+f2);
        end
        
        hu2 = hu2 + s .* g;
        d   = d + 2 .* g;
    end
    
    if(abs(1-hu0./hu2) < 1.0d-7)
        break
    end
    
    hu0 = hu2;
end
ga  = gamma(a);
hu2 = hu2./ga;
hu  = hu1 + hu2;
return
end
%% FUNCTION GAMMA
function ga = gamma(x)
%     ==================================================
%     Purpose: Compute gamma function ג(x)
%     Input :  x  --- Argument of ג(x)
%(x is not equal to 0,-1,-2,תתת)
%     Output:  GA --- ג(x)
%     ==================================================
g=zeros(1,26);

if(x == fix(x));
    if(x > 0.0d0);
        ga=1.0d0;
        m1=x-1;
        for  k=2:m1;
            ga=ga.*k;
        end
    else
        ga=1.0d+300;
    end
else
    if(abs(x)> 1.0d0);
        z=abs(x);
        m=fix(z);
        r=1.0d0;

        for  k=1:m;
            r=r.*(z-k);
        end

        z=z-m;
    else
        z=x;
    end

    g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1, ...
        0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2, ...
        .72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3, ...
        .1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5, ...
        -.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9, ...
        .77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
    gr=g(26);

    for  k=25:-1:1;
        gr=gr.*z+g(k);
    end

    ga=1.0d0./(gr.*z);
    if(abs(x)> 1.0d0);
        ga=ga.*r;

        if(x < 0.0d0)
            ga=-pi./(x.*ga.*sin(pi.*x)); 
        end
    end
end

return
end

%% FUNCTION PSI
function ps=psi(x)
%     ======================================
%     Purpose: Compute Psi function
%     Input :  x  --- Argument of psi(x)
%     Output:  PS --- psi(x)
%     ======================================
xa=abs(x);
pi=3.141592653589793d0;
el=.5772156649015329d0;
s=0.0d0;

if(x == fix(x) && x <= 0.0);
    ps=1.0d+300;
    return
elseif(xa == fix(xa));
    n=xa;

    for  k=1 :n-1;
        s=s+1.0d0./k;
    end

    ps=-el+s;
    elseif(xa+.5 == fix(xa+.5));
    n=xa-.5;
    for  k=1:n;
        s=s+1.0./(2.0d0.*k-1.0d0);
    end  
    ps=-el+2.0d0.*s-1.386294361119891d0;

else
    if(xa < 10.0);
        n=10-fix(xa);

        for  k=0:n-1;
            s=s+1.0d0./(xa+k);
        end

        xa=xa+n;
    end
    
    x2=1.0d0./(xa.*xa);
    a1=-.8333333333333d-01;
    a2=.83333333333333333d-02;
    a3=-.39682539682539683d-02;
    a4=.41666666666666667d-02;
    a5=-.75757575757575758d-02;
    a6=.21092796092796093d-01;
    a7=-.83333333333333333d-01;
    a8=.4432598039215686d0;
    ps=log(xa)-.5d0./xa+x2.*(((((((a8.*x2+a7).*x2+a6).*x2+a5).*x2+a4).*x2+a3).*x2+a2).*x2+a1);
    ps=ps-s;
end

if(x < 0.0)
    ps=ps-pi.*cos(pi.*x)./sin(pi.*x)-1.0d0./x; 
end
return
end