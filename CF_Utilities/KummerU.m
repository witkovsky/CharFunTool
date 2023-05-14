function [y,md,id] = KummerU(a,b,z)
%KummerU Computes the confluent hypergeometric function of the second
%  kind, Kummer U(a,b,z), for complex arguments a, b, z, where Re(a)>0,
%  Re(z)>=0. Here, a, b, and z must be arrays of the same size. 
% 
%  Here, the input variables a, b, and z must be arrays of the same size.
%  Alternatively, some of the input variables can be scalars (and the other
%  scalars or arays of the same size). The output variables, y (the values
%  of the hypergeometric U function), md (indicator of the algorithm used
%  for evaluation), and id (indicator of precision - the number of valid
%  figures in the result) are arrays of the size equal to the input
%  variables. 
% 
% KummerU (part of CharFunTool) is using MATLAB built-in arrayfun algorithm
% and it is an alternative implementation to the HypergeomU (part of
% CharFunTool). 
% 
% SEE e.g., http://dlmf.nist.gov/13.
%
% SYNTAX
%   [y,md,id] = KummerU(a,b,z)
%
% EXAMPLE 1
%  a = 3;
%  b = 2.5;
%  z = 1i*(0:0.05:1)';
%  f =  KummerU(a,b,z)
%
% EXAMPLE 2 (CF of the Fisher-Snedecor F(d1,d2) distribution)
%  % See e.g. https://en.wikipedia.org/wiki/F-distribution
%  d1 = 5;
%  d2 = 3;
%  t  = linspace(-30,30,301)';
%  cf = gamma((d1+d2)/2) * ...
%       KummerU(d1/2,1-d2/2,-(d2/d1)*1i*t) / gamma(d2/2);
%  figure; plot(t,real(cf),t,imag(cf)); grid on
%  title('CF of the Fisher-Snedecor F(5,3) distribution')
%  xlabel('t')
%  ylabel('CF')
%
% EXAMPLE 3 (CF of a linear combination of the Fisher-Snedecor RVs)
%  d1   = [5 7 4];
%  d2   = [3 4 5];
%  coef = [1 1 1]/3;
%  t    = linspace(-10,10,301)';
%  o    = ones(size(t));
%  a    = o * (d1/2);
%  b    = o * (1-d2/2);
%  z    = -1i * t * (coef .* d2 ./d1);
%  cf = prod(gamma((d1+d2)/2) ./ gamma(d2/2)) * prod(KummerU(a,b,z),2);
%  cf(t==0)  = 1;
%  figure; plot(t,real(cf),t,imag(cf)); grid on
%  title('CF of a linear combination of the Fisher-Snedecor RVs')
%  xlabel('t')
%  ylabel('CF')
%
% CREDENTIALS:
%  The function U(a,b,z) is computed by using MATLAB function CHGU.
%  CHGU.m is a direct conversion of the corresponding Fortran program
%  developed by S. Zhang and J. Jin, "Computation of Special Functions"
%  (Wiley 1996). 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Feb-2023 11:25:10

%% ALGORITHM CALL
% [y,md,id] = KummerU(a,b,z)

%% ALGORITHM

[errorcode,a,b,z] = distchck(3,a,b,z);
if errorcode > 0
    error(message('InputSizeMismatch'));
end

[y,md,id] = arrayfun(@(a,b,z)chgu_1(a,b,z),a,b,z);
end

%% FUNCTION CHGU_1
function [hu, md, id] =  chgu_1( a, b, x )
%  CHGU_L Calculate the confluent hypergeometric function of the second
%  kind,  U(a, b, x) for a single argument x.

%% ALGORITHM
id  = NaN;
aa  = a - b + 1;
il1 = a==fix(a) & a<=0;
il2 = aa==fix(aa) & aa<=0;
il3 = abs(a*(a-b+1)./x)<=2;
bl1 = abs(x)<=5 | (abs(x)<=10 & a<=2);
bl2 = (abs(x)>5 & abs(x)<=12.5) & (a>=1 & b>=a+4);
bl3 = abs(x)>12.5 & a>=5 & b>=a+5;
bn  = b==fix(b) & b~=0;
id1 = -100;

if b~=fix(b)
    [hu,id1] = chgus(a,b,x);
    md  = 1;
    if id1>=12 
        id  = id1;
        return
    end    
    hu1 = hu;
end

if il1 || il2 || il3
    [hu,id] = chgul(a,b,x);
    md = 2;
    if id>=12 
        return
    end    
    if id1>id
        md = 1;
        id = id1;
        hu = hu1;
    end
end

if a>=0
    if bn && (bl1 || bl2 || bl3)
        [hu,id] = chgubi(a,b,x);
        md = 3;
    else
        [hu,id] = chguit(a,b,x);
        md = 4;
    end
else
    if b<=a
        b00 = b;
        a   = a - b + 1;
        b   = 2 - b;
        [hu,id] = chguit(a,b,x);
        hu  = x.^(1-b00) .* hu;
        md  = 4;
    elseif bn && ~il1
        [hu,id] = chgubi(a,b,x);
        md  = 3;
    end
end

% if(id < 6)
%     warning('CharFunTool:KummerU','chgu: no accurate result obtained');
% end
return
end
%% FUNCTION CHGUS
function [hu,id]=chgus(a,b,x)
%  CHGUS Compute confluent hypergeometric function U(a,b,x) for small
%  argument x

h0   = 0;
ga   = GammaZX(a);
gb   = GammaZX(b);
xg1  = 1 + a - b;
gab  = GammaZX(xg1);
xg2  = 2 - b;
gb2  = GammaZX(xg2);
hu0  = pi / sin(pi*b);
r1   = hu0 / (gab*gb);
r2   = hu0 * x.^(1-b) / (ga*gb2);
hu   = r1 - r2;
hmax = 0;
hmin = 1 + 300;

for  j=1:150
    r1=r1.*(a+j-1.0d0)./(j.*(b+j-1.0d0)).*x;
    r2=r2.*(a-b+j)./(j.*(1.0d0-b+j)).*x;
    hu=hu+r1-r2;
    hua=abs(hu);    
    if hua>hmax
        hmax = hua;
    end    
    if hua<hmin
        hmin = hua;
    end    
    if abs(hu-h0) < abs(hu)*1e-15
        break
    end
    h0 = hu;
end

d1 = log10(hmax);
if hmin ~= 0
    d2 = log10(hmin);
end
id = fix(15-abs(d1-d2));

return

end
%% FUNCTION CHGUL
function [hu,id]=chgul(a,b,x)
%  CHGUL Compute the confluent hypergeometric function U(a,b,x)for large
%  argument x 

aa  = a - b + 1;
il1 = a==fix(a) & a<=0.0;
il2 = aa==fix(aa) & aa<=0.0;
if il1
    nm = abs(a);
end
if il2
    nm = abs(aa);
end
if il1 || il2
    hu = 1;
    r  = 1;
    for k = 1:nm
        r = -r .* (a+k-1) .* (a-b+k) ./ (k*x);
        hu = hu + r;
    end    
    hu = x .^(-a) .* hu;
    id = 10;    
else
    hu = 1;
    r  = 1;
    for k  = 1:25
        r  = -r .* (a+k-1) .* (a-b+k) ./ (k*x);
        ra = abs(r);
        if(k>5 && ra>=r0 || ra<1e-15)
            break
        end
        r0 = ra;
        hu = hu + r;
    end
    id = fix(abs(log10(ra)));
    hu = x.^(-a) .* hu;
end

return

end
%% FUNCTION CHGUBI
function [hu,id]=chgubi(a,b,x)
%  CHGUBI Compute confluent hypergeometric function U(a,b,x) with integer
%  b (b = -1,-2,...) 

el  = 0.5772156649015329;
n   = abs(b-1);
rn1 = 1;
rn  = 1;

for j  = 1:n
    rn = rn*j;
    if j==n-1
        rn1 = rn;
    end
end

ps = psi(a);
ga = gamma(a);
if b>0
    a0  = a;
    a1  = a-n;
    a2  = a1;
    ga1 = gamma(a1);
    ua  = (-1)^(n-1) / (rn*ga1);
    ub  = rn1 ./ ga.*x.^(-n);
else
    a0  = a+n;
    a1  = a0;
    a2  = a;
    ga1 = gamma(a1);
    ua  = (-1)^(n-1) ./ (rn*ga).*x.^n;
    ub  = rn1/ga1;
end

hm1  = 1;
h0   = hm1; 
r    = 1;
hmax = 0;
hmin = 1 + 300;

for k   = 1:150
    r   = r.*(a0+k-1.0d0).*x./((n+k).*k);
    hm1 = hm1+r;
    hu1 = abs(hm1);
    if hu1>hmax
        hmax = hu1;
    end
    if hu1<hmin
        hmin = hu1;
    end    
    if abs(hm1-h0) < abs(hm1)*1e-15
        break
    end    
    h0 = hm1;
end
da1 = log10(hmax);

da2 = 0;
if hmin~=0
    da2 = log10(hmin);
end
id  = fix(15-abs(da1-da2));
hm1 = hm1 .* log(x);
s0  = 0;

for m =1:n
    if b>=0
        s0 = s0 - 1/m;
    end
    if b<0
        s0 = s0 + (1-a) / (m*(a+m-1));
    end
end
hm2  = ps + 2*el + s0;
r    = 1;
hmax = 0;
hmin = 1 + 300;

for k  = 1:150
    s1 = 0;
    if b>0
        for m  = 1:k
            s1 = s1 - (m+2*a-2) / (m*(m+a-1));
        end
        m  = 1:n;
        s2 = sum(1./(k+m));        
    else
        for m  = 1:k+n
            s1 = s1 + (1-a) / (m*(m+a-1));
        end
        m  = 1:k;
        s2 = sum(1./m);
    end
    hw  = 2*el + ps + s1 - s2;
    r   = r .* (a0+k-1) .* x/((n+k)*k);
    hm2 = hm2 + r*hw;
    hu2 = abs(hm2);
    if hu2>hmax
        hmax = hu2;
    end    
    if hu2<hmin
        hmin = hu2;
    end   
    if abs((hm2-h0)./hm2) < 1e-15
        break
    end    
    h0 = hm2;
end

db1 = log10(hmax);
db2 = 0;
if hmin~=0
    db2 = log10(hmin);
end
id1 = 15 - abs(db1-db2);

if id1<id
    id = id1;
end

hm3 = 1;
if n==0
    hm3 = 0;
end

r = 1;
for k   = 1:n-1
    r   = r .* (a2+k-1) ./ ((k-n)*k).*x;
    hm3 = hm3 + r;
end

sa = ua * (hm1+hm2);
sb = ub .* hm3;
hu = sa + sb;
if sa~=0
    id1 = fix(log10(abs(sa)));
end

if hu~=0.0
    id2 = fix(log10(abs(hu)));
end

if sa*sb<0
    id = id - abs(id1-id2);
end

return

end
%% FUNCTION CHGUIT
function [hu,id]=chguit(a,b,x)
%  CHGUIT Compute hypergeometric function U(a,b,x)by using
%  Gaussian-Legendre integration (n=60) 

% Gauss - Legendre abscissae and weight
t = zeros(1,30);
w = zeros(1,30);
t(:) = [.259597723012478d-01,.778093339495366d-01,...
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

w(:) = [.519078776312206d-01,.517679431749102d-01,.514884515009810d-01,...
    .510701560698557d-01,.505141845325094d-01,.498220356905502d-01,...
    .489955754557568d-01,.480370318199712d-01,.469489888489122d-01,...
    .457343797161145d-01,.443964787957872d-01,.429388928359356d-01,...
    .413655512355848d-01,.396806954523808d-01,.378888675692434d-01,...
    .359948980510845d-01,.340038927249464d-01,.319212190192963d-01,...
    .297524915007890d-01,.275035567499248d-01,.251804776215213d-01,...
    .227895169439978d-01,.203371207294572d-01,.178299010142074d-01,...
    .152746185967848d-01,.126781664768159d-01,.100475571822880d-01,...
    .738993116334531d-02,.471272992695363d-02,.202681196887362d-02];
t = [t, -t];
w = [w, w];

id  = 8;
a1  = a - 1;
b1  = b - a - 1;
c   = 12/x;
hu0 = 0;

% Loop over number of setments
for m   = 10:5:100    
    hu1 = 0;
    g   = 0.5 .* c/m;
    d   = g;
    for j   = 1:m
        t1  = d + g .* t;
        f1  = exp(-x * t1) .* t1.^a1 .* (1+t1).^b1;
        s   = sum(w.*f1);
        hu1 = hu1 + s .* g;
        d   = d + 2*g;
    end
    if abs(1- hu0/hu1) < 1e-7
        break
    end
    hu0 = hu1;
end

hu0 = 0;
for m   = 2:2:10
    hu2 = 0;
    g   = 0.5/m;
    d   = g;    
    for j   = 1:m
        t1  = d + g * t;
        t3  = c ./ (1-t1);
        f1  = exp(-x*t3) .* t3.^(a1+2) .* (1+t3).^b1 ./ c;
        s   = sum( w .* f1 );        
        hu2 = hu2 + s.*g;
        d   = d + 2*g;
    end
    if abs(1 - hu0./hu2) < 1e-8
        break
    end    
    hu0 = hu2;
end
hu = (hu1+hu2) / gamma(a);

return

end
%% GammaZX  
function g = GammaZX(z,funmode)
%  GAMMA function valid in the entire complex plane, the argument z
%  may be complex and of any size. 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 25-Sep-2019 13:56:14

if nargin < 2
    funmode = [];
end

if isempty(funmode)
    funmode = 1;
end

if funmode == 0
    g = GammaLog(z);
else
    g = exp(GammaLog(z));
end
end
%% GammaLog
function f = GammaLog(z)
% GammaLog  Natural Log of the Gamma function valid in the entire complex
%           plane. This routine uses an excellent Lanczos series
%           approximation for the complex ln(Gamma) function.
%
% SYNTAX:
%  f = GammaLog(z)
%             z may be complex and of any size.
%             Also  n! = prod(1:n) = exp(gammalog(n+1))
%
% REFERENCES: 
%  C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%  Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%  Y. Luke, "Algorithms ... functions", 1977
%  J. Spouge,  SIAM JNA 31, 1994. pp. 931
%  W. Press,  "Numerical Recipes"
%  S. Chang, "Computation of special functions", 1996
%
% AUTHOR:
%  Paul Godfrey, pgodfrey@conexant.com, 07-13-01

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 10:06:48

%% 
siz = size(z);
z   = z(:);
zz  = z;

%f = 0.*z; % reserve space in advance

p = find(real(z)<0);
if ~isempty(p)
    z(p) = -z(p);
end

g=607/128; % best results when 4<=g<=5

c = [0.99999999999999709182;
    57.156235665862923517;
    -59.597960355475491248;
    14.136097974741747174;
    -0.49191381609762019978;
    0.33994649984811888699e-4;
    0.46523628927048575665e-4;
    -0.98374475304879564677e-4;
    0.15808870322491248884e-3;
    -0.21026444172410488319e-3;
    0.21743961811521264320e-3;
    -0.16431810653676389022e-3;
    0.84418223983852743293e-4;
    -0.26190838401581408670e-4;
    0.36899182659531622704e-5];

s = 0;
for k = size(c,1):-1:2
    s = s + c(k)./(z+(k-2));
end

zg   = z+g-0.5;
s2pi = 0.9189385332046727417803297;

f = (s2pi + log(c(1)+s)) - zg + (z-0.5).*log(zg);

f(z==1 | z==2) = 0.0;

if ~isempty(p)
    lpi  = 1.14472988584940017414342735 + 1i*pi;
    f(p) = lpi-log(zz(p))-f(p)-log(sin(pi*zz(p)));
end

p = find(round(zz)==zz & imag(zz)==0 & real(zz)<=0);
if ~isempty(p)
    f(p) = Inf;
end

f = reshape(f,siz);
end