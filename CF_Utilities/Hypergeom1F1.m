function [f, method, loops] = Hypergeom1F1(a, b, z, n)
%Hypergeom1F1 Computes the confluent hypergeometric function 1F1(a;b;z)
%  also known as the Kummer's function M(a,b,z), for the parameters a and b
%  (here assumed to be real scalars), and the complex argument z (scalar,
%  vector or array).
%
%  For more details on confluent hypergeometric function 1F1(a;b;z) or the
%  Kummer's (confluent hypergeometric) function M(a, b, z) see WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Confluent_hypergeometric_function.
%
%  The present algorithm was adapted from the MATLAB version of the FORTRAN
%  program suggested by S.Zhang and J.Jin (1996). For more details see also
%  http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html.
%
%  Computation of 1F1(a;b;z) for large arguments z (with |Im(z)| >=
%  |Re(z)|) and for b > a > 0) is based on using the steepest descent
%  integration method as suggested by G. Navas Palencia and A.A. Arratia
%  Quesada (2016).
%
% SYNTAX
%   f = Hypergeom1F1(a,b,z)
%   [f, method, loops] = Hypergeom1F1(a, b, z, n)
%
% INPUTS
%  z      - parameter a (real scalar),
%  b      - parameter b (real scalar),
%  z      - complex argument (scalar, vector or array),
%  n      - number of the Gauss-Laguerre nodes and weights on (0,Inf). If
%           empty, default value is n = 64.
%
% OUTPUTS
%  f      - calculated 1F1(a,b,z) of the same dimension as z,
%  method - indicator of the used method:
%           -1 is the method used for special cases,
%            0 is the method used for z == 0,
%            1 is the method based on the series expansion,
%            2 is the method based on the steepest descent integration,
%            3 is the method based on the asymptotic expansion,
% loops   - indicator of the used number of recursive loops.
%
% EXAMPLE 1
%  a = 10;
%  b = 15;
%  z = [0; 1; 10i; 50i; 10+50i; 100+50i];
%  n = 64;
%  [f,method,loops] = Hypergeom1F1(a,b,z,n)
%
% EXAMPLE 2 (CF of Beta(1/2,1/2) distribution)
%  a  = 1/2;
%  b  = 1/2;
%  t  = linspace(-100,100,1001)';
%  cf = Hypergeom1F1(a,a+b,1i*t);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of Beta(1/2,1/2) distribution')
%  xlabel('t')
%  ylabel('CF')
%
% REFERENCES
% [1] Jin, J. M., and Zhang Shan Jjie. Computation of Special Functions.
%     Wiley, 1996.
% [2] Navas Palencia, G. and Arratia Quesada, A.A., 2016. On the
%     computation of confluent hypergeometric functions for large imaginary
%     part of parameters b and z. In Mathematical Software - ICMS 2016: 5th
%     International Conference, Berlin, Germany, July 11-14, 2016:
%     proceedings, pp. 241-248. Springer.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 30-Mar-2018 13:45:00

%% FUNCTION
%  [f, method] = Hypergeom1F1(a, b, z, n)

%% CHECK THE INPUT PARAMETERS

if nargin < 4
    n = [];
end

if isempty(n)
    n = 64;
end

%%
done   = false;
transf = false;
szz    = size(z);
z      = z(:);
sz     = size(z);
f      = NaN(sz);
method = -ones(sz);
loops  = zeros(sz);

% 1F1(a,b,z) for special cases of the parameters a and b and any argument z
if(b == 0 || b == -fix(abs(b)))
    f = Inf;
    done = true;
elseif (a == b)
    f = exp(z);
    done = true;
elseif(a-b == 1)
    f = (1 + z ./ b) .* exp(z);
    done = true;
elseif(a == 0)
    f = 1;
    done = true;
elseif(a == -1)
    f = 1 - z./b;
    done = true;
elseif(a == 1 && b == 2)
    f = (exp(z) - 1) ./ z;
    done = true;
elseif(a == fix(a) && a < 0)
    m   = -a;
    cr  = 1;
    f = 1;
    for k  = 1:m
        cr = cr .* (a + k - 1) ./ k ./ (b + k - 1) .* z;
        f  = f + cr;
    end
    done = true;
end

% 1F1(a,b,z) for other cases of the parameters a ,b and the argument z
if ~done
    % If b < a  set 1F1(a,b,z) = exp(z)*1F1(b-a,b,-z)
    if (b < a)
        transf = true;
        a      = b - a;
        z      = -z;
    end
    im_z = imag(z);
    re_z = real(z);
    ind0 = z==0;
    ind1 = ((abs(z)<10 & abs(im_z)<abs(re_z)) | abs(z)<20+abs(b) | a<0)...
           & ~ind0; 
    ind2 = (abs(z)>=10) & (abs(im_z)>=abs(re_z)) & (a>0) & (b>a);
    ind3 = (~ind0 & ~ind1 & ~ind2);
    % If z == 0 set 1F1(a,b,0) = 1
    if any(ind0)
        f(ind0) = 1;
        method(ind0) = 0;        
    end
    % 1F1(a,b,z) for small abs(z) or negative a:
    % abs(z) < 10 & abs(im_z) < abs(re_z) & ~ind0, OR
    % abs(z) < 20 + abs(b) & ~ind0, OR
    % a < 0
    if any(ind1)
        chg = 1;
        crg = 1;
        chw = 0;
        zz  = z(ind1);
        for  j = 1:500
            crg = crg .* (a+j-1) ./ (j * (b + j - 1)) .* zz;
            chg = chg + crg;
            if all(abs((chg - chw) ./ chg) < 1e-15)
                break
            end
            chw = chg;
        end
        method(ind1) = 1;
        loops(ind1)  = j;
        f(ind1)      = chg;
    end
    % 1F1(a,b,z) for large abs(z) such that abs(imag(z)) >= abs(re_z)
    % and a > b > 0 by using the steepest descent integration: 
    % abs(z) >= 10 & abs(im_z) >= abs(re_z) & a > 0 & b > a    
    if any(ind2)
        [x,w] = GaussLaguerre(n);
        gba   = gammaln(b) - (gammaln(a) + gammaln(b - a));
        rez   = re_z(ind2);
        imz   = im_z(ind2);
        ewa   = 1 ./ imz;
        ewb   = exp(imz * 1i) ./ imz;
        a1    = a - 1;
        ba1   = b - a - 1;
        r1    = 0;
        r2    = 0;
        for j = 1:n
            x_i  = x(j) ./ imz * 1i;
            aux1 = rez .* x_i + a1 * log(x_i) + ba1 * log(1 - x_i);
            aux2 = rez .* (1 + x_i) + a1 * log(1+x_i) + ba1 * log(-x_i);
            r1   = r1 + w(j) * exp(gba + aux1);
            r2   = r2 + w(j) * exp(gba + aux2);
        end
        method(ind2) = 2;
        loops(ind2)  = j;
        f(ind2)      = (ewa .* r1 - ewb .* r2) * 1i;
    end
    % 1F1(a,b,z) for (otherwise) large z by using asymptotic expansion
    if any(ind3)
        zz  = z(ind3);
        g1  = gamma(a);
        g2  = gamma(b);
        ba  = b - a;
        g3  = gamma(ba);
        cs1 = 1;
        cs2 = 1;
        cr1 = 1;
        cr2 = 1;
        for  j  = 1:500
            cr1 = -cr1 .* (a + j - 1) .* (a - b + j) ./ (zz * j);
            cr2 =  cr2 .* (b - a + j - 1) .* (j - a) ./ (zz * j);
            cs1 =  cs1 + cr1;
            cs2 =  cs2 + cr2;
            if all(abs(cr1+cr2) < 1e-15)
                break
            end
        end
        x    = real(zz);
        y    = imag(zz);
        phi  = atan(y ./ x);
        phi(x == 0 & y > 0)  = 0.5 * pi;
        phi(x == 0 & y <= 0) = -0.5 * pi;
        ns   = ones(size(x));
        ns(phi > -1.5*pi & phi <= -0.5*pi) = -1;
        cfac = exp(1i * pi * a * ns);
        cfac(y == 0) = cos(pi * a);
        chg1 = (g2 / g3) * zz.^(-a) .* cfac .* cs1;
        chg2 = (g2 / g1) .* exp(zz) .* zz.^(a-b) .* cs2;
        chg  = chg1 + chg2;
        method(ind3) = 3;
        loops(ind3)  = j;
        f(ind3)      = chg;
    end
end

if transf
    f = exp(-z) .* f;
end
f      = reshape(f,szz);
method = reshape(method,szz);
end

%% function GaussLaguerre
function [x, w] = GaussLaguerre(n,alpha)
% GaussLaguerre evaluates the Gauss-Laguerre Nodes and Weights on the
% interval (alpha,Inf).

if nargin < 2
    alpha = [];
end

if isempty(alpha)
    alpha = 0;
end

if n == 64
    x = [2.241587414670593e-02; 1.181225120967662e-01; 2.903657440180303e-01; ...
        5.392862212279714e-01; 8.650370046481124e-01; 1.267814040775241e+00; ...
        1.747859626059435e+00; 2.305463739307505e+00; 2.940965156725248e+00; ...
        3.654752650207287e+00; 4.447266343313093e+00; 5.318999254496396e+00; ...
        6.270499046923656e+00; 7.302370002587399e+00; 8.415275239483027e+00; ...
        9.609939192796107e+00; 1.088715038388638e+01; 1.224776450424431e+01; ...
        1.369270784554751e+01; 1.522298111152473e+01; 1.683966365264873e+01; ...
        1.854391817085919e+01; 2.033699594873023e+01; 2.222024266595088e+01; ...
        2.419510487593325e+01; 2.626313722711848e+01; 2.842601052750102e+01; ...
        3.068552076752596e+01; 3.304359923643782e+01; 3.550232389114120e+01; ...
        3.806393216564646e+01; 4.073083544445863e+01; 4.350563546642153e+01; ...
        4.639114297861618e+01; 4.939039902562468e+01; 5.250669934134629e+01; ...
        5.574362241327837e+01; 5.910506191901708e+01; 6.259526440015138e+01; ...
        6.621887325124754e+01; 6.998098037714681e+01; 7.388718723248294e+01; ...
        7.794367743446311e+01; 8.215730377831930e+01; 8.653569334945649e+01; ...
        9.108737561313303e+01; 9.582194001552071e+01; 1.007502319695140e+02; ...
        1.058845994687999e+02; 1.112392075244396e+02; 1.168304450513065e+02; ...
        1.226774602685386e+02; 1.288028787692377e+02; 1.352337879495258e+02; ...
        1.420031214899315e+02; 1.491516659000494e+02; 1.567310751326712e+02; ...
        1.648086026551505e+02; 1.734749468364243e+02; 1.828582046914315e+02; ...
        1.931511360370729e+02; 2.046720284850595e+02; 2.180318519353285e+02; ...
        2.348095791713262e+02];
    w = [5.625284233902887e-02; 1.190239873124205e-01; 1.574964038621475e-01; ...
        1.675470504157746e-01; 1.533528557792381e-01; 1.242210536093313e-01; ...
        9.034230098648389e-02; 5.947775576835545e-02; 3.562751890403607e-02; ...
        1.948041043116659e-02; 9.743594899382018e-03; 4.464310364166234e-03; ...
        1.875359581323119e-03; 7.226469815750032e-04; 2.554875328334960e-04; ...
        8.287143534397105e-05; 2.465686396788568e-05; 6.726713878829501e-06; ...
        1.681785369964073e-06; 3.850812981546759e-07; 8.068728040991898e-08; ...
        1.545723706757564e-08; 2.704480147613762e-09; 4.316775475431567e-10; ...
        6.277752541794292e-11; 8.306317376250609e-12; 9.984031787119531e-13; ...
        1.088353887008957e-13; 1.074017402970290e-14; 9.575737246084761e-16; ...
        7.697028063946171e-17; 5.564881054436309e-18; 3.609756216814263e-19; ...
        2.095095239662055e-20; 1.084792493734732e-21; 4.994712583627291e-23; ...
        2.037932077329677e-24; 7.340603648778086e-26; 2.324586950075985e-27; ...
        6.464528714804253e-29; 1.582906573670680e-30; 2.881154588676925e-32; ...
        5.412535994048359e-34; 2.241048506440640e-34; 5.283742844838896e-36; ...
        1.338202299148180e-34; 1.586340564468588e-35; 1.215192241351559e-34; ...
        1.217775335792122e-34; 1.673365556974291e-35; 2.735714461640009e-34; ...
        2.185380020634853e-34; 5.648495554594729e-35; 9.997398610925997e-36; ...
        1.500478177990158e-36; 1.416076744376295e-37; 5.444799396304293e-39; ...
        1.153008451969226e-40; 2.474260963687568e-42; 9.293338889710336e-45; ...
        2.974573897074668e-47; 1.941796748832940e-50; 5.776547415033449e-54; ...
        1.794991571658772e-58];
else
    idx = 1:n;
    a   = (2*idx-1) + alpha;
    b   = sqrt( idx(1:n-1) .* ((1:n-1) + alpha) );
    CM  = diag(a) + diag(b,1) + diag(b,-1);
    
    [V,L]   = eig(CM);
    [x,ind] = sort(diag(L));
    V       = V(:,ind)';
    w       = gamma(alpha+1) .* V(:,1).^2;
end
end