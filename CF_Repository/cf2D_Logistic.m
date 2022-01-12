function cf = cf2D_Logistic(t,mu,sigma,coef,niid)
%% cf2D_Logistic  
%  Characteristic function of a linear combination (resp. convolution) of
%  independent BIVARIATE LOGISTIC random variables, X = (X1,X2), with
%  location parameters mu = (mu1,mu2) (real), and scale parameters sigma =
%  (sigma1,sigma2) > 0.  
%
%  That is, cf2D_Logistic evaluates the characteristic function cf(t1,t2)
%  of Y = (Y1,Y2) =  sum_{i=1}^N (coef1_i * X1_i,coef2_i *X2_i), where X_i
%  = (X1_i,X2_i) ~ BivariateLogistic ((mu1_i,nu2_i),(sigma1_i,sigma2_i))
%  are inedependent bivariate random vectors, with location parameters mu_i
%  = (mu1_i,nu2_i) and the scale parameters sigma_i = (sigma1_i,sigma2_i) >
%  0, for i = 1,...,N.
% 
%  In fact, the characteristic function of Y, cf_Y(t1,t2), is specified as
%  cf_Y(t1,t2) = cf_X1(coef11*t1,coef21*t2)*...*cf_XN(coef1N*t1,coef2N*t2).
%  
%
% SYNTAX:
%  cf = cf2D_Logistic(t,mu,sigma,coef,niid)
%
% INPUTS:
%  t      - (n x 2)-matrix t = [t1 t2] or a cell t = {t1 t2} of real
%           values, where the CF is evaluated. If t is cell of two vectors,
%           then cf is a (n1 x n2)-matrix where n1 = length(t1), n2 =
%           length(t2). 
%  mu     - matrix mu = [mu1 mu2] of real location parameters. If empty,
%           default value is mu = [0 0].   
%  sigma  - matrix sigma = [sigma1 sigma2] of the scale parameters sigma >
%           0. If empty, default value is sigma = [1 1].  
%  coef  - matrix of the coefficients, coef = [coef1 coef2], of the linear
%          combination of the BIVARIATE LOGISTIC random vectors. If coef is
%          (1x2)-vector, it is assumed that all coefficients are equal. If
%          empty, default value is coef = [1 1].
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Logistic_distribution.
%
% EXAMPLE 1:
% % CF of the Bivariate Logistic RV
%   mu   = [0 0];
%   sigma = [1 1];
%   cf   = @(t) cf2D_Logistic(t,mu,sigma);
%   tt   = linspace(-10,10,51);
%   [t1 t2] = meshgrid(tt,tt);
%   t = [t1(:) t2(:)];
%   figure; 
%   contour(t1,t2,real(reshape(cf(t),51,51)),'ShowText','on')
%   hold on
%   contour(t1,t2,imag(reshape(cf(t),51,51)),'--','ShowText','on')
%   hold off
%   title('Contour Plot of the CF of the Bivariate Logistic Distribution')
%
% EXAMPLE 2:
% % PDF/CDF of Bivariate Logistic Distribution by numerical inversion of CF
%   cf = @(t) cf2D_Logistic(t);
%   clear options;
%   options.isInterp = true;
%   result = cf2Dist2D(cf,[],options)
%
% REFERENCES
% [1] Popovic, B.V., Mijanovic, A. and Genç, A.I., 2020. On linear
%     combination of generalized logistic random variables with an
%     application to financial returns. Applied Mathematics and
%     Computation, 381, p.125314.   

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 12-Jan-2022 16:29:29

%% ALGORITHM
% cf = cf2D_Logistic(t,mu,sigma,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

%%
if isempty(mu), mu = [0,0]; end
if isempty(sigma), sigma = [1,1]; end
if isempty(coef), coef = [1,1]; end


%% Equal size of the parameters   

[err1,m1,b1,c1] = distchck(3,mu(:,1)',sigma(:,1)',coef(:,1)');
[err2,m2,b2,c2] = distchck(3,mu(:,2)',sigma(:,2)',coef(:,2)');
if err1 + err2 > 0
    error(message('InputSizeMismatch'));
end
mu   = [m1; m2];
sigma = [b1; b2];
coef = [c1; c2];
N = length(c1);

%% Characteristic function

if iscell(t)
    sz = [length(t{1}),length(t{2})];    
    [t2,t1] = meshgrid(t{2},t{1});
    t  = [t1(:),t2(:)];
else
    szt = size(t);
    sztMin = min(szt(1),szt(2));
    sztMax = max(szt(1),szt(2));    
    switch sztMin
        case 1 % if t is N-vector then it is assumed that t1 = t and t2 = t
               % and we create new t = [t1 t2]
            if sztMax > 2
                sz = size(t);
                t = t(:);
                t = [t t];
            elseif sztMax == 1
                sz = [1,1];
                t = [t t];
            else
                sz = [1,1];
                t = t(:)';
            end
        case 2 % t is 2xN or Nx2 matrix
            if sztMax > 2
                if szt(1) == 2 % if t is 2xN matrix transpose it to the Nx2 matrix
                    t = t';
                    sz = [1,szt(2)];
                else
                    sz = [szt(1),1];
                end
            end
        otherwise
            error(message('InputSizeMismatch'));
    end
end

cf = 0;
for i = 1:N
    cf = cf + cfLog(t(:,1),t(:,2),mu(1,i),mu(2,i),...
        sigma(1,i),sigma(2,i),coef(1,i),coef(2,i));
end
cf = exp(cf);
cf = reshape(cf,sz);

%cf(t(:,1)==0,t(:,2)==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end
%% Function cfLog
function cf = cfLog(t1,t2,mu1,mu2,sigma1,sigma2,coef1,coef2)

cf = 1i*mu1*coef1*t1 + 1i*mu2*coef2*t2 + ...
    GammaLogVW(1-1i*sigma1*coef1*t1) + ...
    GammaLogVW(1-1i*sigma2*coef2*t2) + ...
    GammaLogVW(1+1i*sigma1*coef1*t1+1i*sigma2*coef2*t2);
end
%% Function GammaLogVW
function f = GammaLogVW(z)
% GammaLogVW  Natural Log of the Gamma function valid in the entire complex
%           plane. This routine uses an excellent Lanczos series
%           approximation for the complex ln(Gamma) function.
%
% SYNTAX:
%  f = GammaLogVW(z)
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

%% FUNCTION
%  f = GammaLogVW(z)

%% CHECK THE INPUT PARAMETERS
siz = size(z);
z   = z(:);
zz  = z;

%f = 0.*z; % reserve space in advance

p = find(real(z)<0);
if ~isempty(p)
    z(p) = -z(p);
end

%% ALGORITHM
%Lanczos approximation for the complex plane

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