function cf = cf_LogGammaRV(t,alpha,beta,coef,n)
%%cf_LogGammaRV Characteristic function of a linear combination (resp.
%  convolution) of independent log-transformed random variables (RVs)
%  log(X), where X ~ Gamma(alpha,beta), and alpha and  beta represent the
%  'shape' and the 'rate' parameters of the GAMMA distribution.
%  
%  That is, cf_LogGammaRV evaluates the characteristic function cf(t) of
%  Y = coef(1)*log(X_1)+ ... + coef(N)*log(X_N), where X_i ~
%  Gamma(alpha_i,beta_i), with the parameters alpha_i > 0 and beta_i > 0.
%
%  The characteristic function of Y = log(X), with X ~ Gamma(alpha,beta is
%  defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)). 
%  That is, the characteristic function can be derived from expression for
%  the h-th moment of X, E(X^h) by using (1i*t) instead of h. In
%  particular, thecharacteristic function of Y = log(X) is
%   cf_Y(t) = gamma(alpha + 1i*t) /(beta^(1i*t)*gamma(alpha).
%
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters alpha(i) and beta(i).
%
% SYNTAX
%  cf = cf_LogGammaRV(t,alpha,beta,coef,n)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'rate' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          logGamma random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  n     - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of n iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i, with X_i ~ logGamma(alpha(i),beta(i)))
%          independently and identically distributed random variables. If
%          empty, default value is n = 1.    
%
% EXAMPLE 1:
% % CF of a weighted linear combination of independent log-Gamma RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   alpha  = 5/2;
%   beta   = 3/2;
%   t      = linspace(-20,20,1001);
%   cf     = cf_LogGammaRV(t,alpha,beta,-weight);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a linear combination of minus log-Gamma RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a linear combination of independent log-Gamma RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   alpha  = 5/2;
%   beta   = 3/2;
%   cf     = @(t) cf_LogGammaRV(t,alpha,beta,weight);
%   clear options
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,[],[],options);
%   disp(result)
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Gamma_distribution
% Wolfram MathWorld:
%  http://mathworld.wolfram.com/GammaDistribution.html
%
% REFERENCES:
% MATHAI, A.M. (1973). A review of the different techniques used for
% deriving the exact distributions of multivariate test criteria. Sankhya:
% The Indian Journal of Statistics, Series A, 39-60.  

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24

%% ALGORITHM
% cf = cf_LogGammaRV(t,alpha,beta,coef,n)

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

%% Check size of the parameters
[errorcode,coef,alpha,beta] = distchck(3,coef(:)',alpha(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of linear combination 
szt = size(t);
t   = t(:);
aux = 1i*bsxfun(@times,t,coef);
aux = gammalog(bsxfun(@plus,aux,alpha))-ones(length(t),1)*gammalog(alpha);
aux = bsxfun(@plus,aux,-(1i*t)*log(beta));
cf  = prod(exp(aux),2);
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
%% Function gammalog
function [f] = gammalog(z)
% GAMMALOG  Natural Log of the Gamma function valid in the entire complex
%           plane. This routine uses an excellent Lanczos series
%           approximation for the complex ln(Gamma) function.
%
%usage: [f] = gammalog(z)
%             z may be complex and of any size.
%             Also  n! = prod(1:n) = exp(gammalog(n+1))
%
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996

% Paul Godfrey, pgodfrey@conexant.com, 07-13-01

siz = size(z);
z   = z(:);
zz  = z;

%f = 0.*z; % reserve space in advance

p = find(real(z)<0);
if ~isempty(p)
    z(p) = -z(p);
end

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