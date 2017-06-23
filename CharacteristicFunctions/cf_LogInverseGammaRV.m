function cf = cf_LogInverseGammaRV(t,alpha,beta,coef,niid)
%% cf_LogInverseGammaRV 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED INVERSE-GAMMA random variables (RVs) log(X),
%  where X ~ InvGamma(alpha,beta) is INVERSE-GAMMA distributed RV with the
%  shape parameter alpha > 0 and the rate parameter beta > 0.
%  
%  That is, cf_LogGammaRV evaluates the characteristic function cf(t) of Y
%  = coef_1*log(X_1) +...+ coef_N*log(X_N), where X_i ~ InvGamma(alpha_i,
%  beta_i), with the parameters alpha_i > 0 and beta_i > 0, for i =
%  1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ InvGamma(alpha,beta)
%  is defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) =
%  E(X^(1i*t)).  That is, the characteristic function can be derived from
%  expression for the r-th moment of X, E(X^r) by using (1i*t) instead of
%  r. In particular, the characteristic function of Y = log(X) is
%   cf_Y(t) = cf_Z(-t) =  beta^(1i*t).*(gamma(alpha-1i*t) / gamma(alpha)),
%  where cf_Z(t)  denotes the characteristic function of the log
%  transformed gamma random variable.
%
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters alpha(i) and beta(i).
%
% SYNTAX
%  cf = cf_LogInverseGammaRV(t,alpha,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'rate' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          IGamma random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.    
%
% EXAMPLE 1:
% % CF of a weighted linear combination of independent log-InverseGamma RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   alpha  = 5/2;
%   beta   = 3/2;
%   t      = linspace(-20,20,1001);
%   cf     = cf_LogInverseGammaRV(t,alpha,beta,weight);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a linear combination of log-InverseGamma RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a linear combination of independent  log-InverseGamma RVs
%   coef   = [1 2 3 4 5];
%   weight = coef/sum(coef);
%   alpha  = 5/2;
%   beta   = 3/2;
%   cf     = @(t) cf_LogInverseGammaRV(t,alpha,beta,weight);
%   clear options
%   options.N = 2^12;
%   options.SixSigmaRule = 8;
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Inverse-gamma_distribution

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Jun-2017 12:08:24

%% ALGORITHM
% cf = cf_LogGammaRV(t,alpha,beta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
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
aux = 1i*t*coef;
aux = gammalog(bsxfun(@plus,-aux,alpha))-ones(length(t),1)*gammalog(alpha);
aux = aux + 1i*t*log(beta);
cf  = prod(exp(aux),2);
cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
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