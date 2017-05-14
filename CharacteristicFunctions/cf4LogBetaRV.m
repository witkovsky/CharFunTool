function cf = cf4LogBetaRV(t,alpha,beta,coef,n)
%%cf4LogBetaRV Characteristic function of a linear combination (resp.
%  convolution) of independent log-transformed random variables (RVs)
%  log(X), where X ~ Beta(alpha,beta). Notice that log(X)<=0 if X ~
%  Beta(alpha,beta).  
%  
%  That is, cf4LogBetaRV evaluates the characteristic function cf(t) of
%  Y = coef(1)*log(X_1)+ ... + coef(N)*log(X_N), where X_i ~
%  Beta(alpha_i,beta_i), with the parameters alpha_i > 0 and beta_i > 0.
%
%  The characteristic function of Y = log(X), with X ~ Beta(alpha,beta) is
%  defined by  
%   cf_Y(t) = (gamma(beta)*gamma(alpha + 1i*t)) / ...
%                          (beta(alpha,beta)*gamma(alpha + beta + 1i*t)).
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters alpha(i) and beta(i).
%
% SYNTAX
%  cf = cf4LogBetaRV(t,alpha,beta,coef,n)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'shape' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  n     - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of n iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically
%          distributed random variable. If empty, default value is n = 1.  
%
% EXAMPLE 1:
% % CF of a linear combination of K=100 independent log-Beta RVs
%   coef = 1./(((1:50) - 0.5)*pi).^2;
%   figure; plot(idx,coef,'.-'); grid on;
%   title('Coefficients of the linear combination of log-Beta RVs')
%   alpha = 5/2;
%   beta  = 3/2;
%   t     = linspace(-100,100,201);
%   cf    = cf4LogBetaRV(t,alpha,beta,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('Characteristic function of a linear combination of log-Beta RVs')
%
% EXAMPLE 2:
% % PDF/CDF from the CF by cf2DistGP
%   alpha = 5/2;
%   beta  = 3/2;
%   coef  = 1./(((1:50) - 0.5)*pi).^2;
%   cf    = @(t) cf4LogBetaRV(t,alpha,beta,coef);
%   clear options
%   options.xMax = 0;
%   result = cf2DistGP(cf,[],[],options);
%
% EXAMPLE 3: 
% % Distribution of log(R), where R = geometric/arithmetic mean of Gamma RVs
% % Let X_1,...,X_n are independent X_j ~ Gamma(alpha,beta), and let
% % R = geometricmean(X)/mean(X). According to Glaser (JASA 1976) we get 
% % log(R) ~ (1/n) * sum_{j=1}^{n-1} log(Y_j), Y_j ~ Beta(alpha,j/n), for j
% % = 1,...,n-1. That is, log(R) is distributed as linear combination of
% % independent logBeta random variables log(Y_j).
%   n = 10;
%   alpha = 1; % i.e. Exponnetial distribution
%   beta  = (1:n-1)'/n;
%   coef  = -1/n;
%   cf = @(t) cf4LogBetaRV(t,alpha,beta,coef);
%   t = linspace(-25,25,201);
%   figure; plot(t,real(cf(t)),t,imag(cf(t))); grid on;
%   prob = [ 0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Beta_distribution
%
%  REFERENCES:
%  GLASER, R.E. (1976). The ratio of the geometric mean to the arithmetic
%  mean for a random sample from a gamma distribution. Journal of the
%  American Statistical Association, 71(354), 480-487.

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 14-May-2017 12:08:24

%% ALGORITHM
% cf = cf4LogBetaRV(t,alpha,beta,coef,n)

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

%% Equal size of the parameters
[errorcode,coef,alpha,beta] = distchck(3,coef,alpha,beta);
if errorcode > 0
    error(message('InputSizeMismatch'));
end
alpha = alpha(:);
beta  = beta(:);
coef  = coef(:);

%% Characteristic function of linear combination of noncentral chi-squares
szt = size(t);
t   = t(:);
aux = 1i*bsxfun(@times,t,coef');
aux = gammalog(bsxfun(@plus,aux,alpha')) - ...
      gammalog(bsxfun(@plus,aux,(alpha+beta)'));
aux = bsxfun(@plus,aux,ones(length(t),1) * ...
      (gammalog(alpha+beta)-gammalog(alpha))');
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