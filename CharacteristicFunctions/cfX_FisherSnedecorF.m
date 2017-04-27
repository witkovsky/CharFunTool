function cf = cfX_FisherSnedecorF(t,df1,df2,tol) 
% cfX_FisherSnedecorF(t,df1,df2,tol)  evaluates the characteristic function
% cf(t) of the Fisher-Snedecor F-distribution with the parameters df1 and
% df2 (degrees of freedom, df1>=1, df2>=1) computed for real vector
% argument t, i.e. 
%   cf(t) = cfX_FisherSnedecorF(t,df1,df2)  
%         = (gamma(df1/2+df2/2)/gamma(df2/2)) * ...
%           U(df1/2, 1-df2/2, -1i*(df2/df1)*t),
%  where U(a,b,z) denotes the confluent hypergeometric function of the
%  second kind. 
%  For more details see WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/F-distribution 
%
% SYNTAX:
%  cf = cfX_FisherSnedecorF(df1,df2,t,tol)    
%
% EXAMPLE1 (CF of the F-distribution with df1 = 3, df2 = 5)
%  df1 = 3;
%  df2 = 5;
%  t   = linspace(-30,30,2^10+1)';
%  cf  = cfX_FisherSnedecorF(t,df1,df2);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the F-distribution')
%
% EXAMPLE2 (PDF/CDF of the F-distribution with df1 = 3, df2 = 5)
%  df1 = 3;
%  df2 = 5;
%  x = linspace(0,25,101);
%  prob = [0.9 0.95 0.99];
%  cf = @(t) cfX_FisherSnedecorF(t,df1,df2);
%  clear options
%  options.xMin = 0;
%  options.xMax = 500;
%  options.N  = 2^15;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-Fisher-Snedecor distribution)
%  n = 25;  
%  p = 0.3;
%  df1 = 3;
%  df2 = 5;
%  cfX = @(t) cfX_FisherSnedecorF(t,df1,df2);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,80,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [2] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  
% [3] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Working Paper.
%     Insurance: Mathematics and Economics. 
% [4] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models.  Working Paper.
%     Journal of Statistical Software.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfX_FisherSnedecorF(t,df1,df2,tol);

%% CHECK THE INPUT PARAMETERS
if nargin < 1
    error(message('VW:cfX_FisherSnedecorF:TooFewInputs'));
end
if nargin < 2, df1 = []; end
if nargin < 3, df2 = []; end
if nargin < 4, tol = []; end

if isempty(df1), df1 = 1; end
if isempty(df2), df2 = 1; end
if isempty(tol), tol = 1e-6; end
reltol = tol;

%% WARNING: The result is not reliable for large and comparable df1 and df2
if ((df1+df2)>150 && min(df1,df2)>50)
    warning('cfF:LossOfPrecision',...
        ['CF is not reliable for LARGE DF1 = ',num2str(df1), ...
        ' and LARGE DF2 = ',num2str(df2)]);
end

%% EVALUATE THE CHARACTERISTIC FUNCTION: cfF(df1,df2,t)

sz = size(t);
t  = t(:);
cf = ones(size(t));
id = t~=0;
cf(id) = integral(@(x) bsxfun(@times,funCF(df1,df2,t(id),(x/(1-x))^2), ...
    2*x/(1-x)^3),0,1,'ArrayValued',true,'RelTol',reltol);
cf = reshape(cf,sz);

end

%% Function funCF
function f = funCF(df1,df2,t,x)
%FUNCF Integrand function of the integral representation of the
%  characteristic function CF of the F-distribution with df1 and df2
%  degrees of freedom at the real argument t.
%
%  The characteristic function of F-distribution with df1 and df2 degrees 
%  of freedom, evaluated at the real t from (-inf,+inf), is defined by
%   CF(t) = (gamma(df1/2+df2/2)/gamma(df2/2)) * ...
%           HypergeometricU(df1/2, 1-df2/2, -1i*(df2/df1)*t),
%  where HypergeometricU(a,b,z) denotes the confluent hypergeometric
%  function of the second kind: U(a,b,z). 
%  
%  Here we use an integral representation of the hypergeometric function 
%  U(a,b,z), defined for purly complex argument z as (VW2016):
%   U(a,b,z) = gamma(1-b)/gamma(a-b+1) * Integral_0^inf cfFun(a,b,x) dx. 
%
% SYNTAX:
%   f = funCF(df1,df2,t,x)
%
% EXAMPLE
%  df1 = 5;
%  df2 = 4;
%  t   = 1:5;
%  x   = linspace(0,1);
%  f   = funCF(df1,df2,t,x);
%  plot(x,real(f),x,imag(f))

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Apr-2016 13:40:09

%% ALGORITHM FUNCF
t  = t(:);
nt = length(t);
o1 = ones(nt,1);
x  = x(:)';
nx = length(x);
o2 = ones(1,nx);
a  = df1/2;
b  = 1 - df2/2;
c  = (gammaln(a-b+1) - gammaln(1-b) - gammaln(a));
z  = -(df2/df1)*t;
f  = zeros(nt,nx);

% z == 0
id = (z==0);
if any(id)
    f(id,:) = exp(c + (a-1).*log(x) + (b-a-1).*log(1+x));
end

% abs(z) >= 1
id = (abs(z)>=1);
if any(id)
    zi = -1i./z(id);
    f(id,:) = exp(c + log(zi*o2) + (a-1)*log(zi*x) + ...
        (b-a-1)*log(1+zi*x) - o1(id)*x);
end

% z > 0 & z < 1
id = (z>0 & z<1);
if any(id)
    f(id,:) = exp(c + log(-1i) + (a-1)*log(-1i*o1(id)*x) + ...
        (b-a-1)*log(1-1i*o1(id)*x) - (z(id))*x);
end

% z < 0 & z > -1
id = (z<0 & z>-1);
if any(id)
    f(id,:) = exp(c + log(1i) + (a-1)*log(1i*o1(id)*x) + ...
        (b-a-1)*log(1+1i*o1(id)*x) + (z(id))*x);
end

end