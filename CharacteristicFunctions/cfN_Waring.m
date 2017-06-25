function cf = cfN_Waring(t,a,b,r,cfX)
%cfN_Waring(t,a,b,r) evaluates the characteristic function cf(t) of the
% Waring distribution, with the parameters a (a > 0), b (b > 0), and r
% (r > 0), i.e.
%   cf(t) = cfN_Waring(t,a,b,r)
%         = ((gamma(a+r)*gamma(a+b)) / (gamma(a)*gamma(a+b+r)))
%           * 2F1(r,b,a+b+r,e^(1i*t));
% where 2F1 denotes the Gauss hypergeometric function. The Waring
% distribution is also known as beta negative binomial distribution. For
% more details see [4], p. 643, and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Beta_negative_binomial_distribution
%
% SYNTAX
%  cf = cfN_Waring(t,a,b,r)
%  cf = cfN_Waring(t,a,b,m,cfX)
%
% EXAMPLE1 (CF of the Waring distribution with a = 2.2, b = 3.3, r = 4)
%  % The CF is not computed correctly!! Because the hypergeomtric function
%  % Hypergeom2F1(r,b,a+b+r,z) does not converege abs(z)>=1. Here z =
%  % exp(1i*t), ans abs(exp(1i*t)) = 1.
%  a = 2.2;
%  b = 3.3;
%  r = 4;
%  t = linspace(-5,5,1001);
%  cf = cfN_Waring(t,a,b,r);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('(CF of the Waring distribution with a = 2.2, b = 3.3, r = 4')
%
% EXAMPLE2 (CF of the compound Waring-Exponential distribution)
%  a = 2.2;
%  b = 3.3;
%  r = 4;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  t = linspace(-10,10,501);
%  cf = cfN_Waring(t,a,b,r,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the compound Waring-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Waring-Exponential distribution)
%  a = 2.2;
%  b = 3.3;
%  r = 4;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Waring(t,a,b,r,cfX);
%  x = linspace(0,35,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Preprint submitted
%     to Insurance: Mathematics and Economics.
% [2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models. Preprint submitted
%     to Journal of Statistical Software.
% [3] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
% [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
%     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN
%     3-87773-025-6.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfN_Waring(t,a,b,r,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(4, 5);
if nargin < 5, cfX = []; end

%% Characteristic function of the (compound) Waring distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf    = exp(gammaln(a+r) + gammaln(a+b) - gammaln(a) - gammaln(a+b+r));
cf    = cf * Hypergeom2F1(r,b,a+b+r,expit);
cf    = reshape(cf,szt);
cf(t==0) = 1;

end