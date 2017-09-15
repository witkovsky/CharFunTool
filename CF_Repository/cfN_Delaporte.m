function cf = cfN_Delaporte(t,a,b,c,cfX)
%cfN_Delaporte(t,a,b,c)  evaluates the characteristic function cf(t) of the
% Delaporte distribution with the parameters a (parameter of variable mean,
% a > 0), b (parameter of variable mean, b > 0 ), and c (fixed mean, c >
% 0), i.e.  
%   cf(t) = cfN_Binomial(t,n,p) = (1 - p + p*e^(1i*t))^n
%         = (b/(1+b))^a * (1-e^(1i*t)/(b+1))^(-a) * exp(-c*(1-e^(1i*t)));
% For more details see [4], and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Delaporte_distribution
%
% SYNTAX
%  cf = cfN_Delaporte(t,a,b,c)
%  cf = cfN_Delaporte(t,a,b,c,cfX)
%
% EXAMPLE1 (CF of the Delaporte distribution with a=2.2, b=3.3, c=4)
%  a = 2.2;
%  b = 3.3;
%  c = 4; 
%  t = linspace(-15,15,1001);
%  cf = cfN_Delaporte(t,a,b,c);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Delaporte distribution with a=2.2, b=3.3, c=4')
%
% EXAMPLE2 (CF of the compound Delaport-Exponential distribution)
%  a = 2.2;
%  b = 3.3;
%  c = 4; 
%  lambda = 10;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  t = linspace(-10,10,501);
%  cf = cfN_Delaporte(t,a,b,c,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the compound Delaport-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Delaport-Exponential distribution)
%  a = 2.2;
%  b = 3.3;
%  c = 4; 
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Delaporte(t,a,b,c,cfX);
%  x = linspace(0,4,101);
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
%cf = cfN_Delaporte(t,a,b,c,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(4, 5);
if nargin < 5, cfX = []; end

%% Characteristic function of the (compound) Delaporte distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf = (b/(1+b))^a * (1-expit/(b+1)).^(-a) .* exp(-c*(1-expit));
cf = reshape(cf,szt);
cf(t==0) = 1;

end