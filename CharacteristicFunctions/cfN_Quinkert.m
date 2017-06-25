function cf = cfN_Quinkert(t,a,b,cfX)
%cfN_Quinkert(t,a,b) evaluates the characteristic function cf(t) of the
% Quinkert distribution, with the parameters a (a > 0) and b (b > 0), i.e. 
%   cf(t) = cfN_Quinkert(t,a,b) = 1F1(a,a+b,e^(1i*t)-1);
% where 1F1 denotes the confluent hypergeometric (Kummer's) function. For
% more details see [4], p. 564.
%
% SYNTAX
%  cf = cfN_Quinkert(t,a,b)
%  cf = cfN_Quinkert(t,a,b,cfX)
%
% EXAMPLE1 (CF of the Quinkert distribution with the parameter a=3, b=5)
%  a = 3;  
%  b = 5;
%  t = linspace(-15,15,501);
%  cf = cfN_Quinkert(t,a,b);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the Quinkert distribution with a = 3, b = 5')
%
% EXAMPLE2 (CF of the compound Quinkert-Exponential distribution)
%  a = 3;  
%  b = 5;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Quinkert(t,a,b,cfX)
%  t = linspace(-15,15,501);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('CF of the compound Quinkert-Exponential distribution')
%
% EXAMPLE3 (PDF/CDF of the compound Quinkert-Exponential distribution)
%  a = 3;  
%  b = 5;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Quinkert(t,a,b,cfX);
%  x = linspace(0,1.5,101);
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
%cf = cfN_Quinkert(t,a,b,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, a = []; end
if nargin < 3, b = []; end
if nargin < 4, cfX = []; end

%%
if isempty(a), a = 1; end
if isempty(b), b = 1; end

%% Characteristic function of the Quinkert distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf = Hypergeom1F1(a,a+b,expit-1);
cf = reshape(cf,szt);
cf(t==0) = 1;

end
