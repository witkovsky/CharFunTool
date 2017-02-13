function cf = cfN_Geometric(t,p,type,cfX)
%cf = cfN_Geometric(t,p,type,cfX) evaluates the characteristic function
% cf(t) of the Geometric distribution. 
% The standard Geometric distribution (type = 'standard' or 'zero') is
% defined on non-negative integers k = 0,1,... . The shifted Geometric
% distribution (type = 'shifted') is defined on positive integers k =
% 1,2,... . Both types are parametrized by the success probability
% parameter p in [0,1]), i.e. the 'standard' CF is defined by 
%  cf(t) = cfN_Geometric(t,p,'standard')
%        = p ./ (1 - (1-p) * exp(1i*t)),
% and the 'shifted' CF is defined by
%   cf(t) = cfN_Geometric(t,p,'shifted') 
%         = exp(1i*t) .* (p ./ (1 - (1-p) * exp(1i*t))).
% For more details see [4], and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Geometric_distribution.
%
% cfN_Geometric(t,p,type,cfX) evaluates the compound characteristic
% function 
%   cf(t) = cfN_Geometric(-1i*log(cfX(t)),p),
% where cfX is function handle of the characteristic function cfX(t) of a
% continuous distribution and/or random variable X.  
% 
% Note that such CF is characteristic function of the compound distribution,
% i.e. distribution of the random variable Y = X_1 + ... + X_N, where X_i ~
% F_X are i.i.d. random variables with common CF cfX(t), and N ~ F_N is
% independent RV with its CF given by cfN(t). 
% 
% SYNTAX
%  cf = cfN_Geometric(t,p)
%  cf = cfN_Geometric(t,p,type,cfX)
%
% EXAMPLE1 (CF of the Geometric distribution with the parameter p = 0.5)
%  p = 0.5;  
%  t = linspace(-10,10,501);
%  cf = cfN_Geometric(t,p);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('Characteristic function of the Geometric distribution')
%
% EXAMPLE2 (CF of the Geometric distribution with the parameter p = 0.5)
%  p = 0.5;  
%  type = 'shifted';
%  t = linspace(-10,10,501);
%  cf = cfN_Geometric(t,p,type);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('Characteristic function of the shifted Geometric distribution')
%
% EXAMPLE3 (CF of the compound Geometric-Exponential distribution)
%  p = 0.5;
%  type = 'standard';
%  lambda = 5;
%  t = linspace(-10,10,501);
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = cfN_Geometric(t,p,type,cfX);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF of the compound Geometric-Exponential distribution')
%
% EXAMPLE4 (PDF/CDF of the compound Geometric-Exponential distribution)
%  p = 0.5;
%  type = 'standard';
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_Geometric(t,p,type,cfX);
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
%     output quantity in linear measurement models. Acta IMEKO, 5(3),32-44.  
% [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
%     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN
%     3-87773-025-6. 

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfN_Geometric(t,p,type,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, p = []; end
if nargin < 3, type = []; end
if nargin < 4, cfX = []; end

%%
if isempty(p)
    p = 1;
end

if isempty(type)
    type = 'standard';
end

%% Characteristic function of the (compound) Geometric distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

switch lower(type)
    case lower({'standard','zero'})
        cf  = 1;     
    case lower({'shifted'})
        cf  = expit;
    otherwise
        cf  = 1;
end

cf = cf .* (p ./ (1 - (1-p) * expit));
cf = reshape(cf,szt);
cf(t==0) = 1;

end