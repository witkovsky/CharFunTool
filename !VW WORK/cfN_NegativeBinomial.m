function cf = cfN_NegativeBinomial(t,r,p,cfX)
%cfN_NegativeBinomial(t,r,p) evaluates the characteristic function cf(t) of
% the Negative-Binomial distribution with the parameters r ( number of
% failures until the experiment is stopped, r in N) and p (success
% probability in each experiment, p in [0,1]), i.e.  
%   cf(t) = cfN_NegativeBinomial(t,r,p) 
%         = p^r * (1 - (1-p) * e^(1i*t))^(-r);
% For more details see [4], and also WIKIPEDIA:
% https://en.wikipedia.org/wiki/Negative_binomial_distribution
%
% SYNTAX
%  cf = cfN_NegativeBinomial(t,r,p)
%  cf = cfN_NegativeBinomial(t,r,p,cfX)
%
% EXAMPLE1 (CF of the Negative Binomial distribution with r = 5, p = 0.3)
%  r = 5;  
%  p = 0.3;
%  t = linspace(-15,15,1001);
%  cf = cfN_NegativeBinomial(t,r,p);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('CF the Negative Binomial distribution with r = 5, p = 0.3')
%
% EXAMPLE2 (PDF/CDF of the compound NegativeBinomial-Exponential distribution)
%  r = 5;  
%  p = 0.3;
%  lambda = 5;
%  cfX = @(t) cfX_Exponential(t,lambda);
%  cf = @(t) cfN_NegativeBinomial(t,r,p,cfX);
%  x = linspace(0,10,101);
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
%cf = cfN_NegativeBinomial(t,r,p,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, r = []; end
if nargin < 3, p = []; end
if nargin < 4, cfX = []; end

%%
if isempty(r)
    r = 10;
end

if isempty(p)
    p = 0.5;
end
%% Characteristic function of the (compound) Negative-Binomial distribution
szt = size(t);
t   = t(:);

if isempty(cfX)
    expit = exp(1i*t);
else
    expit = cfX(t);
end

cf  = p.^r .* (1 - (1-p).*expit).^(-r);
cf = reshape(cf,szt);
cf(t==0) = 1;

end