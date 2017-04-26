function cf = cfE_Empirical(t,data,cfX)
%cfE_Empirical(t,coefs,cfX) evaluates the characteristic function cf(t) of
%  the Empirical distribution, based on the observed data. In particular,
%    cf(t) = cfE_Empirical(t,data) 
%          = (1/n) * sum_{j=1}^n cf_Dirac(data(j)*t), 
%  where cf_Dirac(t) represents the CF of the Dirac random variable
%  concentrated at the constant 1.
%
%  cfE_Empirical(t,Data,cf_X) evaluates the compound characteristic
%  function
%   cf(t) = cfE_Empirical(-1i*log(cfX(t)),data)
%         = (1/n) * sum_{j=1}^n cfX(t)^data(j);
% where cfX is function handle of the characteristic function cfX(t) of the
% random variable X (as e.g. another empirical CF based on observed data of
% X).   
%
% SYNTAX
%  cf = cfE_Empirical(t,data)
%  cf = cfE_Empirical(t,data,cfX)
%
% EXAMPLE1 (Empirical CF - a weighted mixture of independent Dirac variables)
%  rng(101);
%  n = 1000;
%  data = [normrnd(5,0.2,3*n,1); trnd(3,n,1); chi2rnd(1,n,1)];
%  t = linspace(-50,50,2^10);
%  cf = cfE_Empirical(t,data);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('Empirical CF - CF of the mixture of Dirac random variables')
%
% EXAMPLE2 (Convolution of the ECF and the Gaussian kernel)
%  rng(101);
%  n = 1000;
%  data = [normrnd(5,0.2,3*n,1); trnd(3,n,1); chi2rnd(1,n,1)];
%  bandwidth = 0.25;
%  cf_DATA   = @(t) cfE_Empirical(t,data)
%  cf_KERNEL = @(t) exp(-(bandwidth*t).^2/2);
%  cf = @(t) cf_DATA(t) .* cf_KERNEL(t);
%  t = linspace(-50,50,2^10);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('Smoothed Empirical CF')
%  result = cf2DistGP(cf)
%
% EXAMPLE3 (PDF/CDF of the compound Empirical-Empirical distribution)
%  rng(101);
%  lambda = 25; nN = 10; Ndata = poissrnd(lambda,1,nN);
%  mu = 0.1; sigma = 2; nX = 1500; Xdata = lognrnd(mu,sigma,1,nX);
%  cfX = @(t) cfE_Empirical(t,Xdata);
%  cf  = @(t) cfE_Empirical(t,Ndata,cfX);
%  t = linspace(-0.2,0.2,2^10);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('Compound Empirical CF')
%  x = linspace(0,1000,501);
%  prob = [0.9 0.95];
%  clear options
%  options.N = 2^12;
%  options.xMin = 0;
%  options.SixSigmaRule = 10;
%  

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
%cf = cfE_Empirical(t,data,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 2, data = []; end
if nargin < 3, cfX = []; end

%%
if isempty(data)
    data = 1;
end

weights = 1 / length(data);
data    = data(:)';

% Special treatment for mixtures with large number of variables
szcoefs  = size(data);
szcoefs  = szcoefs(1)*szcoefs(2);
szt      = size(t);
sz       = szt(1)*szt(2);
szcLimit = ceil(1e3 / (sz/2^16));
idc = 1:fix(szcoefs/szcLimit)+1;

%% Characteristic function of a weighted mixture of Dirac variables
t     = t(:);
idx0  = 1;
cf    = 0;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    if isempty(cfX)
        aux  = exp(1i * t * data(idx));
    else
        aux  = bsxfun(@power,cfX(t),data(idx));
    end      
    if isscalar(weights)
        cf   = cf + sum(weights*aux,2);
    else
        cf   = cf + sum(bsxfun(@times,aux,weights(idx)),2);
    end
end
cf = reshape(cf,szt);

end