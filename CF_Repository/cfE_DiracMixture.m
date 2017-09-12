function cf = cfE_DiracMixture(t,d,weight,cfX)
%% cfE_DiracMixture
%  Characteristic function of the weighted mixture distribution of
%  independent DIRAC random variables D_1,...,D_N, concentrated at the
%  fixed constants (data) given by the vector d = [d_1,...,d_N].
%  
%  That is, cf(t) = weight_1*cfD(d_1*t) +...+ weight_N*cfD(d_N*t), where
%  cfD(t) represents the characteristic function of the DIRAC RV
%  concentrated at the constant d=1, i.e. cfD(t) = exp(1i*t).
%
%  cfE_DiracMixture(t,d,weight,cfX) evaluates the compound characteristic
%  function  
%   cf(t) = cfE_DiracMixture(-1i*log(cfX(t)),d,weight)
%         = weight_1*cfX(t)^d_1 +...+ weight_N*cfX(t)^d_N
%  where cfX denotes the function handle of the characteristic function
%  cfX(t) of the random variable X.   
%
% SYNTAX
%  cf = cfE_DiracMixture(t,d,weight)
%  cf = cfE_DiracMixture(t,d,weight,cfX)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  d      - vector of constants (data) where the DIRAC RVs are concentrated.
%           If empty, default value is d = 1. 
%  weight - vector of weights of the distribution mixture. If empty,
%           default value is weight = 1/length(d). 
%  cfX    - function handle of the characteristic function of a random
%           variable X. If cfX is non-empty, a compound CF is evaluated as
%           cf(t) = cf(-1i*log(cfX(t)),d,weight).
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
%
% EXAMPLE1 (Empirical CF - a weighted mixture of independent Dirac variables)
%  rng(101);
%  n = 1000;
%  data = [normrnd(5,0.2,3*n,1); trnd(3,n,1); chi2rnd(1,n,1)];
%  t = linspace(-50,50,2^10);
%  weights = 1/length(data);
%  cf = cfE_DiracMixture(t,data,weights);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('Empirical CF - CF of the mixture of Dirac random variables')
%
% EXAMPLE2 (Convolution of the ECF and the Gaussian kernel)
%  rng(101);
%  n = 1000;
%  data = [normrnd(5,0.2,3*n,1); trnd(3,n,1); chi2rnd(1,n,1)];
%  bandwidth = 0.25;
%  cf_DATA   = @(t) cfE_DiracMixture(t,data,weights)
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
%  cfX = @(t) cfE_DiracMixture(t,Xdata,1/nX);
%  cf  = @(t) cfE_DiracMixture(t,Ndata,1/nN,cfX);
%  t = linspace(-0.2,0.2,2^10);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('Compound Empirical CF')
%  x = linspace(0,1000,501);
%  prob = [0.9 0.95];
%  clear options
%  options.N = 2^12;
%  options.xMin = 0;
%  options.SixSigmaRule = 10;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
%  WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate
%  loss distribution based on numerical inversion of the compound empirical
%  characteristic function of frequency and severity. arXiv preprint
%  arXiv:1701.08299.   

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Jun-2017 10:00:49

%% ALGORITHM
%cf = cfE_DiracMixture(t,d,weight,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, d = []; end
if nargin < 3, weight = []; end
if nargin < 4, cfX = []; end

%%
if isempty(d)
    d = 1;
end

if isempty(weight)
    weight = 1 / length(d);
end

if isscalar(weight)
    d = d(:)';
else
    [errorcode,d,weight] = distchck(2,d(:)',weight(:)');
    if errorcode > 0
        error(message('InputSizeMismatch'));
    end
end

% Special treatment for mixtures with large number of variables
szcoefs  = size(d);
szcoefs  = szcoefs(1)*szcoefs(2);
szt      = size(t);
sz       = szt(1)*szt(2);
szcLimit = ceil(1e3 / (sz/2^16));
idc = 1:fix(szcoefs/szcLimit)+1;

%% Characteristic function
t     = t(:);
idx0  = 1;
cf    = 0;
for j = 1:idc(end)
    idx1 = min(idc(j)*szcLimit,szcoefs);
    idx  = idx0:idx1;
    idx0 = idx1+1;
    if isempty(cfX)
        aux  = exp(1i * t * d(idx));
    else
        aux  = bsxfun(@power,cfX(t),d(idx));
    end      
    if isscalar(weight)
        cf   = cf + sum(weight*aux,2);
    else
        cf   = cf + sum(bsxfun(@times,aux,weight(idx)),2);
    end
end
cf = reshape(cf,szt);

end