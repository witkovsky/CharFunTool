function cf = cfE_Empirical(t,data,cfX)
%% cfE_Empirical
%  Characteristic function of the EMPIRICAL distribution, based on the
%  observed data. 
%
%  That is, cf(t) = (1/N)*(cfD(data_1*t) +...+ cfD(data_N*t)), where
%  cfD(t) represents the characteristic function of the DIRAC RV
%  concentrated at the constant d=1, i.e. cfD(t) = exp(1i*t). 
%
%  cfE_Empirical(t,data,cfX) evaluates the compound characteristic function
%   cf(t) = cfE_Empirical(-1i*log(cfX(t)),data)
%         = (1/N) * sum_{j=1}^N cfX(t)^data_j;
%  where cfX is function handle of the characteristic function cfX(t) of the
%  random variable X.   
%
% SYNTAX
%  cf = cfE_Empirical(t,data)
%  cf = cfE_Empirical(t,data,cfX)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  data   - vector of data, i.e. constants where the DIRAC RVs are
%           concentrated. If empty, default value is data = 1. 
%  cfX    - function handle of the characteristic function of a random
%           variable X. If cfX is non-empty, a compound CF is evaluated as
%           cf(t) = cf(-1i*log(cfX(t)),data).
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
%
% EXAMPLE 1 (Empirical CF - a weighted mixture of independent Dirac variables)
%  rng(101);
%  n = 1000;
%  data = [normrnd(5,0.2,3*n,1); trnd(3,n,1); chi2rnd(1,n,1)];
%  t = linspace(-50,50,2^10);
%  cf = cfE_Empirical(t,data);
%  figure; plot(t,real(cf),t,imag(cf)),grid
%  title('Empirical CF - CF of the mixture of Dirac random variables')
%
% EXAMPLE 2 (Convolution of the ECF and the Gaussian kernel)
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
% EXAMPLE 3 (PDF/CDF of the compound Empirical-Empirical distribution)
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
%  result = cf2DistGP(cf,x,prob,options)
%  
% REFERENCES:
%  WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate
%  loss distribution based on numerical inversion of the compound empirical
%  characteristic function of frequency and severity. arXiv preprint
%  arXiv:1701.08299.   

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Sep-2018 12:56:00
%% ALGORITHM
%cf = cfE_Empirical(t,data,cfX);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 2, data = []; end
if nargin < 3, cfX = []; end

%%

weights = [];
cf = cfE_DiracMixture(t,data,weights,cfX);

end