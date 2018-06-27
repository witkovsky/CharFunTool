function cf = cf_RandomMean(t,cfX,n,p)
%% cf_RandomMean
%  Characteristic function of a RANDOM MEAN of N i.i.d. random variables X,
%  specified by the characteristic function cfX, where the random number N
%  is specified by the probability distribution (n is vector of possible
%  values and p is vector of its probabilities).
%
%  Distribution of RANDOM MEAN is specified as a mixture of the convolved
%  distributions of X_n = X/n +...+ X/n. That is, its CF is given by
%   cf(t) = sum_{n=1}^n_max  p_n * cfX(t/n)^n   
%  
% SEE ALSO: CHRISTOPH, MONAKHOV and ULYANOV (2018).
%
% SYNTAX:
%  cf = cf_RandomMean(t,cfX,n)
%  cf = cf_RandomMean(t,cfX,n,p)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  cfX    - function handle of the characteristic function of a random
%           variable X. 
%  n      - vector or possible integer values n = 0,1,2,... of the random
%           variable N. 
%  p      - vector or probabilities (with sum(p) = 1) of the values
%           specified by n. If empty, default value is p = 1/length(n) for
%           all n.
% 
% REMARK:
%  If n is scalar (and p = 1), we get CF of the regular (non-random) mean
%  of n iid RVs. 
%
% EXAMPLE1 (CF of RANDOM MEAN of N chi-squared RVs, N~Bino(10,0.5))
%  df   = 1;
%  cfX  = @(t) cf_ChiSquare(t,df);
%  N    = 10;
%  pr   = 0.5;
%  n    = 0:N;
%  p    = binopdf(n,N,pr);
%  cf   = @(t) cf_RandomMean(t,cfX,n,p);
%  t = linspace(-20,20,501);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('CF of the RANDOM MEAN of chi-squared RVs with N~Bino(10,0.5)')
%
% EXAMPLE2 (PDF/CDF/QF of RANDOM MEAN of N chi-squared RVs, N~Bino(10,0.5))
%  df   = 1;
%  cfX  = @(t) cf_ChiSquare(t,df);
%  N    = 10;
%  pr   = 0.5;
%  n    = 0:N;
%  p    = binopdf(n,N,pr);
%  cf   = @(t) cf_RandomMean(t,cfX,n,p);
%  x    = linspace(0,5)';
%  prob = [0.9 0.95 0.99];
%  options.xMin = 0;
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE2 (Continued)
% % PDF/CDF/QF of the standardized RANDOM MEAN Z = sqrt(m)*(X-E(X))/STD(X)
%  EX   = result.xMean;
%  STDX = result.xStd;
%  m    = N*pr;
%  cfZ  = @(t) cf(t) .* exp(-1i*t*EX);
%  cfZ  = @(t) cfZ(sqrt(m)*t/STDX);
%  resultZ = cf2DistGP(cfZ,[],prob)
%
% EXAMPLE3 (See also CHRISTOPH, MONAKHOV and ULYANOV (2018))
% % PDF/CDF/QF of the RANDOM MEAN of N chi-squared RVs
%  df   = 1;
%  cfX  = @(t) cf_ChiSquare(t,df);
%  P    = @(k,s,n) (k./(s+k)).^n - ((k-1)./(s+k-1)).^n;
%  N    = 10;
%  S    = 2;
%  kMax = 1000; % This should LARGE!, e.g. kMax = 10000
%  K    = (1:kMax)';
%  pr   = P(K,S,N);
%  n    = [0;K];
%  p    = [1-sum(pr);pr];
%  cf   = @(t) cf_RandomMean(t,cfX,n,p);
%  x    = linspace(0,3,2^10)';
%  prob = [0.9 0.95 0.99];
%  options.xMin = 0;
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES: 
% [1] CHRISTOPH, G., MONAKHOV  M.M., and ULYANOV, V.V. (2018). Second Order
%     Expansions for Distributions of Statistics and Its Quantiles Based on
%     Random Size Samples. Preprint March 2018. 
%     DOI:10.13140/RG.2.2.28836.17281.  
%     http://www.math.sci.hiroshima-u.ac.jp/stat/TR/TR17/TR17-15.pdf   

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 7-Jun-2018 16:14:36

%% ALGORITHM
%cf = cf_RandomMean(t,cfX,n,p)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, p   = []; end
if nargin < 3, n   = []; end
if nargin < 2, cfX = []; end

if isempty(cfX)
    cfX = @(t) cfS_Gaussian(t);
end

if isempty(n), n = 1; end
if isempty(p), p = 1/length(n); end

%% Characteristic function
szt   = size(t);
t     = t(:);
cf    = 0;

for k = 1:length(n)
    cf = cf + p(k) * cfX(t/n(k)).^n(k);
end

cf = reshape(cf,szt);
cf(t==0) = 1;

end