function cf = cf_RandomMean(t,cfX,n,p)
%% cf_RandomMean
%  Characteristic function of the RANDOM MEAN (RM) Y = 1/N Sum_(k=1)^N X_k  
%  of N (random number) i.i.d. random variables X_k, with their common
%  distribution specified by the characteristic function cfX. Probability
%  distribution of the random number N is given by its PMF, here specified
%  by the K-dimensional vector n (the vector of possible observed values of
%  N, say n = (n_1,...,n_K)) and by the vector p of its probabilities, say
%  p = (p_1,...,p_K), such that sum_{k=1}^K p_k = 1.
%
%  Probability distribution of the RANDOM MEAN is specified as a mixture of
%  n-times convolved distributions of X/n. That is, its CF of Y is given by
%   cf(t) = sum_{k=1}^K  p_k * cfX(t/n_k)^n_k.
%  For more details see also CHRISTOPH, MONAKHOV and ULYANOV (2018).
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
% REMARKS:
%  - If n is scalar (and p = 1), we get CF of the regular (non-random) mean
%    of n iid RVs. 
%  - The algorithm allows to consider also the situation when one of the
%    possible observed values of N is 0, i.e. with n_1 = 0 with p_1 > 0. By
%    definition we set cfX(t/n_1)^n_1 = 1 if n_1 = 0. That is, proportion
%    p_1 of the distribution is concentrated at 0.   
%  - Notice that CF of the RANDOM SUM, Y = Sum_(k=1)^N X_k, is specified by
%    cf(t) = sum_{k=1}^K  p_k * cfX(t)^n_k. If CF of the discrete random
%    variable N is known as cfN(t), this can be simplified to the following
%    expression: cf(t) = cfN(-1i*log(cfX(t)). For more details see the
%    implementation of algorithms for computing CFs of the discrete
%    distributions.
%
% EXAMPLE1 (CF of RANDOM MEAN of N chi-squared RVs, N~Bino(10,0.5))
%  df    = 1;
%  cfX    = @(t) cf_ChiSquare(t,df);
%  nBino  = 10;
%  pBino  = 0.5;
%  n      = 0:nBino;
%  p      = binopdf(n,nBino,pBino);
%  cf     = @(t) cf_RandomMean(t,cfX,n,p);
%  t      = linspace(-20,20,501);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('CF of the RANDOM MEAN of chi-squared RVs with N~Bino(10,0.5)')
%
% EXAMPLE2 (PDF/CDF/QF of RANDOM MEAN of N chi-squared RVs, N~Bino(10,0.5))
%  df     = 1;
%  cfX    = @(t) cf_ChiSquare(t,df);
%  nBino  = 10;
%  pBino  = 0.5;
%  n      = 0:nBino;
%  p      = binopdf(n,nBino,pBino);
%  cf     = @(t) cf_RandomMean(t,cfX,n,p);
%  x      = linspace(0,5)';
%  prob   = [0.9 0.95 0.99];
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF/QF of Z = sqrt(E(N))*RandomMean((X-E(X))/STD(X)))
%  df     = 1;
%  cfX    = @(t) cf_ChiSquare(t/sqrt(2*df),df) .* exp(-1i*t/sqrt(2));
%  nBino  = 10;
%  pBino  = 0.5;
%  n      = 0:nBino;
%  p      = binopdf(n,nBino,pBino);
%  EN     = nBino * pBino;
%  cf     = @(t) cf_RandomMean(sqrt(EN)*t,cfX,n,p);
%  result = cf2DistGP(cf,[],prob)
%
% EXAMPLE4 (PDF/CDF/QF of the RANDOM MEAN of N chi-squared RVs)
%  % See also CHRISTOPH, MONAKHOV and ULYANOV (2018)
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
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES: 
% [1] CHRISTOPH, G., MONAKHOV  M.M., and ULYANOV, V.V. (2018). Second Order
%     Expansions for Distributions of Statistics and Its Quantiles Based on
%     Random Size Samples. Preprint March 2018. 
%     DOI:10.13140/RG.2.2.28836.17281.  
%     http://www.math.sci.hiroshima-u.ac.jp/stat/TR/TR17/TR17-15.pdf   
%
% SEE ALSO: cf_RandomSum

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Jun-2018 10:16:08

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
    if n(k)~=0
        cf = cf + p(k) * cfX(t/n(k)).^n(k);
    else
        cf = cf + p(k);
    end
end

cf = reshape(cf,szt);
cf(t==0) = 1;

end