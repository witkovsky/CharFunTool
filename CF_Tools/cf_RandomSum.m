function cf = cf_RandomSum(t,cfX,n,p)
%% cf_RandomSum
%  Characteristic function of the RANDOM SUM (RM) Y = Sum_(k=1)^N X_k  
%  of N (random number) i.i.d. random variables X_k, with their common
%  distribution specified by the characteristic function cfX. Probability
%  distribution of the random number N is given by its PMF, here specified
%  by the K-dimensional vector n (the vector of possible observed values of
%  N, say n = (n_1,...,n_K)) and by the vector p of its probabilities, say
%  p = (p_1,...,p_K), such that sum_{k=1}^K p_k = 1.
%
%  Probability distribution of the RANDOM SUM is specified as a mixture of
%  n-times convolved distributions of X. That is, its CF of Y is given by
%   cf(t) = sum_{k=1}^K  p_k * cfX(t)^n_k.
%
% SYNTAX:
%  cf = cf_RandomSum(t,cfX,n)
%  cf = cf_RandomSum(t,cfX,n,p)
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
%  - The algorithm allows to consider also the situation when one of the
%    possible observed values of N is 0, i.e. with n_1 = 0 with p_1 > 0. By
%    definition we set cfX(t)^n_1 = 1 if n_1 = 0. That is, proportion
%    p_1 of the distribution is concentrated at 0.   
%  - Notice that if CF of the discrete random variable N is known as
%    cfN(t), this can be simplified to the following expression: cf(t) =
%    cfN(-1i*log(cfX(t)). For more details see the implementation of
%    algorithms for computing CFs of the discrete distributions.
%
% EXAMPLE1 (PDF/CDF/QF of RANDOM MEAN of N chi-squared RVs, N~Bino(10,0.5))
%  % Evaluation by using the function RandomSum
%  df     = 1;
%  cfX    = @(t) cf_ChiSquare(t,df);
%  nBino  = 10;
%  pBino  = 0.5;
%  n      = 0:nBino;
%  p      = binopdf(n,nBino,pBino);
%  cf     = @(t) cf_RandomSum(t,cfX,n,p);
%  x      = linspace(0,30)';
%  prob   = [0.9 0.95 0.99];
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE2 (PDF/CDF/QF of RANDOM MEAN of N chi-squared RVs, N~Bino(10,0.5))
%  % Alternative evaluation by using the compound CF function cfN_Binomial
%  df     = 1;
%  cfX    = @(t) cf_ChiSquare(t,df);
%  nBino  = 10;
%  pBino  = 0.5;
%  cf     = @(t) cfN_Binomial(t,nBino,pBino,cfX);
%  x      = linspace(0,30)';
%  prob   = [0.9 0.95 0.99];
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF/QF of the RANDOM SUM of N chi-squared RVs)
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
%  cf   = @(t) cf_RandomSum(t,cfX,n,p);
%  x    = linspace(0,800,2^10)';
%  prob = [0.9 0.95 0.99];
%  options.xMin = 0;
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% SEE ALSO: cf_RandomMean

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Jun-2018 10:10:45

%% ALGORITHM
%cf = cf_RandomSum(t,cfX,n,p)

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
        cf = cf + p(k) * cfX(t).^n(k);
    else
        cf = cf + p(k);
    end
end

cf = reshape(cf,szt);
cf(t==0) = 1;

end