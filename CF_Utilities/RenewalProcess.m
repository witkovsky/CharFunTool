function [pmf,n,tMin,tMax,cdfFun,pmfFun] = ... 
    RenewalProcess(t,cf,nMin,nMax,tol)
%RenewalProcess 
%  Evaluates the renewal process distribution function, say PMF(t),
%  specified by the characteristic function (CF) of the holding time
%  distribution.
%
%  For each specified time t the renewal process distribution function
%  pmf(t) gives the probability mass function defined on the integers from
%  the specified interval(nMin, nMax). The pmf(t) can be evaluated at
%  arbitrary time t >= 0. However, for correct result t must be from the
%  associated interval (tMin,tMax), which can be specified by the
%  RenewalProcess algorithm.
%
%  If time t is specified as empty set, t = [], then the RenewalProcess
%  returns pmf as anonymous function of time t, i.e. pmf = @(t) pmf(t),
%  which can be evaluated for arbitrary time t.
%
% SYNTAX:
%  pmf = RenewalProcess(t,cf)
%  [pmf,n] = RenewalProcess(t,cf,nMin,nMax,tol)
%  [pmf,n,tMin,tMax,cdfFun,pmfFun] = RenewalProcess(t,cf,nMin,nMax,tol)
% 
% INPUTS:
%  t      - scalar value of time t >= 0, where the probability mass
%           function pmf(t) is evaluated. If empty, t = [], the
%           RenewalProcess returns pmf as anonymous function of time t,
%           pmf = @(t) pmf(t);
%  cf     - function handle of the characteristic function (CF) of the
%           holding time distribution (a continuous non-negative
%           distribution). If empty, cf = [], the default value is the
%           characteristic function of the standard exponential
%           distribution, i.e. cf = @(t)cf_Exponential(t);
%  nMin   - minimum value of the integers n, where the probability mass
%           function pmf(t) is calculated. If empty, nMin = [], the default
%           value is nMin = 0;
%  nMax   - maximum value of the integers n, where the probability mass
%           function pmf(t) is calculated. If empty, nMax = [], the default
%           value is nMax = 25;
%  tol    - tolerance factor for the numerical errors of the calculated 
%           probabilities, i.e. such that for Prob(n=k) < tol we set
%           Prob(n=k) = 0; If empty, default value is tol = 1e-6. 
% 
% OUTPUTS:
%  pmf    - vector of calculated probabilities evaluated for each integer n
%           in (nMin,nMax) at specified time t, i.e. the probability mass
%           function evaluated at specified time t, i.e. pmf = pmf(t).
%           Alternatively, if t = [], the result is an anonymous function
%           of time, i.e. pmf = @(t) pmf(t);
%  n      - vector of integers where the pmf(t) is evaluated;
%  tMin   - the suggested lower bound of time t, where pmf(t) can be safely
%           evaluated;
%  tMax   - the suggested upper bound of time t, where pmf(t) can be safely
%           evaluated;
%  cdfFun - auxiliary cell array of function handles of the n-times
%           convolved CDFs;
%  pmfFun - auxiliary cell array of function handles of the probabilities
%           for each specified integer n, expressed as a function of time
%           t. 
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Renewal_theory
%
% NOTES (from WIKIPEDIA)
%  A renewal process is a generalization of the Poisson process. In
%  essence, the Poisson process is a continuous-time Markov process on the
%  positive integers (usually starting at zero) which has independent
%  identically distributed holding times at each integer i (exponentially
%  distributed) before advancing (with probability 1) to the next integer:
%  i+1. In the same informal spirit, we may define a renewal process to be
%  the same thing, except that the holding times take on a more general
%  distribution. (Note however that the independence and identical
%  distribution (IID) property of the holding times is retained).
%
% EXAMPLE 1:
% % The renewal process with the gamma distributed holding times
%  alpha   =  2;
%  beta    =  5;
%  cf      = @(u)cf_Gamma(u,alpha,beta);
%  nMin    = 0;
%  nMax    = 70;
%  [pmf,n] = RenewalProcess([],cf,nMin,nMax);
%  plot(n,pmf(1),'o--',n,pmf(5),'o--',n,pmf(10),'o--',n,pmf(15),'o--')
%  title('Renewal process & gamma-distributed holding times t=[1,5,10,15]') 
%  xlabel('n')
%  ylabel('probability')

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Oct-2018 15:37:55

%% ALGORITHM
%  [pmf,n,tMin,tMax,cdfFun,pmfFun] = RenewalProcess(t,cf,nMin,nMax,tol)

%% CHECK THE INPUT PARAMETERS
narginchk(0, 5);
if nargin < 1, t    = []; end
if nargin < 2, cf  = []; end
if nargin < 3, nMin = []; end
if nargin < 4, nMax = []; end
if nargin < 5, tol  = []; end

if isempty(tol)
    tol = 1e-6;
end

if isempty(nMax)
    nMax = 25;
end

if isempty(nMin)
    nMin = 0;
end

if isempty(cf)
    cf = @(t) cf_Exponential(t);
end

n      = (nMin:nMax)';
N      = length(n);
cdfFun = cell(N+1,1);
pmfFun = cell(N,1);

options.xMin = 0;
options.isPlot = false;
options.isInterp = true;

tMin = inf;
tMax = 0;
if nargout > 2
    tprob = tol;
else
    tprob = [];
end
id   = 0;
for k = nMin:nMax
    id = id + 1;
    if k == 0
        cdfFun(id) = {@(xnew) 1};
    else
        result = cf2DistGP(@(t)cf(t).^k,[],tprob,options);
        cdfFun(id) = {result.CDF};
        tMinNew = result.xMin;
        tMaxNew = result.qf;
        options.xMin = tMinNew;
        tMin = min(tMin,tMinNew);
        tMax = max(tMax,tMaxNew);
    end
end
cdfFun(N+1) = {@(xnew) 0};

id   = 0;
for k = nMin:nMax
    id = id + 1;
    pmfFun(id) = {@(t)cdfFun{id}(t)-cdfFun{id+1}(t)};
end

if isscalar(t)
    pmf = pmfDist(t,pmfFun,N,tol);
else
    pmf = @(t) pmfDist(t,pmfFun,N,tol);
end

end
%% Function pmfDist
function pmf = pmfDist(t,pmfFun,n,tol)
%pmfDist 
%  Evaluates the PMF probabilities at the specified integers n, for given
%  tolerance tol
%
% SYNTAX:
%  pmf = pmfDist(t,pmfFun,n,tol)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 18-Oct-2018 15:37:55

pmf = zeros(n,1);
for i=1:n
    pmf(i) = tol*round(pmfFun{i}(t)/tol);
end

pmf = max(0,pmf);
pmf = pmf/sum(pmf);

end