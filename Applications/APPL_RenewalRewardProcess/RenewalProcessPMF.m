function [pmf,n,tMax,cdfFun,pmfFun,cfN,cfY] = ...
    RenewalProcessPMF(t,cfX,cfX1,nMax,options)
%RenewalProcessPMF
%  Evaluates the probability mass function (PMF) of the (delayed) renewal
%  process at time t.
%
%  The renewal process N(t) is a continuous-time process on the
%  non-negative integers n = 0,1,2,... which has independent inter-arrival
%  times (holding times) at each integer i before advancing to the next
%  integer, i+1. N(t) represents the number of arrivals to the system in
%  the interval (0,t]. The inter-arrival times may have any distribution on
%  the positive numbers, so long as they are independent and identically
%  distributed (IID) and have finite mean.
%
%  A delayed renewal process is just like an ordinary renewal process,
%  except that the first arrival time is allowed to have a different
%  distribution than the other interarrival times. Delayed renewal
%  processes arise naturally in applications and are also found embedded in
%  other random processes.
%
%  In RenewalProcessPMF, the delayed renewal process N(t) is specified by
%  the characteristic function cfX1 of the delay time (first arrival time
%  X1) and the characteristic function cfX of the other interarrival times
%  represented by the IID random variables X2,X3,... . For the ordinary
%  renewal processes we assume that cfX1 = cfX.
%
%  For each specified time t, RenewalProcessPMF evaluates the probability
%  mass function pmf(t) of N(t) defined on the non-negative integers n =
%  0,...,nMax. The pmf(t) can be evaluated at arbitrary time 0 < t <= tMax,
%  tMax is associated with the considered maximum number of arrivals nMax,
%  and is specified automatically by the algorithmPMF RenewalProcessPMF.
%
%  If time t is specified as empty set, t = [], then the RenewalProcessPMF
%  returns pmf as anonymous function of time t, i.e. pmf = @(t) pmf(t),
%  which can be evaluated for arbitrary time t from the interval (0,tMax].
%  Note that pmf(t=0) is concentrated with probability 1 at n = 0.
%
% SYNTAX:
%  pmf = RenewalProcessPMF(t,cfX)
%  pmf = RenewalProcessPMF(t,cfX,cfX1)
% [pmf,n,tMax,cdfFun,pmfFun,cfN,cfY] = ...
%   RenewalProcessPMF(t,cfX,cfX1,nMax,options)
%
% INPUTS:
%  t       - scalar value of time t >= 0, where the probability mass
%            function pmf(t) is evaluated. If empty, t = [], the
%            RenewalProcessPMF returns pmf as anonymous function of
%            time t, pmf = @(t) pmf(t);
%  cfX     - function handle of the characteristic function (CF) of the
%            holding time distribution (a continuous non-negative
%            distribution). If empty, cfX = [], the default value is the
%            characteristic function of the standard exponential
%            distribution, i.e. cfX = @(t)cf_Exponential(t);
%  cfX1    - function handle of the characteristic function (CF) of the
%            non-negative delay distribution (initial holding time). If
%            empty, cfX1 = [], the renewal process is a standard renewal
%            process (without delay) and its holding time specified by the
%            characteristic function cf;
%  nMax    - maximum value of the non-negative integers n = 0,...,nMax,
%            where the probability mass function pmf(t) is calculated. If
%            empty, nMax = [], the default value is nMax = 25;
% options  - structure for setting optional parameters for the used
%            numerical inversion algorithmPMF. For more details see cf2DistGP.
%            Moreover, options structure is used to set the following
%            parameters:
%            - tol : tolerance factor for the numerical errors of the
%            calculated probabilities, i.e. such that for Prob(n=k) < tol we
%            set Prob(n=k) = 0; If empty, default value is options.tol =
%            1e-6.
%            - cfW : function handle of the characteristic function (CF) of
%            the reward distribution. If options.cfW is specified, it is
%            used to construct the characteristic function cfY of the
%            renewal-reawrd process Y(t) = sum_{i=1}^N(t) W_i.
%            - algorithmPMF - selected algorithm for numerical inversion of
%            the characteristic function used computing PMF. Currently
%            available inversion algorithmPMFs include:
%             'cf2DistBTAV'  (Bromwich-Talbot-Abate-Valko method),
%             'cf2DistBV'    (Bakhvalov-Vasilieva method)
%             'cf2DistFFT'   (Fast Fourier Transform algorithmPMF method)
%             'cf2DistGP'    (Gil-Pelaez with Trapezoidal quadrature)
%             'cf2DistGPT'   (Gil-Pelaez with Trapezoidal quadrature)
%             If empty, default value is algorithmPMF = 'cf2DistGP'.
%
% OUTPUTS:
%  pmf    - vector of calculated probabilities evaluated for each integer n
%           = 0,...,nMax, at specified time t, i.e. the probability mass
%           function evaluated at specified time t, i.e. pmf = pmf(t).
%           Alternatively, if t = [], the result is an anonymous function
%           of time, i.e. pmf = @(t) pmf(t);
%  n      - vector of integers where the pmf(t) is evaluated, n = (0:nMax)
%  tMax   - the suggested upper bound of time t, where pmf(t) can be safely
%           evaluated;
%  cdfFun - auxiliary cell array of function handles of the n-times
%           convolved CDFs;
%  pmfFun - auxiliary cell array of function handles of the probabilities
%           for each specified integer n, expressed as a function of time
%           t.
%  cfN    - function handle of the characteristic function of the
%           distribution of the renewal process N(t) at time t, cfN =
%           @(t) cfE_DiracMixture(u,n,pmf). If t = [], then cfN = @(t,u)
%           cfE_DiracMixture(u,n,pmf(t));
%  cfY    - function handle of the characteristic function of the
%           distribution of the renewal-reward process Y(t) =
%           sum_{i=1}^N(t) W_i, at time t, cfY = @(u)
%           cfE_DiracMixture(u,n,pmf(t),cfW); If t = [], then cfY = @(t,u)
%           cfE_DiracMixture(u,n,pmf(t),cfW).
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
% % The renewal process with the exponentially distributed holding times
% % with the rate parameter lambda = 3.
% % This is equivalent with Poisson process with the parameter lambda = 3.
%  lambda  = 3;
%  cfX     = @(u)cf_Exponential(u,lambda);
%  cfX1    = [];
%  nMax    = 70;
%  [pmf,n] = RenewalProcessPMF([],cfX,cfX1,nMax);
%  plot(n,pmf(1),'o--',n,pmf(5),'o--',n,pmf(10),'o--',n,pmf(15),'o--')
%  title('Renewal process with exponential holding times t=[1,5,10,15]')
%  xlabel('n')
%  ylabel('probability')
%
% EXAMPLE 2:
% % The delayed renewal process with the exponentially distributed holding
% % times with the rate parameter lambda1 = 3 and with the exponentially
% % distributed delay time with the rate parameter lambda2 = 0.5.
%  lambda1 = 3;
%  lambda2 = 0.5;
%  cfX     = @(u)cf_Exponential(u,lambda1);
%  cfX1    = @(u)cf_Exponential(u,lambda2);
%  nMax    = 70;
%  [pmf,n] = RenewalProcessPMF([],cfX,cfX1,nMax);
%  plot(n,pmf(1),'o--',n,pmf(5),'o--',n,pmf(10),'o--',n,pmf(15),'o--')
%  title('Delayed renewal process with exponential holding times t=[1,5,10,15]')
%  xlabel('n')
%  ylabel('probability')
%
% EXAMPLE 3:
% % The renewal process with the gamma distributed holding times
% % with the shape parameter alpha = 2 and rate parameter beta = 5
%  alpha   = 2;
%  beta    = 5;
%  cfX     = @(u)cf_Gamma(u,alpha,beta);
%  cfX1    = @(u)cf_Exponential(u);
%  nMax    = 70;
%  ltxt    = cell(1,8);
%  clear options
%  options.algorithmPMF = 'cf2DistFFT';
%  [pmf,n,tMax] = RenewalProcessPMF([],cfX,cfX1,nMax,options);
%  figure;hold on
%  for t = 0:floor(tMax)
%      plot(n,pmf(t),'o--');
%  end
%  legend;
%  title('Delayed renewal process PMF(t) with Gamma(2,5) holding times t = 0:1:18')
%  xlabel('n')
%  ylabel('probability')

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 09-Sep-2020 14:05:07
%
% Revision history:
% Ver.: 04-Sep-2020 10:29:44
% Ver.: 18-Oct-2018 19:49:41

%% algorithmPMF
%  [pmf,n,tMax,cdfFun,pmfFun] = RenewalProcessPMF(t,cfX,cfX1,nMax,options)

%% CHECK THE INPUT PARAMETERS
narginchk(0, 5);
if nargin < 1, t         = []; end
if nargin < 2, cfX       = []; end
if nargin < 3, cfX1      = []; end
if nargin < 4, nMax      = []; end
if nargin < 5, options   = []; end

if ~isfield(options,'SixSigmaRule')
    options.SixSigmaRule = 6;
end

if ~isfield(options,'cfW')
    options.cfW = [];
end

if ~isfield(options,'isInterp')
    options.isInterp = true;
end

if ~isfield(options,'isPlot')
    options.isPlot = false;
end

if ~isfield(options,'xMin')
    options.xMin = 0;
end

if ~isfield(options,'tol')
    options.tol = 1e-6;
end

if ~isfield(options,'tailprob')
    options.tailprob = 1e-4;
end

if ~isfield(options,'maxiter')
    options.maxiter = 10;
end

if ~isfield(options,'algorithmPMF')
    options.algorithmPMF = [];
end

algorithmPMF = options.algorithmPMF;

if isempty(algorithmPMF)
    algorithmPMF = 'cf2DistGP';
end

if strcmp(algorithmPMF,'cf2DistFFT')
    if ~isfield(options,'N')
        options.N = 2^8;
    end
end

if strcmp(algorithmPMF,'cf2DistBTAV')
    if ~isfield(options,'trapezoidal')
        options.trapezoidal = 'trapezoidal';
    end
end

if strcmp(algorithmPMF,'cf2DistBV')
    if ~isfield(options, 'Limits')
        options.Limits = [1e-300 10.^linspace(-15,15,11) 1e+300];
    end
end

if isempty(nMax)
    nMax = 25;
end

if isempty(cfX)
    cfX = @(t) cf_Exponential(t);
end

if isempty(cfX1)
    cfX1 = @(t) cfX(t);
end

n      = (0:nMax)';
cdfFun = cell(nMax+2,1);
pmfFun = cell(nMax+1,1);

tMax = 0;
if nargout > 2
    tailprob = options.tailprob;
else
    tailprob = [];
end

cdfFun(1)      = {@(xnew) 1};
cdfFun(nMax+2) = {@(xnew) 0};
id = 1;
for k = 1:nMax
    id = id + 1;
    if k == nMax
        result = cf2Dist(@(t)cfX1(t).*cfX(t).^(k-1),[],tailprob, ...
            options,algorithmPMF);
        cdfFun(id) = {result.CDF};
        tMaxNew = result.qf;
        tMax = max(tMax,tMaxNew);
    else
        result = cf2Dist(@(t)cfX1(t).*cfX(t).^(k-1),[],[], ...
            options,algorithmPMF);
        cdfFun(id) = {result.CDF};
    end
end

if tMax == 0
    tMax = (result.xMean-result.xMin)/10;
end

id = 1;
pmfFun(id) = {@(t)cdfFun{id}(t)-cdfFun{id+1}(t)};
for k = 1:nMax
    id = id + 1;
    pmfFun(id) = {@(t)cdfFun{id}(t)-cdfFun{id+1}(t)};
end

if isscalar(t)
    pmf = pmfDist(t,pmfFun,nMax+1,options.tol);
    cfN = @(u) cfE_DiracMixture(u,n,pmf);
    cfY = @(u) cfE_DiracMixture(u,n,pmf,options.cfW);
else
    pmf = @(t) pmfDist(t,pmfFun,nMax+1,options.tol);
    cfN = @(t,u) cfE_DiracMixture(u,n,pmf(t));
    cfY = @(t,u) cfE_DiracMixture(u,n,pmf(t),options.cfW);
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
% Ver.: 29-Aug-2020 12:52:57

pmf = zeros(n,1);
for i=1:n
    pmf(i) = tol*round(pmfFun{i}(t)/tol);
end

pmf = max(0,pmf);
pmfsum = sum(pmf);
pmf = pmf/pmfsum;

end