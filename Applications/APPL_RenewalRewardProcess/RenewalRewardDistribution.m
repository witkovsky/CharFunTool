function [pdf,cdf,qf,y,n,pmf,tMax,cfY,result] = ... 
    RenewalRewardDistribution(t,y,prob,cfW,cfX,cfX1,nMax,options)
%RenewalRewardDistribution 
%  Evaluates the distribution functions PDF/CDF/Qf, the probability density
%  function, the cumulative distribution function, and the quantile
%  function of the renewal-reward process at time t for specified values y
%  (PDF and CDF) and QF for specified probabilities prob.
%  
%  The renewal-reward process is defined as Y(t) = sum_{i=1}^ND(t) Wi,
%  where ND(t) is the deleyed renewal process and Wi are IID random
%  variables representing the rewards. Note that Wi may take negative
%  values as well as positive values. In fact, Y(t) is specified by the
%  characteristic function cfW of the independent rewards Wi, by the
%  characteristic function cfX1 of the delay time (first arrival time X1),
%  and by the characteristic function cfX of the other inter-arrival
%  (holding) times, represented by the IID random variables X2,X3,...,Xn.
%
%  For each specified time t the algorithm evaluates the PDF/CDF/QF of the
%  compound renewal-reward random variable Y(t). The distribution can be
%  evaluated at arbitrary time t >= 0. However, for correct result t must
%  be from the associated interval (tMin,tMax), which can be specified by
%  the algorithm. 
%
% SYNTAX:
%  pdf = RenewalRewardDistribution(t,y,[],cfW,cfX,cfX1,nMax)
%  [pdf,cdf,qf] = RenewalRewardDistribution(t,y,prob,cfW,cfX,cfX1,nMax)
%  [pdf,cdf,qf,y,n,pmf,tMax,cfY,result] = ...
%    RenewalRewardDistribution(t,y,prob,cfW,cfX,cfX1,nMax,options)
% 
% INPUTS:
%  t        - scalar value of time t >= 0, where the probability mass
%             function pmf(t) is evaluated. If empty, t = [], the
%             RenewalRewardDistribution returns pmf as anonymous function
%             of time t, pmf = @(t) pmf(t);
%  y        - vector of y values where the CDF/PDF is computed. If y = [],
%             the values of y are specified automatically;
%  prob     - vector of values from [0,1] for which the quantiles function
%             QF is evaluated;
%  cfW      - function handle of the characteristic function (CF) of the
%             distribution of the renewal-reawrd random variable Y(t) =
%             sum_{i=1}^N(t) W_i at time t. If empty, cfW = [], the default
%             value is the characteristic function of the standard normal
%             distribution, i.e. cfW = @(u) cf_Normal(u);
%  cfX      - function handle of the characteristic function (CF) of the
%             holding time distribution (a continuous non-negative
%             distribution). If empty, cfX = [], the default value is the
%             characteristic function of the standard exponential
%             distribution, i.e. cfX = @(u)cf_Exponential(u);
%  cfX1     - function handle of the characteristic function (CF) of the
%             non-negative delay distribution (initial holding time). If
%             empty, cfX1 = [], the renewal process is a standard renewal
%             process (without delay) and its holding time specified by the
%             characteristic function cf;
%  nMax     - maximum value of the integers n, where the probability mass
%             function pmf(t) is calculated. If empty, nMax = [], the
%             default value is nMax = 25;
% options   - structure for setting optional parameters for the used
%             numerical inversion algorithm. 
%             - algorithmPMF - selected inversion algorithm used for
%              computing PMF at specified time t. 
%             - algorithm - selected inversion algorithm used for computing
%             renewal-reward distribution Y(t) at specified time t. 
%             Currently available inversion algorithms include:
%             'cf2DistBTAV'  (Bromwich-Talbot-Abate-Valko method),
%             'cf2DistBV'    (Bakhvalov-Vasilieva method)
%             'cf2DistFFT'   (Fast Fourier Transform algorithm method)
%             'cf2DistGP'    (Gil-Pelaez with Trapezoidal quadrature)
%             'cf2DistGPT'   (Gil-Pelaez with Trapezoidal quadrature)
%             If empty, default value is algorithm = 'cf2DistGP'.
% 
% OUTPUTS:
%  pdf      - vector of PDF values of Y(t) evaluated at specified y,
%  cdf      - vector of CDF values of Y(t) evaluated at specified y,
%  qf       - vector of QF values evaluated at specified values prob.
%  y        - vector of y values where the CDF/PDF was computed;
%  n        - vector of integers where the pmf(t) was evaluated;
%  pmf      - vector of calculated probabilities evaluated for each integer
%             n in (nMin,nMax) at specified time t, i.e. the probability
%             mass function evaluated at specified time t, i.e. pmf =
%             pmf(t). Alternatively, if t = [], the result is an anonymous
%             function of time, i.e. pmf = @(t) pmf(t);
%  tMax     - the suggested upper bound of time t, where pmf(t) can be
%             safely evaluated;
%  cfY      - function handle of the characteristic function of the
%             distribution of the renewal-reward process Y(t) =
%             sum_{i=1}^N(t) W_i, at time t, cfY = @(u)
%             cfE_DiracMixture(u,n,pmf(t),cfW); If t = [], then cfY =
%             @(t,u) cfE_DiracMixture(u,n,pmf(t),cfW). 
%  result   - structure output of the specified inversion algorithm with
%             CDF/PDF/QF and further details.
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Renewal_theory
%
% EXAMPLE 1:
% % The renewal-reward process with normally distributed rewards, Wi ~
% % N(1,1), and renewal process defined by the exponentially distributed
% % holding Xi ~ Exp(3) for all i = 1,2,3,... . 
%  mu      = 1;
%  sigma   = 1;
%  cfW     = @(u)cf_Normal(u,mu,sigma);
%  lambda  = 3;
%  cfX     = @(u)cf_Exponential(u,lambda);
%  cfX1    = []; 
%  nMin    = 0;
%  nMax    = 70;
%  t       = 10;
%  y       = [];
%  prob    = [0.9 0.95 0.99];
%  [pdf,cdf,qf,y,n,pmf,tMax,cfY,result] = ... 
%    RenewalRewardDistribution(t,y,prob,cfW,cfX,cfX1,nMax);
% figure; stem(n,pmf(t));title('PMF of the renewal process N(t) at time t');
% figure; plot(y,pdf);grid on;title('PDF of the renewal-reward process Y(t) at time t');
% figure; plot(y,cdf);grid on;title('CDF of the renewal-reward process Y(t) at time t');

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 03-Sep-2020 15:23:32

%% ALGORITHM
%  [pdf,cdf,qf,y,n,pmf,tMax,cfY,result] = ... 
%   RenewalRewardDistribution(t,y,prob,cfW,cfX,cfX1,nMax,options)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 8);
if nargin < 2,  y            = []; end
if nargin < 3,  prob         = []; end
if nargin < 4,  cfW          = []; end
if nargin < 5,  cfX          = []; end
if nargin < 6,  cfX1         = []; end
if nargin < 7,  nMax         = []; end
if nargin < 8,  options      = []; end

if isempty(cfW)
    cfW = @(t) cf_Normal(t);
end

if isempty(cfX)
    cfX = @(t) cf_Exponential(t);
end

if isempty(nMax)
    nMax = 50;
end

if ~isfield(options,'isInterp')
    options.isInterp = true;
end

if ~isfield(options,'isPlot')
    options.isPlot = false;
end
 
if ~isfield(options,'algorithmPMF')
    options.algorithmPMF = [];
end

if ~isfield(options,'algorithm')
    options.algorithm = 'cf2DistGP';
end

% Evaluate the renewal process PMF at time t and the renewal-reward
% characteristic function cfY of the compound distribution of Y(t),
% specified by pmf(t)and cfW, at time t 
% cfY = @(u) cfE_DiracMixture(u,n,pmf(t),cfW);
options.cfW = cfW;
[pmf,n,tMax,~,~,~,cfY] = RenewalProcessPMF([],cfX,cfX1,nMax,options);

options.isInterp = false;
result = cf2Dist(@(u)cfY(t,u),y,prob,options,options.algorithm);
pdf    = result.pdf;
cdf    = result.cdf;
qf     = result.qf;
y      = result.x;

end