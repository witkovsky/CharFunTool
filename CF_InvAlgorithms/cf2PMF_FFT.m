function [result,pmf,cdf] = cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
% cf2PMF_FFT Evaluates the distribution functions (PMF and CDF) of a
%  DISCRETE random variable (DRV) with lattice distribution defined on the
%  Support = [xMin:xDelta:xMax], and fully specified by its characteristic
%  function (CF). PMF and CDF is evaluated at the specified support points
%  x = [xMin:xDelta:xMax] by numerical inversion of CF, based on using the
%  inverse FFT algorithm.
%
%  We specifically focus on distributions with the support points specified
%  by n*xDelta for n = 0,±1,±2,... and xDelta > 0. This type of
%  distribution has been termed a lattice or an arithmetic distribution.
%  These include most of the well studied discrete distributions such as
%  the binomial, Poisson, negative binomial and many others, including the
%  empirical distributions defined on integers or any other lattices, and
%  also their well-defined combinations (e.g. convolutions).
%
% The cf2PMF_FFT algorithm considers only finite discrete support, so it
% should be specified carefully and correctly (i.e., it should include a
% probability mass equal to 1 - epsilon, where epsilon is an accepted
% truncation error). See Warr (2014) for more information on the method
% used.    
% 
% SYNTAX:
% [result,pmf,cdf] = cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
%
% INPUT:
%  cf       - function handle of the characteristic function of the
%             discrete lattice distribution, defined on a subset of
%             integers.
%  xMin     - minimum value of the support of the (lattice) discrete RV X.
%  xMax     - maximum value (finite number) of the support of the (lattice)
%             DRV X.  
%  xDelta   - xDelta > 0 is the minimum difference of the support values of
%             the discrete RV X with lattice distribution given as Supp =
%             xMin:xDelta:xMax. If empty, the default value is xDelta = 1.  
%  options  - structure with the following default parameters:
%             options.xMin = 0       % minimum value of the support of X
%             options.xMax = 100     % maximum value of the support of X
%             options.xDelta = 1     % minimum difference of the support
%                                    % values of X 
%             options.isPlot = true  % logical indicator for plotting 
%                                    % PMF and CDF
%
% OUTPUT:
%  result   - structure with all details.
%  pmf      - vector of the PMF values evaluated at x = xMin:xMax.
%  cdf      - vector of the CDF values evaluated at x = xMin:xMax.
%
% EXAMPLE 1
% % PMF/CDF of the binomial RV specified by its CF
%  n  =  10;
%  p  =  0.25;
%  cf = @(t) cfN_Binomial(t,n,p);
%  xMin   = 0;
%  xMax   = n;
%  xDelta = 1;
%  options.isPlot = true;
%  result = cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
%
% EXAMPLE 2
% % PMF/CDF of the convolved discrete RV specified by its CF
%  n  =  10;
%  p  =  0.25;
%  cf_Bino = @(t) cfN_Binomial(t,n,p);
%  N  = 5;
%  cf = @(t) cf_Bino(t).^N;
%  xMin   = 0;
%  xMax   = N*n;
%  xDelta = 1;
%  options.isPlot = true;
%  [result,pmf,cdf] = cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
%
% EXAMPLE 3
% % PMF/CDF of the mean of IID RV with Poisson distribution 
%  lambda = 5;
%  cf_Pois = @(t) cfN_Poisson(t,lambda);
%  N  = 10;
%  cf = @(t) cf_Pois(t/N).^N;
%  xMin   = 0;
%  xMax   = 10;
%  xDelta = 1/N;
%  options.isPlot = true;
%  result = cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
%
% EXAMPLE 4
% % PMF/CDF of the convolved discrete RV specified by its CF
% % Here we consider convolutions of discrete RV defined on {0,1,2}
% % with probabilities p = [0.2,0.5,0.3]
%  supp   = [0,1,2];
%  prob   = [0.2,0.5,0.3];
%  cf_X   = @(t) cfE_DiracMixture(t,supp,prob)
%  N      = 5;
%  cf     = @(t) cf_X(t).^N;
%  xMin   = 0;
%  xMax   = max(supp)*N;
%  xDelta = 1;
%  options.isPlot = true;
%  result = cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
%
% EXAMPLE 5
% % PMF/CDF of the exact bootrap mean distribution specified by its CF
%  data   = [-2,0,0,0,0,0,0,0,1,4];
%  N      = length(data);
%  cf     = @(t) cfE_Empirical(t/N,data).^N
%  xMin   = min(data);
%  xMax   = max(data);
%  xDelta = 1/N;
%  options.isPlot = true;
%  result = cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
%
% REFERENCES:
% [1] Warr, Richard L. Numerical approximation of probability mass
%     functions via the inverse discrete Fourier transform. Methodology
%     and Computing in Applied Probability 16, no. 4 (2014): 1025-1038.   

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 11-Oct-2019 23:38:25
% Revised: 28-Aug-2021 15:23:16

%% CHECK/SET THE INPUT PARAMETERS
StartTime = cputime;
narginchk(1, 5);
if nargin < 5, options = []; end
if nargin < 4, xDelta  = []; end
if nargin < 3, xMax    = []; end
if nargin < 2, xMin    = []; end

if ~isfield(options, 'xMin')
    options.xMin = 0;
end

if ~isfield(options, 'xMax')
    options.xMax = 100;
end

if ~isfield(options, 'xDelta')
    options.xDelta = 1;
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

if isempty(xMax)
    xMax  = options.xMax;
end

if isempty(xMin)
    xMin  = options.xMin;
end

if isempty(xDelta)
    xDelta = options.xDelta;
end

%% ALGORITHM
x     = xMin:xDelta:xMax;
N     = length(x);

if xMin == 0
    omega = (x/xDelta)/N;
    DFTfun = cf(-2*pi*omega/xDelta);
else
    xx    = x - xMin;
    omega = (xx/xDelta)/N;
    DFTfun = cf(-2*pi*omega/xDelta) .* exp(1i*2*pi*omega*xMin/xDelta);
end

% PMF
pmf   = real(ifft(DFTfun));
pmf   = max(0,pmf);
pmf(pmf < 100*eps) = 0;

% CDF
cdf    = cumsum(pmf);
tictoc = cputime - StartTime;

%% RESULTS
result.Description     = 'PMF/CDF of a discrete distribution from its CF';
result.inversionMethod = 'inverse FFT algorithm';
result.pmf = pmf;
result.cdf = cdf;
result.x   = x;
result.cf  = cf;
result.DFTfun  = DFTfun;
result.options = options;
result.tictoc  = tictoc;

%% PLOT the PDF / CDF 
if length(x)==1
    options.isPlot = false;
end
if options.isPlot
    % PMF
    figure
    stem(x,pmf)
    xlabel('x')
    ylabel('pmf')
    title('PMF Specified by the CF')
    % CDF
    figure
    stairs(x,cdf)
    grid
    title('CDF Specified by the CF')
    xlabel('x')
    ylabel('cdf')        
end
end