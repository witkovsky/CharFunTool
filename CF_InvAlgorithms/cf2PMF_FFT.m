function [pmf,cdf,result] = cf2PMF_FFT(cf,xMin,xMax,options)
% cf2PMF_FFT Evaluates the distribution functions, PMF and CDF, of a
% DISCRETE distribution defined on a (subset) of non-negative integers
% specified by its CF, at values x = xMin : xMax, by using numerical
% inversion of the characteristic function, based on the inverse FFT
% algorithm.  
% 
% SYNTAX:
% [pmf,cdf,result] = cf2PMF_FFT(cf,xMin,xMax,options)
%
% INPUT:
%  cf       - function handle of the characteristic function
%  xMin     - minimum value (integer) of the support of the discrete RV X
%  xMax     - maximum value (integer) of the support of the discrete RV X
%  options  - structure with the following parameters:
%             options.xMin = 0       % minimum value of the support of X
%             options.xMax = 100     % maximum value of the support of X
%             options.isPlot = false % logical indicator for plotting 
%                                    % PMF and CDF
%
% OUTPUT:
%  pmf      - vector of the PMF values evaluated at x = xMin:xMax.
%  cdf      - vector of the PMF values evaluated at x = xMin:xMax.
%  result   - structure with further details.
%
% EXAMPLE 1
% % PMF/CDF of the binomial RV specified by its CF
%  n  =  10;
%  p  =  0.25;
%  cf = @(t) cfN_Binomial(t,n,p);
%  xMin = 0;
%  xMax = n;
%  options.isPlot = true;
%  [pmf,cdf,result] = cf2PMF_FFT(cf,xMin,xMax,options)
%
% EXAMPLE 2
% % PMF/CDF of the convolved discrete RV specified by its CF
%  n  =  10;
%  p  =  0.25;
%  cf_Bino = @(t) cfN_Binomial(t,n,p);
%  N  = 5;
%  cf = @(t) cf_Bino(t).^N;
%  xMin = 0;
%  xMax = N*n;
%  options.isPlot = true;
%  [pmf,cdf,result] = cf2PMF_FFT(cf,xMin,xMax,options)
%
% EXAMPLE 3
% % PMF/CDF of the convolved RV with Poisson distribution 
%  lambda = 5;
%  cf_Pois = @(t) cfN_Poisson(t,lambda);
%  N  = 10;
%  cf = @(t) cf_Pois(t).^N;
%  xMin = 0;
%  xMax = 3*lambda*N;
%  options.isPlot = true;
%  [pmf,cdf,result] = cf2PMF_FFT(cf,xMin,xMax,options)
%
% EXAMPLE 4
% % PMF/CDF of the convolved discrete RV specified by its CF
% % Here we consider convolutions of discrete RV defined on {0,1,2}
% % with probabilities p = [0.2,0.5,0.3]
%  supp  =  [0,1,2];
%  prob  =  [0.2,0.5,0.3];
%  cf_X  = @(t) cfE_DiracMixture(t,supp,prob)
%  N     = 5;
%  cf    = @(t) cf_X(t).^N;
%  xMin  = 0;
%  xMax  = N*2;
%  options.isPlot = true;
%  [pmf,cdf,result] = cf2PMF_FFT(cf,xMin,xMax,options)
%
% REFERENCES:
% [1] Warr, Richard L. Numerical approximation of probability mass
%     functions via the inverse discrete Fourier transform. Methodology
%     and Computing in Applied Probability 16, no. 4 (2014): 1025-1038.   

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 22-Feb-2019 09:23:36
%% CHECK/SET THE INPUT PARAMETERS
StartTime = cputime;
narginchk(1, 4);
if nargin < 4, options = []; end
if nargin < 3, xMax = []; end
if nargin < 2, xMin = []; end

if ~isfield(options, 'xMin')
    options.xMin = 0;
end

if ~isfield(options, 'xMax')
    options.xMax = 100;
end

if ~isfield(options, 'isPlot')
    options.isPlot = false;
end

if isempty(xMax)
    xMax  = options.xMax;
end

if isempty(xMin)
    xMin  = options.xMin;
end

%% ALGORITHM
range = xMax-xMin;
N     = range +1;
x     = xMin:xMax;
omega = x/N;
ftFun = cf(-2*pi*omega);
% PMF
pmf   = real(ifft(ftFun));
pmf   = max(0,pmf);
pmf(pmf < 100*eps) = 0;
% CDF
cdf   = cumsum(pmf);
tictoc = cputime - StartTime;

%% RESULTS
if nargout > 1
    result.Description     = 'PMF/CDF of a discrete distribution from its CF';
    result.inversionMethod = 'inverse FFT algorithm';
    result.pmf = pmf;
    result.cdf = cdf;
    result.x   = x;
    result.cf  = cf;
    result.ftFun   = ftFun;
    result.options = options;
    result.tictoc  = tictoc; 
end

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