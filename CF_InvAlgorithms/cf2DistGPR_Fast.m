function [result,pdf,cdf,qf,x,prob] = cf2DistGPR_Fast(cf,x,prob,options)
% cf2DistGPR_Fast is a FAST and highly SIMPLIFIED version of the GPR
%  algorithm (Gil-Pelaez inversion formulae with Riemann sum used for
%  integration) to evaluate the APPROXIMATE values of PDF, CDF and the QF
%  of the continuous univariate distribution specified by its
%  characteristic function CF based on using barycentric INTERPOLATION of
%  the fitted CDF/PDF.
%  
% SYNTAX:
%   result = cf2DistGPR_Fast(cf,x,prob,options)
%   [result,pdf,cdf,qf,x,prob] = cf2DistGPR_Fast(cf,x,prob,options)
%   [~,pdf,cdf,qf,x,prob] = cf2DistGPR_Fast(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), 
%  x       - vector of x values where the CDF/PDF is evaluated, 
%  prob    - vector of values from [0,1] for which the quantiles
%            function is evaluated,
%  options - structure with the following default parameters:
%  xMin     - the lower limit of X
%  xMax     - the lower limit of X
%  N        - number of points used for trapezoidal integration. If empty,
%             by default N = 2^10. 
%  chebyPts - number of Chebyshev points to create the vector of x values,
%             x = [xMin, xMax], where the CDF of the distribution specified
%             by its CF is evaluated (xMin and xMax are specified
%             automatically from the CF). If empty, by default chebyPts =
%             2^7.  
%  SixSigmaRule - multiple of standard deviation which specifies the rule
%             for estimating the principal domain of the distribution. If
%             empty, by default SixSigmaRule = 6. 
%  tolDiff  - parameter used for numerical differentiation. If empty, by
%             default tolDiff = 1e-4.  
%
% OUTPUT:
%  result  - structure with CDF/PDF/QF and further details,
%  cdf     - vector of CDF values evaluated at x,
%  pdf     - vector of PDF values evaluated at x,
%  qf      - vector of QF values evaluated at prob,
%  x       - vector of x values where the CDF/PDF were evaluated, 
%  prob    - vector of values from [0,1] for which the quantiles were
%            evaluated. 
%
% REMARK:
%  cf2DistGPR_Fast was suggested as a fast version of the other inversion
%  algorithm, in particular cf2DistGPR, in order to calculate the
%  approximately the required PDF/CDF and QF values. This can be used for
%  fast random sampling from the distribution specified its characteristic
%  function CF.
%
% EXAMPLE (Calculate PDF/CDF/QF from distribution specified by CF)
%  cf = @(t) exp(-t.^2/2);
%  result = cf2DistGPR_Fast(cf)
%
% EXAMPLE (Calculate specified quantiles from the CF)
%  cf = @(t) cf_Normal(t);
%  prob = (0:0.01:1)'; 
%  [~,~,~,qf,~,prob] = cf2DistGPR_Fast(cf,[],prob)
%
% EXAMPLE (Calculate PDF/CDF/QF from distribution specified by CF)
%  cf = @(t) cf_Chi(t,5);
%  x = linspace(0,5,51)';
%  prob = (0:0.05:1)'; 
%  clear options
%  options.xMin = 0;
%  options.SixSigmaRule = 6;
%  options.chebyPts = 201;
%  [result,pdf,cdf,qf,x,prob] = cf2DistGPR_Fast(cf,x,prob,options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 05-Jan-2022 17:03:20

%% CHECK THE INPUT PARAMETERS
timeVal = tic;
narginchk(1, 4);

if nargin < 4, options = []; end
if nargin < 3, prob = []; end
if nargin < 2, x = []; end

if ~isfield(options, 'xMin')
    options.xMin = -Inf;
end

if ~isfield(options, 'xMax')
    options.xMax = Inf;
end

if ~isfield(options, 'N')
        options.N = 2^10;
end

if ~isfield(options, 'SixSigmaRule')
        options.SixSigmaRule = 6;
end

if ~isfield(options, 'tolDiff')
    options.tolDiff = 1e-4;
end

if ~isfield(options, 'xN')
    options.xN = 101;
end

if ~isfield(options, 'chebyPts')
    options.chebyPts = 2^7;
end

if isempty(prob)
    prob = linspace(0,1,options.xN)'; 
end

%% ALGORITHM
cft      = cf(options.tolDiff*(1:4));
cftRe    = real(cft);
cftIm    = imag(cft);

xMean    = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 - ...
            2*cftIm(4)/280) / options.tolDiff;
xM2      = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 - ...
            16*cftRe(3)/315 + 2*cftRe(4)/560) / options.tolDiff^2;
xStd     = sqrt(xM2 - xMean^2);

xMin     = max(options.xMin,xMean - options.SixSigmaRule * xStd);
xMax     = min(options.xMax,xMean + options.SixSigmaRule * xStd);

if isempty(x) 
    x = linspace(xMin,xMax,options.xN)'; 
end

range    = xMax - xMin;
dt       = 2*pi / range;
t        = (0.5+(0:options.N-1))' * dt;
cft      = cf(t);
cft(options.N) = cft(options.N)/2;

xCheby   = range * (-cos(pi*(0:(options.chebyPts-1)) / ...
    (options.chebyPts-1)) + 1) / 2 + xMin;
xCheby   = xCheby(:);

ExpCheby = exp(-1i*xCheby*t');
pdfCheby = max(0,(real(ExpCheby * cft) * dt) / pi);
cdfCheby = max(0,min(1,0.5 - (imag(ExpCheby * (cft ./ t)) * dt) / pi));
cdfMin   = min(cdfCheby);
cdfMax   = max(cdfCheby);

prob(prob<cdfMin) = NaN;
prob(prob>cdfMax) = NaN;
szp  = size(prob);
prob = prob(:);
[probSort,probID] = sort(prob);

x(x<xMin) = NaN;
x(x>xMax) = NaN;
szx = size(x);
x   = x(:);
[xSort,xID] = sort(x);

pdf(xID)   = max(0,InterpBarycentric(xCheby,pdfCheby,xSort));
cdf(xID)   = max(0,min(1,InterpBarycentric(xCheby,cdfCheby,xSort)));
qf(probID) = InterpBarycentric(cdfCheby,xCheby,probSort);

pdf  = reshape(pdf,szx);
cdf  = reshape(cdf,szx);
qf   = reshape(qf,szp);
x    = reshape(x,szx);
prob = reshape(prob,szp); 

%% RESULT
result.Description         = 'FAST GPR algorithm for CDF/PDF/QF from CF';
result.inversionMethod     = 'Gil-Pelaez';
result.quadratureMethod    = 'Riemann sum quadrature';
result.pdf                 = pdf;
result.cdf                 = cdf;
result.qf                  = qf;
result.x                   = x;
result.prob                = prob;
result.pdfCheby            = pdfCheby;
result.cdfCheby            = cdfCheby;
result.xCheby              = xCheby;
result.cf                  = cf;
result.SixSigmaRule        = options.SixSigmaRule;
result.N                   = options.N;
result.dt                  = dt;
result.T                   = t(end);
result.xMean               = xMean;
result.xStd                = xStd;
result.xMin                = xMin;
result.xMax                = xMax;
result.options             = options;
result.tictoc              = toc(timeVal);
end