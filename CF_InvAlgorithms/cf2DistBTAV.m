function [result,cdf,pdf,qf] = cf2DistBTAV(cf,x,prob,options)
%cf2DistBTAV Calculates the CDF/PDF/QF from the NON-NEGATIVE DISTRIBUTION
%  specified by the given characteristic function CF by using the
%  Bromwich-Talbot-Abate-Valko (BTAV) inversion method (originally
%  suggested as numerical inversion for the Laplace transform function).
%  Here we assume that the specified CF is a characteristic function of
%  nonnegative distribution which is well defined for complex valued
%  arguments.
%
% SYNTAX:
%  result = cf2DistBTAV(cf,x)
%  [result,cdf,pdf,qf] = cf2DistBTAV(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), 
%  x       - vector of x values where the CDF/PDF is computed, 
%  prob    - vector of values from [0,1] for which the quantiles
%            function is evaluated,
% options  - structure with the following parameters:
%             options.quadrature       % quadrature method. Default method
%                  = 'trapezoidal'     % is quadrature = 'trapezoidal'.
%                                      % Alternatively use quadrature =
%                                      % 'matlab' (for the MATLAB built in
%                                      % adaptive Gauss-Kronrod quadrature) 
%                                      % Gauss-Kronrod integral. 
%             options.tol = 1e-12      % absolute tolerance for the MATLAB
%                                      % Gauss-Kronrod integral. 
%             options.crit = 1e-12     % value of the criterion limit for
%                                      % stopping rule.
%             options.nTerms = 50      % number of terms used in the
%                                      % trapezoidal quadrature
%             options.maxiter = 100    % indicator of the maximum number of
%                                      % Newton-Raphson iterations.
%             options.qf0              % starting values of the quantiles.
%                                      % By default, the algorithm starts
%                                      % from the mean of the distribution
%                                      % estimated from the specified CF
%             options.Mpar_BTAV = 10   % parameter M for the deformed
%                                      % Bromwich curve.
%             options.isPlot = true    % plot the graphs of PDF/CDF
%
% OUTPUT:
%  result   - structure with CDF/PDF/QF and further details,
%  cdf      - vector of CDF values evaluated at x,
%  pdf      - vector of PDF values evaluated at x,
%  qf       - vector of QF values evaluated at prob.
%
% EXAMPLE 1
% % CDF/PDF/QF of the Chi-squared distribution with DF = 1
%  df = 1;
%  cf = @(t) (1 - 2i*t).^(-df/2);
%  prob = [0.9 0.95 0.99];
%  result = cf2DistBTAV(cf,[],prob)
%
% EXAMPLE 2 
% % CDF/PDF/QF of the linear combination (convolution) of Chi-squared RVs
%  df   = [1 2 3];
%  cf = @(t) cf_ChiSquare(t,df) ;
%  prob = [0.9 0.95 0.99];
%  result = cf2DistBTAV(cf,[],prob)
%
% EXAMPLE 3
% % CDF/PDF/QF of Bartlett null distribution 
%  k    = 5;
%  df   = [1 2 3 4 5 6 7 8 9 10];
%  cf = @(t) cfTest_Bartlett(t,df);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.quadrature = 'trapezoidal';
%  result = cf2DistBTAV(cf,[],prob,options)
%
% EXAMPLE 4
% % CDF/PDF/QF of the exact null distribution / Sphericity of CovMatrix
%  n  = 10;     % sample size for each population
%  p  = 5;      % dimension of each of the q populations
%  type = 'modified';
%  cf = @(t) cfTest_Sphericity(t,n,p,type);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.quadrature = 'trapezoidal';
%  result = cf2DistBTAV(cf,[],prob,options)
%
% EXAMPLE 5
% % CDF/PDF/QF of the exact null distribution / Compound Symmetry CovMatrix
%  n  = 10;     % sample size 
%  p  = 5;      % dimension of the populations
%  type = 'modified';
%  cf = @(t) cfTest_CompoundSymmetry(t,n,p,type);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.quadrature = 'trapezoidal';
%  result = cf2DistBTAV(cf,[],prob,options)
%
% EXAMPLE 6
% % CDF/PDF/QF of the exact null distribution / Equality Covariance Matrices
%  n  = 10;     % sample size for each population
%  p  = 5;      % dimension of each of the q populations
%  q  = 3;      % number of populations
%  type = 'modified';
%  cf = @(t) cfTest_EqualityCovariances(t,n,p,q,type);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.quadrature = 'trapezoidal';
%  result = cf2DistBTAV(cf,[],prob,options)
%
% EXAMPLE 7
% % CDF/PDF/QF of the exact null distribution / Equality of Means
%  n  = 10;     % sample size for each population
%  p  = 5;      % dimension of each of the q populations
%  q  = 3;      % number of populations
%  type = 'modified';
%  cf = @(t) cfTest_EqualityMeans(t,n,p,q,type);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.quadrature = 'trapezoidal';
%  result = cf2DistBTAV(cf,[],prob,options)
%
% EXAMPLE 8
% % CDF/PDF/QF of the exact null distribution / Equality of Populations
%  n  = 10;     % sample size for each population
%  p  = 5;      % dimension of each of the q populations
%  q  = 3;      % number of populations
%  type = 'modified';
%  cf = @(t) cfTest_EqualityPopulations(t,n,p,q,type);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.quadrature = 'trapezoidal';
%  result = cf2DistBTAV(cf,[],prob,options)
%
% EXAMPLE 9
% % CDF/PDF/QF of the exact null distribution / Test of Independence
%  n  = 20;     % sample size of the compound vector X = [X_1,...,X_q]
%  p  = 5;      % dimension of each of the q populations
%  q  = 3;      % number of populations
%  type = 'modified';
%  cf = @(t) cfTest_Independence(t,n,p,q,type);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.quadrature = 'trapezoidal';
%  result = cf2DistBTAV(cf,[],prob,options)
%
% NOTE OF CAUTION
%  The method was suggested for inverting proper Laplace tranform
%  functions. Here we use available characteristic functions for
%  creatig Laplace transform function, by using M(s) = cf(1i*s). In
%  general, the implemented algorithms for computing CFs assume that the
%  argument t is real. In specific situations CF is well defined also for
%  complex arguments. However, numerical issues could appear in any step
%  during the calculations. The result and the inner calculations should be
%  chcecked and properly controlled. 
%
% REFERENCES:
% [1] Talbot, A., 1979. The accurate numerical inversion of Laplace
%     transforms. IMA Journal of Applied Mathematics, 23(1), pp.97-120.
% [2] Abate, J. and Valkó, P.P., 2004. Multi-precision Laplace transform
%     inversion. International Journal for Numerical Methods in
%     Engineering, 60(5), pp.979-993.
%
% SEE ALSO: cf2Dist, cf2DistGPA, cf2DistGPT, cf2DistBTAV, cf2DistFFT,
%           cf2DistBV, cf2CDF_BTAV, cf2PDF_BTAV, cf2QF_BTAV

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 13:25:21
%
% Revision history:
% Ver.: 25-Dec-2018 13:59:10

%% ALGORITHM
%[result,cdf,pdf,qf] = cf2DistBTAV(cf,x,prob,options);

%% CHECK THE INPUT PARAMETERS
timeVal = tic;
narginchk(1, 4);

if nargin < 4, options = []; end
if nargin < 3, prob = []; end
if nargin < 2, x = []; end

if ~isfield(options, 'quadrature')
    options.quadrature = 'trapezoidal';
end

if ~isfield(options, 'tol')
    options.tol =1e-12;
end

if ~isfield(options, 'crit')
    options.crit = 1e-12;
end

if ~isfield(options, 'nTerms')
    options.nTerms = 50;
end

if ~isfield(options, 'maxiter')
    options.maxiter = 100;
end

if ~isfield(options, 'qf0')
    options.qf0 = (cf(1e-4)-cf(-1e-4))/(2e-4*1i);
end

if ~isfield(options, 'isCompound')
    options.isCompound = false;
end

if ~isfield(options, 'xMin')
        options.xMin = 0;
end

if ~isfield(options, 'xMax')
    options.xMax = Inf;
end

if ~isfield(options, 'xMean')
    options.xMean = [];
end

if ~isfield(options, 'xStd')
    options.xStd = [];
end

if ~isfield(options, 'SixSigmaRule')
    if options.isCompound
        options.SixSigmaRule = 10;
    else
        options.SixSigmaRule = 6;
    end
end

if ~isfield(options, 'tolDiff')
    options.tolDiff = 1e-4;
end

if ~isfield(options, 'crit')
    options.crit = 1e-10;
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

if ~isfield(options, 'xN')
    options.xN = 101;
end

if ~isfield(options, 'chebyPts')
    options.chebyPts = 2^9;
end

if ~isfield(options, 'isInterp')
    options.isInterp = false;
end

if ~isfield(options, 'tol')
    options.tol = 1e-10;
end

%% GET/SET the DEFAULT parameters and the OPTIONS
xMin    = options.xMin;
xMax    = options.xMax;
xMean   = options.xMean;
xStd    = options.xStd;
tolDiff = options.tolDiff;
cft     = cf(tolDiff*(1:4));
cftRe   = real(cft);
cftIm   = imag(cft);
SixSigmaRule = options.SixSigmaRule;
    
if isempty(xMean)
    xMean = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 ...
             - 2*cftIm(4)/280) / tolDiff;
end

if isfinite(xMean)
    options.xMean = xMean;
end

if isempty(xStd)
    xM2   = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 ...
        - 16*cftRe(3)/315 + 2*cftRe(4)/560) / tolDiff^2;
    xStd  = sqrt(xM2 - xMean^2);
end


if ~isfield(options, 'xN')
    options.xN = 101;
end

isPlot = options.isPlot;
options.isPlot = false;

%% ALGORITHM

if isempty(x)
    xempty = true;
    xMin = max(xMin,xMean - SixSigmaRule * xStd);
    xMax = min(xMax,xMean + SixSigmaRule * xStd);
    x = linspace(xMax,xMin,options.xN);
else
    xempty = false;
end

if options.isInterp
    x0 = x;
    % Chebyshev points
    x = (xMax-xMin) * (-cos(pi*(0:options.chebyPts) / ...
        options.chebyPts) + 1) / 2 + xMin;
else
    x0 = [];
end

cdf = cf2CDF_BTAV(cf,x,options);
pdf = cf2PDF_BTAV(cf,x,options);

if ~isempty(prob)
    qf = cf2QF_BTAV(cf,prob,options);
else
    qf = [];
end
options.isPlot = isPlot;

%% Create the INTERPOLAN functions
if options.isInterp
    id   = isfinite(pdf);
    PDF  = @(xnew)InterpPDF(xnew,x(id),pdf(id));
    id   = isfinite(cdf);
    CDF  = @(xnew)InterpCDF(xnew,x(id),cdf(id));
    QF   = @(prob)InterpQF(prob,x(id),cdf(id));
    RND  = @(dim)InterpRND(dim,x(id),cdf(id));
    try
    if ~xempty
        x   = x0;
        cdf = CDF(x);
        pdf = PDF(x);
    end
    catch
        warning('VW:CharFunTool:cf2DistBTAV', ...
            'Problem using the interpolant function');
    end
else
    PDF  = [];
    CDF  = [];
    QF   = [];
    RND  = [];
end

%% RESULT
result.Description         = 'CDF/PDF/QF from the characteristic function CF';
result.inversionMethod     = 'Bromwich-Talbot-Abate-Valkó';
result.quadratureMethod    = options.quadrature;
result.x                   = x;
result.cdf                 = cdf;
result.pdf                 = pdf;
result.prob                = prob;
result.qf                  = qf;
if options.isInterp
    result.PDF             = PDF;
    result.CDF             = CDF;
    result.QF              = QF;
    result.RND             = RND;
end
result.cf                  = cf;
result.isInterp            = options.isInterp;
result.SixSigmaRule        = options.SixSigmaRule;
result.xMean               = xMean;
result.xStd                = xStd;
result.xMin                = xMin;
result.xMax                = xMax;
result.options             = options;
result.tictoc              = toc(timeVal);

%% PLOT the PDF / CDF 
if length(x)==1
    isPlot = false;
end
if isPlot
    % PDF
    plot(x,pdf,'.-')
    grid
    title('PDF Specified by the CF')
    xlabel('x')
    ylabel('pdf')    
    % CDF
    figure
    plot(x,cdf,'.-')
    grid
    title('CDF Specified by the CF')
    xlabel('x')
    ylabel('cdf')   
end
end