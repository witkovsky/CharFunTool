function [result,cdf,pdf,qf] = cf2DistGPA(cf,x,prob,options)
%cf2DistGPA Calculates the CDF/PDF/QF from the characteristic function CF
% by using the Gil-Pelaez inversion formulae:
%   cdf(x) = 1/2 + (1/pi) * Integral_0^inf imag(exp(-1i*t*x)*cf(t)/t)*dt.
%   pdf(x) = (1/pi) * Integral_0^inf real(exp(-1i*t*x)*cf(t))*dt.
%
%  The required FOURIER INTEGRALs are calculated by using the adaptive
%  Gauss-Kronrod quadrature rule for numerical integration of the
%  oscillatory integrand function divided into sub-intervals (found by a
%  fast root-finding algorithm) and subsequent application of the
%  convergence acceleration techniques for computing the limit of the
%  resulted alternating series, for more details see, e.g., Cohen et al.
%  (2000), Sidi (2011).
%
% SYNTAX:
%  result = cf2DistGPA(cf,x)
%  [result,cdf,pdf,qf] = cf2DistGPA(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), 
%  x       - vector of x values where the CDF/PDF is computed, 
%  prob    - vector of values from [0,1] for which the quantiles
%            function is evaluated,
%  options - structure with the following default parameters:
%             options.isCompound = false % indicator of the compound
%                                        % distributions
%             options.isCircular = false % indicator of the circular
%                                        % distributions on (-pi,pi)
%             options.isInterp   = false % create and use the interpolant
%                                          functions for PDF/CDF/QF/RND
%             options.xMin = -Inf      % set the lower limit of X
%             options.xMax = Inf       % set the lower limit of X
%             options.xMean = []       % set the MEAN value of X
%             options.xStd = []        % set the STD value of X
%             options.SixSigmaRule = 6 % set the rule for automatic
%                                      % computation of the domain of X
%             options.tolDiff = 1e-4   % tolerance for numerical
%                                      % differentiation 
%             options.tol = 1e-10      % tolerance for numerical
%                                      % integration
%             options.tolFindRoots     % tolerance for root-finding
%                                      % procedure. Default value is
%                                      % options.tolFindRoots = 1e-32.
%                                      % Alternatively, set lower levels up
%                                      % to options.tolFindRoots = 1e-321.
%             options.maxiterFindRoots % maximum number of iterations the
%                                      % root-finding procedure. Default
%                                      % value is options.maxiterFindRoots
%                                      % = 100. 
%             options.isPlot = true    % plot the graphs of PDF/CDF
%             options.isAccelerated = true % indicator of activated
%                                      % acceleration algorithm
%             options.nPeriods = 25    % the upper integration limit: UPPER
%                                      % = nPeriods * pi / x. The the basic
%                                      % integration interval [0,UPPER] is
%                                      % devided into two subintervals 
%                                      % [0 A] and [A UPPER]. 
%             options.nPoly            % order of the chebyshev polynomials
%                                      % used for the root-finding
%                                      % procedure. Default value is
%                                      % options.nPoly = 2^5.
%
% OUTPUT:
%  result   - structure with CDF/PDF/QF and further details,
%  cdf      - vector of CDF values evaluated at x,
%  pdf      - vector of PDF values evaluated at x,
%  qf       - vector of QF values evaluated at prob.
%
% EXAMPLE1 (Calculate CDF/PDF of N(0,1) by inverting its CF)
%  cf = @(t) exp(-t.^2/2);
%  result = cf2DistGPA(cf)
%
% EXAMPLE2 (Calculate CDF/PDF of a linear combination of chi-squared RVs)
%  df   = [1 2 3];
%  cf   = @(t) cf_ChiSquare(t,df) ;
%  x    = linspace(0,40,101)';
%  prob = [0.9 0.95 0.99 0.999];
%  [result,cdf,pdf,qf] = cf2DistGPA(cf,x,prob);
%  % Comparison of the calculated and the true CDF values
%  disp([x cdf chi2cdf(x,6)])
%
% REFERENCES:
% [1] Gil-Pelaez, J., 1951. Note on the inversion theorem. Biometrika,
%     38(3-4), pp.481-482.
% [2] Imhof, J.: Computing the distribution of quadratic forms in normal
%     variables. Biometrika 48, 419–426 (1961).
% [3] Cohen, H., Villegas, F.R. and Zagier, D., 2000. Convergence
%     acceleration of alternating series. Experimental mathematics, 9(1),
%     pp.3-12.
% [4] Weniger, E.J., 1989. Nonlinear sequence transformations for the
%     acceleration of convergence and the summation of divergent series.
%     Computer Physics Reports, 10(5-6), pp.189-371. 
% [5] Sidi, A., 2011. A user-friendly extrapolation method for computing
%     infinite range integrals of products of oscillatory functions. IMA
%     Journal of Numerical Analysis, 32(2), pp.602-631.
% [6] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3),
%     32-44.
%
% SEE ALSO: cf2Dist, cf2DistGP, cf2DistGPT, cf2DistBTAV, cf2DistFFT,
%           cf2DistBV, cf2CDF_GPA, cf2PDF_GPA, cf2QF_GPA

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 03-Dec-2018 23:00:20

%% ALGORITHM
%[result,cdf,pdf,qf] = cf2DistGPA(cf,x,prob,options);

%% CHECK THE INPUT PARAMETERS
timeVal = tic;
narginchk(1, 4);

if nargin < 4, options = []; end
if nargin < 3, prob = []; end
if nargin < 2, x = []; end

if ~isfield(options, 'isCompound')
    options.isCompound = false;
end

if ~isfield(options, 'isCircular')
    options.isCircular = false;
end

if ~isfield(options, 'xMin')
    if options.isCompound
        options.xMin = 0;
    else
        options.xMin = -Inf;
    end
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

if ~isfield(options, 'qf0')
    options.qf0 = real((cf(1e-4)-cf(-1e-4))/(2e-4*1i));
end

if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
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

if ~isfield(options, 'nPeriods')
    options.nPeriods = 25;
end

if ~isfield(options, 'division')
    options.division = [5 8 11];
end

if ~isfield(options, 'isAccelerated')
    options.isAccelerated = true;
end

if ~isfield(options, 'firstRootID')
    options.firstRootID = 1;
end

if ~isfield(options, 'nPoly')
    options.nPoly = 2^5;
end

if ~isfield(options, 'tol')
    options.tol = 1e-10;
end

if ~isfield(options, 'tolFindRoots')
    options.tolFindRoots = 1e-32;
end

if ~isfield(options, 'maxiterFindRoots')
    options.maxiterFindRoots = 100;
end

if ~isfield(options, 'verbose')
    options.verbose = false;
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
    if options.isCircular
        % see https://en.wikipedia.org/wiki/Directional_statistics
        xMean = angle(cf(1));
    else
        xMean = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 ...
            - 2*cftIm(4)/280) / tolDiff;
    end
end

if isfinite(xMean)
    options.xMean = xMean;
end

if isempty(xStd)
    if options.isCircular
        % see https://en.wikipedia.org/wiki/Directional_statistics
        xStd  = sqrt(-2*log(abs(cf(1))));
    else
        xM2   = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 ...
            - 16*cftRe(3)/315 + 2*cftRe(4)/560) / tolDiff^2;
        xStd  = sqrt(xM2 - xMean^2);
    end
end

if options.isCircular
    xMin = -pi;
    xMax = pi;
else
    if isfinite(xMin)
        xMax = xMean + SixSigmaRule * xStd;
    elseif isfinite(xMax)
        xMin = xMean - SixSigmaRule * xStd;
    else
        xMin = xMean - SixSigmaRule * xStd;
        xMax = xMean + SixSigmaRule * xStd;
    end
end

if ~isfield(options, 'xN')
    options.xN = 101;
end

isPlot = options.isPlot;
options.isPlot = false;

%% ALGORITHM

if isempty(x)
    xempty = true;
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

cdf = cf2CDF_GPA(cf,x,options);
pdf = cf2PDF_GPA(cf,x,options);

if ~isempty(prob)
    qf = cf2QF_GPA(cf,prob,options);
else
    qf = [];
end
options.isPlot = isPlot;

%% Create the INTERPOLAN functions
if options.isInterp
    id   = isfinite(pdf);
    PDF  = @(xnew)PDFinterp(xnew,x(id),pdf(id));
    id   = isfinite(cdf);
    CDF  = @(xnew)CDFinterp(xnew,x(id),cdf(id));
    QF   = @(prob)QFinterp(prob,x(id),cdf(id));
    RND  = @(dim)RNDinterp(dim,x(id),cdf(id));
    try
    if ~xempty
        x   = x0;
        cdf = CDF(x);
        pdf = PDF(x);
    end
    catch
        warning('VW:CharFunTool:cf2DistGPA', ...
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
result.inversionMethod     = 'Gil-Pelaez';
result.quadratureMethod    = 'adaptive Gauss-Kronrod with acceleration'; 
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