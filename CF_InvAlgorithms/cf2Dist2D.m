function [result,Zcdf,Zpdf] = cf2Dist2D(cf,x,options)
%cf2Dist2D Calculates the CDF/PDF/QF/RND from the BIVARIATE characteristic
%  function CF by using the Gil-Pelaez inversion formulae. The reuired
%  FOURIER INTEGRALs are calculated by using the simple RIEMANN SUM
%  QUADRATURE method, as suggested by Shephard, for more details see [1,2].
% 
%  The numerical inversion of the bivariate characteristic function is a
%  complex and computationally demanding task. In order to achieve adequate
%  calculation efficiency, it was necessary to reduce the demands on the
%  accuracy of the results. The algorithm with default settings provides
%  good results of calculated PDF and CDF for typical values of the
%  distribution, with possibly low accuracy for extreme values of the input
%  variable x = {x1, x2}. If necessary, it is advisable to change the
%  control parameters using the options.
%
%  The algorithm cf2Dist2D is part of the MATLAB toolbox CharFunTool:
%  https://github.com/witkovsky/CharFunTool
%
% SYNTAX:
%  result = cf2Dist2D(cf,x)
%  [result,cdf,pdf] = cf2Dist2D(cf,x,options)
%
% INPUT:
%  cf      - function handle of the bivariate characteristic function (CF),
%  x       - (N x 2)-matrix of x values where the CDF/PDF is computed,
%            or (1 x 2) cell array with x{1} being the vector of x1 values
%            and x{2} being the vector of x2 values. Then we shall
%            construct a meshgrid of x values, x = meshgrid(x{1},x{2}).
%            If x = [], the algorithm sets the automatic values, such that
%            x = meshgrid(x1,x2) covers the support of the distribution.
%  options - structure with the following default parameters:
%             options.isInterp = false % create and use the interpolant
%                                        functions for PDF/CDF/QF/RND
%                                        Selected interpolants - works
%                                        for marginal distributions
%             options.N = 2^8          % N points used quadrature rule
%             options.chebyPts = 2^6+1 % number of Chebyshev points used
%                                      % for interpolants
%             options.xMin = -Inf      % set the lower limit of X
%             options.xMax = Inf       % set the lower limit of X
%             options.xMean = []       % set the MEAN value of X
%             options.xStd = []        % set the STD value of X
%             options.dt = []          % set grid step dt = 2*pi/xRange
%             options.T = []           % set upper limit of (0,T), T = N*dt
%             options.SixSigmaRule = 6 % set the rule for computing domain
%             options.tolDiff = 1e-4   % tol for numerical differentiation
%             options.isPlot = true    % plot the graphs of PDF/CDF
%             options.ContourLevelList = [0.01 0.05 0.1:0.1:0.9 0.95 0.99] 
%                                      % defines the plotted countours at
%                                      % specified (relative) levels 
%  options.DIST - structure with information for future evaluations.
%             options.DIST is created automatically after first call:
%             options.DIST.xMin  = xMin   % the lower limit of X
%             options.DIST.xMax  = xMax   % the upper limit of X
%             options.DIST.xMean = xMean  % the MEAN value of X,
%             options.DIST.cft   = cft    % CF evaluated at t_j : cf(t_j).
%             options.DIST.cft1  = cft1;  % marginal CF evaluated at t1.
%             options.DIST.cft2  = cft2;  % marginal CF evaluated at t2.
%             options.DIST.N     = N;     % N points used quadrature rule
%             options.DIST.t     = t;     % t pairs of values where cf was
%                                         % evaluated
%             options.DIST.t1    = t1;    % t1 values where cf1 was
%                                         % evaluated
%             options.DIST.t2    = t2;    % t2 values where cf2 was
%                                         % evaluated
%             options.DIST.dt    = dt;    % pair of delta difference in t.
%
% REMARKS: 
% If options.DIST is provided, then cf2Dist2D evaluates CDF/PDF based on
% this information only. options.DIST is created automatically after first
% call. This is supposed to be much faster, bacause there is no need for
% further evaluations of the characteristic function. In fact, in such case
% the function handle of the CF is not required, i.e. in such case set cf =
% [];
%
% OUTPUT:
%  result   - structure with CDF/PDF/QF and further details,
%  Zcdf     - vector or matrix of CDF values evaluated at specified x,
%  Zpdf     - vector or matrix of PDF values evaluated at specified x,
%
% REMARKS:
% The required integrals are evaluated approximately by using the simple
% quadrature using Riemann sum rules on the intervals (0,T) x (-T,T), where
% T = N * dt is a sufficienly large integration upper limit in the
% frequency domain.
%
% If the optimum values of N and T are unknown, we suggest, as a simple
% rule of thumb, to start with the application of the six-sigma-rule for
% determining the value of dt = (2*pi)/(xMax-xMin), where xMax = xMean +
% 6*xStd, and xMin = xMean - 6*xStd.
%
% Please note that THIS QUADRATURE METHOD IS AN APPROXIMATE METHOD:
% Frequently, with relatively low numerical precision of the results of the
% calculated PDF/CDF/QF, but frequently more efficient and more precise
% than comparable Monte Carlo methods.
%
% However, the numerical error (truncation error and/or the integration
% error) could be and should be properly controled! See [1,2].
%
% EXAMPLE1 (CDF/PDF of bivariate normal distribution)
%  cf = @(t) exp(-(0.9*t(:,1).^2 + 0.3*t(:,2).^2 +2*0.4*t(:,1).*t(:,2))/2);
%  result = cf2Dist2D(cf)
%  disp([result.x result.cdf])
%
% EXAMPLE2 (CDF/PDF of bivariate normal distribution at specific x)
%  cf = @(t) exp(-(0.9*t(:,1).^2 + 0.3*t(:,2).^2 +2*0.4*t(:,1).*t(:,2))/2);
%  x1 = linspace(-3,3,31);
%  x2 = linspace(-3,3,31);
%  x = {x1 x2};
%  result = cf2Dist2D(cf,x)
%  disp([result.x result.cdf])
%
% EXAMPLE3 (Iterpolated CDF/PDF of the standard bivariate logistic distribution)
%  cf = @(t) cf2D_Logistic(t);
%  clear options;
%  options.isInterp = true;
%  result = cf2Dist2D(cf,[],options)
%
% EXAMPLE4 (Iterpolated COPULA of the standard bivariate logistic distribution)
%  cf = @(t) cf2D_Logistic(t);
%  clear options;
%  options.isInterp = true;
%  options.isPlot = false;
%  result = cf2Dist2D(cf,[],options);
%  COPULAcdf = result.COPULAcdf;
%  COPULApdf = result.COPULApdf;
%  u1 = linspace(0,1,21);
%  u2 = linspace(0,1,21);
%  CopulaCDF = COPULAcdf({u1,u2});
%  CopulaPDF = COPULApdf({u1,u2});
%  [U1,U2]= meshgrid(u1,u2);
%  figure
%  mesh(U1,U2,CopulaCDF)
%  xlabel('u1')
%  ylabel('u2')
%  zlabel('CDF')
%  title('CDF COPULA of the Bivariate Distribution Specified by CF')
%  figure
%  mesh(U1,U2,CopulaPDF)
%  xlabel('u1')
%  ylabel('u2')
%  zlabel('PDF')
%  title('PDF COPULA of the Bivariate Distribution Specified by CF')
%
% EXAMPLE5 (CDF/PDF of bivariate mixture of logistic distributions)
%  mu1    = [1 2];
%  sigma1 = [1 2];
%  cf1    = @(t) cf2D_Logistic(t,mu1,sigma1);
%  mu2    = [-1 0];
%  sigma2 = [2 1];
%  cf2    = @(t) cf2D_Logistic(t,mu2,sigma2);
%  cf     = @(t) 0.25*cf1(t) + 0.75*cf2(t);
%  clear options;
%  options.isInterp = true;
%  options.N = 2^9;          
%  options.chebyPts = 2^7;
%  options.ContourLevelList = [0.01 0.05 0.1 0.2:0.2:0.8 0.9 0.95 0.99];
%  result = cf2Dist2D(cf,[],options)
%
% REFERENCES:
%
% [1] SHEPHARD, N.G., 1991. Numerical integration rules for multivariate
%     inversions. Journal of Statistical Computation and Simulation,
%     39(1-2), pp.37-46.
% [2] SHEPHARD, N.G., 1991. From characteristic function to distribution
%     function: a simple framework for the theory. Econometric theory,
%     pp.519-529.
%
% SEE ALSO: cf2Dist, cf2DistGP, cf2Dist2D, cf2DistGPA, cf2DistFFT,
%           cf2DistBV, cf2CDF, cf2PDF, cf2QF

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 06-Jan-2022 16:55:48
% Ver.: 04-Jan-2022 19:57:42
% Ver.: 09-Dec-2021 18:07:29
% Ver.: 13-May-2021 17:00:00

%% ALGORITHM
%[result,cdf,pdf] = cf2Dist2D(cf,x,options)

%% CHECK THE INPUT PARAMETERS
timeVal = tic;
narginchk(1, 3);

if nargin < 3, options = []; end
if nargin < 2, x = []; end

if ~isfield(options, 'N')
    options.N = 2^8; % Set large N to improve the precision, e.g. N = 2^10
end

if ~isfield(options, 'xMin')
    options.xMin = [-Inf -Inf];
end

if ~isfield(options, 'xMax')
    options.xMax = [Inf Inf];
end

if ~isfield(options, 'xMean')
    options.xMean = [];
end

if ~isfield(options, 'xStd')
    options.xStd = [];
end

if ~isfield(options, 'dt')
    options.dt = [];
end

if ~isfield(options, 'T')
    options.T = [];
end

if ~isfield(options, 'SixSigmaRule')
    options.SixSigmaRule = 6;
end

if ~isfield(options, 'tolDiff')
    options.tolDiff = 1e-4;
end

if ~isfield(options, 'crit')
    options.crit = 1e-12;
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

if ~isfield(options, 'isForcedPlot')
    options.isForcedPlot = false;
end

if ~isfield(options, 'DIST')
    options.DIST = [];
end

if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end

if ~isfield(options, 'xN')
    options.xN = 51;
end

if ~isfield(options, 'chebyPts')
    options.chebyPts = 2^6+1;
end

if ~isfield(options, 'correctedCDF')
    options.correctedCDF = false;
end

if ~isfield(options, 'isInterp')
    options.isInterp = false;
end

if ~isfield(options, 'cftTol')
    options.cftTol = 1e-14;       
end

if ~isfield(options, 'ContourLevelList')
    options.ContourLevelList = [0.01 0.05 0.1000 0.2000 0.3000 ...
        0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 0.95 0.99];       
end

%% GET/SET the DEFAULT parameters and the OPTIONS
cfOld  = [];

if ~isempty(options.DIST)
    N                  = options.N;
    xMean              = options.DIST.xMean;
    cft                = options.DIST.cft;
    cft1               = options.DIST.cft1;
    cft2               = options.DIST.cft2;
    xMin               = options.DIST.xMin;
    xMax               = options.DIST.xMax;
    t                  = options.DIST.t;
    t1                 = options.DIST.t1;
    t2                 = options.DIST.t2;
    dt                 = options.DIST.dt;
    xStd               = [];   
else
    N                  = options.N;
    dt                 = options.dt;
    T                  = options.T;
    xMin               = options.xMin;
    xMax               = options.xMax;
    xMean              = options.xMean;
    xStd               = options.xStd;
    SixSigmaRule       = options.SixSigmaRule;
    tolDiff            = options.tolDiff;
    cfOld              = cf;
    cft                = cf(([tolDiff;0]*(1:4))');
    cftRe1             = real(cft);
    cftIm1             = imag(cft);
    cft2               = cf(([0;tolDiff]*(1:4))');
    cftRe2             = real(cft2);
    cftIm2             = imag(cft2);
    if isempty(xMean)
        xMean(1) = (8*cftIm1(1)/5 - 2*cftIm1(2)/5 + 8*cftIm1(3)/105 ...
            - 2*cftIm1(4)/280) / tolDiff;
        xMean(2) = (8*cftIm2(1)/5 - 2*cftIm2(2)/5 + 8*cftIm2(3)/105 ...
            - 2*cftIm2(4)/280) / tolDiff;
    end
    if isempty(xStd)
        xM2(1)   = (205/72 - 16*cftRe1(1)/5 + 2*cftRe1(2)/5 ...
            - 16*cftRe1(3)/315 + 2*cftRe1(4)/560) / tolDiff^2;
        xStd(1)  = sqrt(xM2(1) - xMean(1)^2);
        xM2(2)   = (205/72 - 16*cftRe2(1)/5 + 2*cftRe2(2)/5 ...
            - 16*cftRe2(3)/315 + 2*cftRe2(4)/560) / tolDiff^2;
        xStd(2)  = sqrt(xM2(2) - xMean(2)^2);
    end
    if all(isfinite(xMin)) && all(isfinite(xMax))
        range          = xMax - xMin;
    elseif ~isempty(dt)
        range = 2*pi ./ dt;
        if all(isfinite(xMin))
            xMax = xMin + range;
        elseif all(isfinite(xMax))
            xMin = xMax - range;
        else
            xMax = xMean + range/2;
            xMin = xMean - range/2;
        end
    elseif ~isempty(T)
        range = 2*pi ./ (T / N);
        if all(isfinite(xMin))
            xMax = xMin + range;
        elseif all(isfinite(xMax))
            xMin = xMax - range;
        else
            xMax = xMean + range/2;
            xMin = xMean - range/2;
        end
    else
        if all(isfinite(xMin))
            xMax = xMean + 2*SixSigmaRule * xStd;
        elseif all(isfinite(xMax))
            xMin = xMean - 2*SixSigmaRule * xStd;
        else
            xMin = xMean - SixSigmaRule * xStd;
            xMax = xMean + SixSigmaRule * xStd;
        end
        range = xMax - xMin;
    end
    dt                 = 2*pi ./ range;
    t1                 = (0.5+(0:N))'*dt(1);
    t2                 = (0.5+(0:N))'*dt(2);
    t3                 = (0.5+(-N:N))'*dt(2);
    cft1               = cf([t1 0*t1]);
    cft2               = cf([0*t2 t2]);
    [tt1,tt2]          = meshgrid(t1,t3);
    t                  = [tt1(:) tt2(:)];
    cft                = cf(t);
    options.DIST.xMin  = xMin;
    options.DIST.xMax  = xMax;
    options.DIST.xMean = xMean;
    options.DIST.cft   = cft;
    options.DIST.cft1  = cft1;
    options.DIST.cft2  = cft2;
    options.DIST.N     = N;
    options.DIST.t     = t;
    options.DIST.t1    = t1;
    options.DIST.t2    = t2;
    options.DIST.dt    = dt;
end

%% ALGORITHM
cftTol = options.cftTol;
isPlot = options.isPlot;
if isempty(x)  
    if options.isInterp
        % Chebyshev points if x = [];
        isMeshed = false;
        x1 = (xMax(1)-xMin(1)) * (-cos(pi*(0:(options.chebyPts-1)) / ...
            (options.chebyPts-1)) + 1) / 2 + xMin(1);
        x2 = (xMax(2)-xMin(2)) * (-cos(pi*(0:(options.chebyPts-1)) / ...
            (options.chebyPts-1)) + 1) / 2 + xMin(2);
        [X1,X2] = meshgrid(x1,x2);
        x = [X1(:) X2(:)];
    else
        % Default values if x = [];
        isMeshed = false;
        x1 = linspace(xMin(1),xMax(1),options.xN);
        x2 = linspace(xMin(2),xMax(2),options.xN);
        [X1,X2] = meshgrid(x1,x2);
        x = [X1(:) X2(:)];
    end    
else
    if iscell(x)
        isMeshed = false;
        x1 = x{1};
        x2 = x{2};
        [X1,X2] = meshgrid(x1,x2);
        x = [X1(:) X2(:)];
    else
        isMeshed = true;
        x1 = x(:,1);
        x2 = x(:,2);
        X1 = [];
        X2 = [];
        isPlot = false;
    end
    if any(x1 < xMin(1)) || any(x1 > xMax(1))
        warning('VW:CharFunTool:cf2Dist2D',['x1 out-of-range: ', ...
            '[x1Min, x1Max] = [',num2str(xMin(1)),...
            ', ',num2str(xMax(1)),'] !']);
    end
    if any(x2 < xMin(2)) || any(x2 > xMax(2))
        warning('VW:CharFunTool:cf2Dist2D',['x2 out-of-range: ', ...
            '[x2Min, x2Max] = [',num2str(xMin(2)),...
            ', ',num2str(xMax(2)),'] !']);
    end
end

% MARGINAL CDF1 and PDF1
x1      = x1(:);
n1      = length(x1);
id      = abs(cft1) > cftTol;
cft1    = cft1(id);
t1      = t1(id);
options.DIST.cft1  = cft1;
options.DIST.t1    = t1;

% MARGINAL PDF1
E1      = exp(-1i*x1*t1');
pdf1    = real(E1 * cft1);
pdf1    = (pdf1 * dt(1)) / pi;
pdf1    = max(0,pdf1);

% MARGINAL CDF1
cftt1   = cft1 ./ t1;
cdf1    = imag(E1 * cftt1);
cdf1    = 0.5 - (cdf1 * dt(1)) / pi;

% MARGINAL CDF2 and PDF2
x2      = x2(:);
n2      = length(x2);
id      = abs(cft2) > cftTol;
cft2    = cft2(id);
t2      = t2(id);
options.DIST.cft2  = cft2;
options.DIST.t2    = t2;

% MARGINAL PDF2
E2      = exp(-1i*x2*t2');
pdf2    = real(E2 * cft2);
pdf2    = (pdf2 * dt(2)) / pi;
pdf2    = max(0,pdf2);

% MARGINAL CDF2
cftt2   = cft2 ./ t2;
cdf2    = imag(E2 * cftt2);
cdf2    = 0.5 - (cdf2 * dt(2)) / pi;

% BIVARIATE CDF and PDF
c       = 2 * dt(1) * dt(2) / (2*pi)^2;
id      = abs(cft) > cftTol;
cft     = cft(id);
t       = t(id,:);
options.DIST.cft   = cft;
options.DIST.t     = t;

% BIVARIATE PDF
E       = exp(-1i*x*t');
pdf     = c * real(E*cft);
pdf     = max(0,pdf);

% BIVARIATE CDF
if ~isMeshed
    [f1,f2] = meshgrid(cdf1,cdf2);
    cdf     = (f1 + f2)/2 - 0.25;
else
    cdf     = (cdf1 + cdf2)/2 - 0.25;
end
cdf     = cdf(:);
cftt    = cft ./ t(:,1) ./ t(:,2);
cdf     = cdf - c * real(E*cftt);
cdf     = max(0,min(1,cdf));

% RESHAPE to the MESHGRID
if ~isMeshed
    Zcdf    = reshape(cdf,n2,n1);
    Zpdf    = reshape(pdf,n2,n1);
else
    Zcdf    = cdf;
    Zpdf    = pdf;
end

% Create the INTERPOLANTS - Interpolation Functions
if options.isInterp
    % INTERPOLANTS of the MARGINAL PDF / CDF / QF / RND
    PDF1  = @(x1new) InterpPDF(x1new,x1,pdf1);
    CDF1  = @(x1new) InterpCDF(x1new,x1,cdf1);
    QF1   = @(prob) InterpQF(prob,x1,cdf1);
    RND1  = @(dim) InterpRND(dim,x1,cdf1);
    PDF2  = @(x2new) InterpPDF(x2new,x2,pdf2);
    CDF2  = @(x2new) InterpCDF(x2new,x2,cdf2);
    QF2   = @(prob) InterpQF(prob,x2,cdf2);
    RND2  = @(dim) InterpRND(dim,x2,cdf2);
    % INTERPOLANTS of the CONDITIONAL PDF / CDF / QF / RND
%     PDF12  = @(x1new,x2fix) InterpPDF12(x1new,x2fix,x1,x2,t,cft,dt,PDF2);
%     CDF12  = @(x1new,x2fix) InterpCDF12(x1new,x2fix,x1,x2,t,cft,dt,cdf1,CDF2);
%     QF12   = @(prob,x2fix) InterpQF(prob,x1,CDF12(x1,x2fix));
%     RND12  = @(dim)InterpRND (dim,x1,CDF12(x1,x2fix));
%     PDF21  = @(x2new,x1fix) InterpPDF21(x1fix,x2new,x1,x2,t,cft,dt,PDF1);
%     CDF21  = @(x2new,x1fix) InterpCDF21(x1fix,x2new,x1,x2,t,cft,dt,cdf2,CDF1);
%     QF21   = @(prob,x1fix)InterpQF(prob,x2,CDF21(x2,x1fix));
%     RND21  = @(dim)InterpRND(dim,x2,CDF21(x2,x1fix));
    % INTERPOLANTS of the BIVARIATE PDF / CDF / RND
    PDF    = @(xyNew) InterpBarycentric2D(x1,x2,Zpdf,xyNew);
    CDF    = @(xyNew) InterpBarycentric2D(x1,x2,Zcdf,xyNew);
%     RND    = @(dim) InterpRND2D(dim,RND1,RND2,PDF1,PDF2,PDF,xMin,xMax);
    COPULAcdf = @(u12) CopFunCDF(u12,x1,x2,Zcdf,QF1,QF2);
    COPULApdf = @(u12) CopFunPDF(u12,x1,x2,Zpdf,QF1,QF2);
else
    PDF1   = [];
    CDF1   = [];
    QF1    = [];
    PDF2   = [];
    CDF2   = [];
    QF2    = [];
%     PDF12  = [];
%     CDF12  = [];
%     QF12   = [];
%     PDF21  = [];
%     CDF21  = [];
%     QF21   = [];
    PDF    = [];
    CDF    = [];
%     RND    = [];
    COPULAcdf = [];
    COPULApdf = [];
end

PrecisionCrit = abs(cft(end)/t(end));
isPrecisionOK = (PrecisionCrit<=options.crit);

%% RESULT
result.Description         = 'CDF/PDF/QF from the characteristic function CF';
result.inversionMethod     = 'Gil-Pelaez';
result.quadratureMethod    = 'Riemann sums quadrature rule';
result.x                   = x;
result.cdf                 = cdf;
result.pdf                 = pdf;
result.x1                  = x1;
result.cdf1                = cdf1;
result.pdf1                = pdf1;
result.x2                  = x2;
result.cdf2                = cdf2;
result.pdf2                = pdf2;
result.X1                  = X1;
result.X2                  = X2;
result.Zcdf                = Zcdf;
result.Zpdf                = Zpdf;
result.cf                  = cfOld;
if options.isInterp
    result.PDF             = PDF;
    result.CDF             = CDF;
%    result.RND             = RND;
    result.COPULAcdf       = COPULAcdf;
    result.COPULApdf       = COPULApdf;
    result.PDF1            = PDF1;
    result.CDF1            = CDF1;
    result.QF1             = QF1;
    result.RND1            = RND1;
    result.PDF2            = PDF2;
    result.CDF2            = CDF2;
    result.QF2             = QF2;
    result.RND2            = RND2;
%     result.PDF12           = PDF12;
%     result.CDF12           = CDF12;
%     result.QF12            = QF12;
%     result.RND12           = RND12;
%     result.PDF21           = PDF21;
%     result.CDF21           = CDF21;
%     result.QF21            = QF21;
%     result.RND21           = RND21;
end
result.isInterp            = options.isInterp;
result.SixSigmaRule        = options.SixSigmaRule;
result.N                   = N;
result.dt                  = dt;
result.T                   = t(end);
result.PrecisionCrit       = PrecisionCrit;
result.myPrecisionCrit     = options.crit;
result.isPrecisionOK       = isPrecisionOK;
result.xMean               = xMean;
result.xStd                = xStd;
result.xMin                = xMin;
result.xMax                = xMax;
result.options             = options;
result.tictoc              = toc(timeVal);

%% PLOT the PDF / CDF
if options.isForcedPlot
    isPlot = true;
end
if length(x)==1
    isPlot = false;
end

if isPlot
    if options.isInterp
        x1 = linspace(xMin(1),xMax(1),101);
        x2 = linspace(xMin(2),xMax(2),101);
        % Marginal PDF1
        figure
        plot(x1,PDF1(x1),'LineWidth',2)
        grid
        title('Marginal PDF1 Specified by the CF')
        xlabel('x1')
        ylabel('pdf1')
        
        % Marginal CDF1
        figure
        plot(x1,CDF1(x1),'LineWidth',2)
        grid
        title('Marginal CDF1 Specified by the CF')
        xlabel('x1')
        ylabel('cdf1')
        
        % Marginal PDF2
        figure
        plot(x2,PDF2(x2),'LineWidth',2)
        grid
        title('Marginal PDF2 Specified by the CF')
        xlabel('x2')
        ylabel('pdf2')
        
        % Marginal CDF2
        figure
        plot(x2,CDF2(x2),'LineWidth',2)
        grid
        title('Marginal CDF2 Specified by the CF')
        xlabel('x2')
        ylabel('cdf2')
        
        % PDF
        [X1,X2] = meshgrid(x1,x2);
        figure
        mesh(X1,X2,PDF({x1,x2}))
        title('PDF Specified by the CF')
        xlabel('x1')
        ylabel('x2')
        zlabel('pdf')
        
        % CDF
        figure
        mesh(X1,X2,CDF({x1,x2}))
        title('CDF Specified by the CF')
        xlabel('x1')
        ylabel('x2')
        zlabel('cdf')
               
        % CDF COPULA
        figure
        u1 = linspace(0,1,51);
        u2 = linspace(0,1,51);
        [U1,U2]= meshgrid(u1,u2);
        mesh(U1,U2,COPULAcdf({u1,u2}))
        xlabel('u1')
        ylabel('u2')
        zlabel('CDF')
        title('COPULA of the Bivariate Distribution Specified by the CF')
        
        % PDF COPULA
        figure
        mesh(U1,U2,COPULApdf({u1,u2}))
        xlabel('u1')
        ylabel('u2')
        zlabel('PDF')
        title('PDF COPULA of the Bivariate Distribution Specified by the CF')
        
        % Contour plot of PDF + CDF
        figure
        [~,c1] = contour(X1,X2,PDF({x1,x2}),'ShowText','on');
        c1.LineWidth = 2;
        c1.LevelList = options.ContourLevelList*max(max(Zpdf));
        hold on
        [~,c] = contour(X1,X2,CDF({x1,x2}),'ShowText','on');
        c.LineStyle = '--';
        c.LineWidth = 2;
        c.LevelList = options.ContourLevelList;
        grid
        hold off
        title('Contour Plot of the PDF and the CDF Specified by the CF')
        xlabel('x1')
        ylabel('x2')
    else
        % Marginal PDF1
        figure
        plot(x1,pdf1,'LineWidth',2)
        grid
        title('Marginal PDF1 Specified by the CF')
        xlabel('x1')
        ylabel('pdf1')
        
        % Marginal CDF1
        figure
        plot(x1,cdf1,'LineWidth',2)
        grid
        title('Marginal CDF1 Specified by the CF')
        xlabel('x1')
        ylabel('cdf1')
        
        % Marginal PDF2
        figure
        plot(x2,pdf2,'LineWidth',2)
        grid
        title('Marginal PDF2 Specified by the CF')
        xlabel('x2')
        ylabel('pdf2')
        
        % Marginal CDF2
        figure
        plot(x2,cdf2,'LineWidth',2)
        grid
        title('Marginal CDF2 Specified by the CF')
        xlabel('x2')
        ylabel('cdf2')
        
        % PDF
        figure
        mesh(X1,X2,Zpdf)
        title('PDF Specified by the CF')
        xlabel('x1')
        ylabel('x2')
        zlabel('pdf')
        
        % CDF
        figure
        mesh(X1,X2,Zcdf)
        title('CDF Specified by the CF')
        xlabel('x1')
        ylabel('x2')
        zlabel('cdf')
        
        % Contour plot of PDF + CDF
        figure
        [~,c1] = contour(X1,X2,Zpdf,'ShowText','on');
        c1.LineWidth = 2;
        hold on
        [~,c] = contour(X1,X2,Zcdf,'ShowText','on');
        c.LineStyle = '--';
        c.LineWidth = 2;
        c.LevelList = [0.01 0.05 0.1000 0.2000 0.3000 0.4000 0.5000 ...
            0.6000 0.7000 0.8000 0.9000 0.95 0.99];
        grid
        hold off
        title('Contour Plot of the PDF and the CDF Specified by the CF')
        xlabel('x1')
        ylabel('x2')
    end
end
end
% %% Function InterpPDF12
% function pdf = InterpPDF12(x1New,x2Given,x1,x2,t,cft,dt,PDF2)
% % InterpPDF12 Auxiliary function evaluates the conditional PDF for x1New
% % and given x2Given.
% 
% % (c) Viktor Witkovsky (witkovsky@gmail.com)
% % Ver.: 01-May-2021 14:20:34
% 
% %% ALGORITHM
% if isempty(x1New)
%     x1New = x1;
% end
% 
% if isempty(x2Given)
%     x2Given = mean(x2);
% end
% 
% % PDF12a = @(x2fix) (2 * dt(1) * dt(2) / (2*pi)^2 / PDF2(x2fix)) ...
% %     * max(0,real(exp(-1i*[x1 x2fix*ones(size(x1))]*t')*cft));
% % PDF12 =@(x1new,x2fix) InterpPDF(x1new,x1,PDF12a(x2fix));
% 
% % WARNING: Out-of-range
% xMin = min(x1);
% xMax = max(x1);
% if any(x1New < xMin) || any(x1New > xMax)
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x1Min, x1Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% xMin = min(x2);
% xMax = max(x2);
% if x2Given < xMin || x2Given > xMax
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x2Min, x2Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% 
% 
% pdf = InterpPDF(x1New,x1,2*dt(1)*dt(2)/(2*pi)^2/PDF2(x2Given)* ...
%     max(0,real(exp(-1i*[x1 x2Given*ones(size(x1))]*t')*cft)));
% end
% %% Function InterpPDF21
% function pdf = InterpPDF21(x1Given,x2New,x1,x2,t,cft,dt,PDF1)
% % InterpPDF21 Auxiliary function evaluates the conditional PDF for x2New
% % and given x1Given.
% 
% % (c) Viktor Witkovsky (witkovsky@gmail.com)
% % Ver.: 01-May-2021 14:20:34
% 
% %% ALGORITHM
% if isempty(x2New)
%     x2New = x2;
% end
% 
% if isempty(x1Given)
%     x1Given = mean(x1);
% end
% % PDF21a = @(x1fix) (2 * dt(1) * dt(2) / (2*pi)^2 / PDF1(x1fix)) ...
% %     * max(0,real(exp(-1i*[x1fix*ones(size(x2)) x2]*t')*cft));
% % PDF21 =@(x2new,x1fix) InterpPDF(x2new,x2,PDF21a(x1fix));
% 
% % WARNING: Out-of-range
% xMin = min(x2);
% xMax = max(x2);
% if any(x2New < xMin) || any(x2New > xMax)
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x2Min, x2Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% xMin = min(x1);
% xMax = max(x1);
% if x1Given < xMin || x1Given > xMax
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x1Min, x1Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% 
% pdf = InterpPDF(x2New,x2,2*dt(1)*dt(2)/(2*pi)^2/PDF1(x1Given)* ...
%     max(0,real(exp(-1i*[x1Given*ones(size(x2)) x2]*t')*cft)));
% end
% %% Function InterpCDF12
% function cdf = InterpCDF12(x1New,x2Given,x1,x2,t,cft,dt,cdf1,CDF2)
% % InterpCDF12 Auxiliary function evaluates the conditional CDF for x1New
% % and given x2Given.
% 
% % (c) Viktor Witkovsky (witkovsky@gmail.com)
% % Ver.: 01-May-2021 14:20:34
% 
% %% ALGORITHM
% if isempty(x1New)
%     x1New = x1;
% end
% 
% if isempty(x2Given)
%     x2Given = mean(x2);
% end
% 
% % WARNING: Out-of-range
% xMin = min(x1);
% xMax = max(x1);
% if any(x1New < xMin) || any(x1New > xMax)
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x1Min, x1Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% xMin = min(x2);
% xMax = max(x2);
% if x2Given < xMin || x2Given > xMax
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x2Min, x2Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% 
% cft = cft ./ t(:,1) ./ t(:,2);
% c   = -2 * dt(1) * dt(2) / (2*pi)^2;
% f   = (cdf1 + CDF2(x2Given))/2 - 0.25;
% f   = f + c * real(exp(-1i*[x1 x2Given*ones(size(x1))]*t')*cft);
% f   = f / CDF2(x2Given);
% cdf = InterpCDF(x1New,x1,max(0,min(1,f)));
% 
% end
% %% Function InterpCDF21
% function cdf = InterpCDF21(x1Given,x2New,x1,x2,t,cft,dt,cdf2,CDF1)
% % InterpCDF21 Auxiliary function evaluates the conditional CDF for x2New
% % and given x1Given.
% 
% % (c) Viktor Witkovsky (witkovsky@gmail.com)
% % Ver.: 01-May-2021 14:20:34
% 
% %% ALGORITHM
% if isempty(x2New)
%     x2New = x2;
% end
% 
% if isempty(x1Given)
%     x1Given = mean(x1);
% end
% 
% % WARNING: Out-of-range
% xMin = min(x2);
% xMax = max(x2);
% if any(x2New < xMin) || any(x2New > xMax)
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x2Min, x2Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% xMin = min(x1);
% xMax = max(x1);
% if x1Given < xMin || x1Given > xMax
%     warning('VW:CharFunTool:cf2Dist2D',['x out-of-range: ', ...
%         '[x1Min, x1Max] = [',num2str(xMin),...
%         ', ',num2str(xMax),'] !']);
% end
% 
% cft = cft ./ t(:,1) ./ t(:,2);
% c   = -2 * dt(1) * dt(2) / (2*pi)^2;
% f   = (cdf2 + CDF1(x1Given))/2 - 0.25;
% %f   = (cdf1 + cdf2)/2 - 0.25;
% f   = f + c * real(exp(-1i*[x1Given*ones(size(x2)) x2]*t')*cft);
% f   = f / CDF1(x1Given);
% cdf = InterpCDF(x2New,x2,max(0,min(1,f)));
% 
% end
% %% Function InterpRND2D
% function xRND = InterpRND2D(N,RND1,RND2,PDF1,PDF2,PDF,xMin,xMax)
% % InterpRND2D Auxiliary function generates the (Nx2)-dimensional matrix of
% % random pairs x = [x1,x2] generated from the bivariate distribution
% % specified by the characteristic function.  
% 
% % (c) Viktor Witkovsky (witkovsky@gmail.com)
% % Ver.: 12-May-2021 16:00:38
% %% ALGORITHM
% 
% range1 = xMax(1) - xMin(1);
% range2 = xMax(2) - xMin(2);
% 
% % Set myN as some multiple of given N 
% myN = ceil(N*2);
% 
% % Generate M = 2*myN candidate samples from the joint bivariate
% % distribution created from the equally proportional mixture of joint
% % distribution generated from independent marginals and the joint uniform
% % distribution  
% xCandidate     = [RND1(myN) RND2(myN); ...
%                   rand(myN,1)*range1 + xMin(1), ...
%                   rand(myN,1)*range2 + xMin(2)];
% 
% % Evaluate the true PDF values for the candidate samples
% pdf            = PDF(xCandidate);
% maxPDF         = max(max(pdf));
% 
% % Evaluate the PDF values for the candidate samples using the candidate PDF
% % created from the independent marginals
% pdf12Mix      = (PDF1(xCandidate(:,1)).*PDF2(xCandidate(:,2)) + ...
%     1/(range1*range2))/2;
% maxPdf12Mix   = max(max(pdf12Mix));
%  
% 
% % Premultiply the candidate PDF by suitable constant such that the
% % candidate PDF are greater than the required PDF values.
% pdf12Mix = maxPDF * pdf12Mix / maxPdf12Mix;
% 
% % Acceptance-rejection method. Select such values that fulfill the
% % required relation of the candidate PDF and the required PDF 
% idHappy        = rand(2*myN,1) .* pdf12Mix <= pdf;
% xRND           = xCandidate(idHappy,:);
% 
% % If the number of correct values is greater than required, save only N.
% % Otherwise print the warning message
% NHappy         = sum(idHappy);
% if NHappy >= N
%     xRND = xRND(randperm(NHappy,N),:);
% else
%     warning(['The number of generated samples ',num2str(NHappy), ...
%         ' is less than required ',num2str(N)])
% end
% 
% % This is the ration of situation where the candidate PDF sample has lower
% % probability than the required PDF sample. This should be 0 if we have
% % good candidate. 
% %UnderRatio = sum(pdf12Mix - pdf < 0)/(2*N)
% %disp([xCandidate pdf12Mix pdf (pdf12Mix-pdf)< 0])
% end
% 
%% Function CopFunCDF
function cdf = CopFunCDF(u12,x1,x2,Zcdf,QF1,QF2)
% CopFunCDF for given u1 in (0,1) and u2 in (0,1) Evaluates the copula
% values of the distribution generated by the bivariate CF  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 12-May-2021 16:00:38
%% ALGORITHM

if iscell(u12)
    u1     = u12{1};
    u2     = u12{2};
    x1New  = QF1(u1);
    x2New  = QF2(u2);
    x12New = {x1New x2New};
else
    x12New = [QF1(u12(:,1)) QF2(u12(:,2))];
end

cdf = InterpBarycentric2D(x1,x2,Zcdf,x12New);

end
%% Function CopFunPDF
function pdf = CopFunPDF(u12,x1,x2,Zpdf,QF1,QF2)
% CopFunPDF for given u1 in (0,1) and u2 in (0,1) Evaluates the "copula"
% values of PDF of the distribution generated by the bivariate CF  

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 12-May-2021 16:00:38
%% ALGORITHM


if iscell(u12)
    u1     = u12{1};
    u2     = u12{2};
    x1New  = QF1(u1);
    x2New  = QF2(u2);
    x12New = {x1New x2New};
else
    x12New = [QF1(u12(:,1)) QF2(u12(:,2))];
end

pdf = InterpBarycentric2D(x1,x2,Zpdf,x12New);

end