function [result,cdf,pdf] = cf2Dist2D(cf,x,options)
%cf2Dist2D Calculates the CDF/PDF/QF from the BIVARIATE characteristic
%  function CF by using the Gil-Pelaez inversion formulae using Riemann sum,
%  as suggested in Shephard (1991).
%
%  The FOURIER INTEGRALs are calculated by using the simple RIEMANN SUM
%  QUADRATURE method. For more details see [1,2].
%
%  The algorithm cf2Dist2D is part of the MATLAB toolboc CharFunTool:
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
%             options.isCompound = false % treat the compound distributions
%                                        % of the RV Y = X_1 + ... + X_N,
%                                        % where N is discrete RV and X>=0
%                                        % are iid RVs from nonnegative
%                                        % continuous distribution.
%             options.isCircular = false % treat the circular
%                                        % distributions on (-pi,pi)
%             options.isInterp   = false % create and use the interpolant
%                                          functions for PDF/CDF/QF/RND
%             options.N = 2^8          % N points used quadrature rule
%             options.xMin = -Inf      % set the lower limit of X
%             options.xMax = Inf       % set the lower limit of X
%             options.xMean = []       % set the MEAN value of X
%             options.xStd = []        % set the STD value of X
%             options.dt = []          % set grid step dt = 2*pi/xRange
%             options.T = []           % set upper limit of (0,T), T = N*dt
%             options.SixSigmaRule = 6 % set the rule for computing domain
%             options.tolDiff = 1e-4   % tol for numerical differentiation
%             options.isPlot = true    % plot the graphs of PDF/CDF
%  options.DIST - structure with information for future evaluations.
%             options.DIST is created automatically after first call:
%             options.DIST.xMin  = xMin   % the lower limit of X
%             options.DIST.xMax  = xMax   % the upper limit of X
%             options.DIST.xMean = xMean  % the MEAN value of X,
%             options.DIST.cft   = cft    % CF evaluated at t_j : cf(t_j).
%
% REMARKS:
% If options.DIST is provided, then cf2Dist2D evaluates CDF/PDF based on
% this information only (it is useful, e.g., for subsequent evaluation of
% the quantiles). options.DIST is created automatically after first call.
% This is supposed to be much faster, bacause there is no need for further
% evaluations of the characteristic function. In fact, in such case the
% function handle of the CF is not required, i.e. in such case set cf = [];
%
% OUTPUT:
%  result   - structure with CDF/PDF/QF and further details,
%  cdf      - vector of CDF values evaluated at x,
%  pdf      - vector of PDF values evaluated at x,
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
% EXAMPLE1 (CDF/PDF of bivariate standard normal distribution)
%  cf = @(t) exp(-(0.9*t(:,1).^2 + 0.3*t(:,2).^2 +2*0.4*t(:,1).*t(:,2))/2);
%  result = cf2Dist2D(cf)
%  disp([result.x result.cdf])
%
% EXAMPLE2 (CDF/PDF of bivariate standard normal distribution)
%  cf = @(t) exp(-(0.9*t(:,1).^2 + 0.3*t(:,2).^2 +2*0.4*t(:,1).*t(:,2))/2);
%  x1 = linspace(-2,2,21);
%  x2 = linspace(-3,3,31);
%  x = {x1 x2};
%  result = cf2Dist2D(cf,x)
%  disp([result.x result.cdf])
%
% EXAMPLE3 (CDF/PDF of bivariate logistic distribution)
%  mu = [0 0];
%  beta = [1 2];
%  cf = @(t) cf2D_Logistic(t,mu,beta);
%  result = cf2Dist2D(cf)
%  disp([result.x result.cdf])
%
% EXAMPLE4 (CDF/PDF of bivariate logistic distribution)
%  mu = [0 2];
%  beta = [1 2];
%  cf = @(t) cf2D_Logistic(t,mu,beta);
%  x = {linspace(-5,5,11), linspace(-5,10,11)};
%  clear options;
%  options.isPlot = 'false';
%  result = cf2Dist2D(cf,x,[],options)
%  disp([result.x result.cdf])
%
% EXAMPLE5 (CDF/PDF of bivariate mixture of logistic distributions)
%  mu1 = [0 2];
%  beta1 = [1 2];
%  cf1 = @(t) cf2D_Logistic(t,mu1,beta1);
%  mu2 = [2 1];
%  beta2 = [2 1];
%  cf2 = @(t) cf2D_Logistic(t,mu2,beta2);
%  cf = @(t) 0.25*cf1(t) + 0.75*cf2(t);
%  clear options;
%  options.xN = 51;
%  result = cf2Dist2D(cf,[],[],options)
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
% Ver.: 13-Apr-2021 15:30:56

%% ALGORITHM
%[result,cdf,pdf] = cf2Dist2D(cf,x,options)

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

if ~isfield(options, 'N')
    if options.isCompound
        options.N = 2^8;
    else
        options.N = 2^6; % Set large N to improve the precision N = 2^10?
    end
end

if ~isfield(options, 'xMin')
    if options.isCompound
        options.xMin = [0 0];
    else
        options.xMin = [-Inf -Inf];
    end
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
    options.xN = 21;
end

if ~isfield(options, 'chebyPts')
    options.chebyPts = 2^9;
end

if ~isfield(options, 'correctedCDF')
    options.correctedCDF = false;
end

if ~isfield(options, 'isInterp')
    options.isInterp = false;
end
%% GET/SET the DEFAULT parameters and the OPTIONS
cfOld  = [];

if ~isempty(options.DIST)
    N                  = options.N;
    xMean              = options.DIST.xMean;
    cft                = options.DIST.cft;
    xMin               = options.DIST.xMin;
    xMax               = options.DIST.xMax;
    t                  = options.DIST.t;
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
    cft               = cf(([tolDiff;0]*(1:4))');
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
    options.DIST.t3    = t3;
    options.DIST.dt    = dt;
end

%% ALGORITHM
isPlot = options.isPlot;
if isempty(x)
    % Default values if x = [];
    %xempty = true;
    isMeshed = false;
    x1 = linspace(xMax(1),xMin(1),options.xN);
    x2 = linspace(xMax(2),xMin(2),options.xN);
    [X1,X2] = meshgrid(x1,x2);
    x = [X1(:) X2(:)];
else
    %xempty = false;
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
end

% Evaluate the required CDF/ PDF
% MARGINAL CDF1 and PDF1
x1      = x1(:);
n1      = length(x1);
E1      = exp(-1i*x1*t1');
cdf1    = imag(E1 * (cft1 ./ t1));
cdf1    = 0.5 - (cdf1 * dt(1)) / pi;
pdf1    = real(E1 * cft1);
pdf1    = (pdf1 * dt(1)) / pi;
pdf1    = max(0,pdf1);

% MARGINAL CDF2 and PDF2
x2      = x2(:);
n2      = length(x2);
E2      = exp(-1i*x2*t2');
cdf2    = imag(E2 * (cft2 ./ t2));
cdf2    = 0.5 - (cdf2 * dt(2)) / pi;
pdf2    = real(E2 * cft2);
pdf2    = (pdf2 * dt(2)) / pi;
pdf2    = max(0,pdf2);

% BIVARIATE CDF and PDF
if ~isMeshed
    [f1,f2] = meshgrid(cdf1,cdf2);
    cdf     = (f1 + f2)/2 - 0.25;
else
    cdf     = (cdf1 + cdf2)/2 - 0.25;
end
cdf     = cdf(:);
c       = -2 * dt(1) * dt(2) / (2*pi)^2;
cftt    = cft ./ t(:,1) ./ t(:,2);
f       = c * real(exp(-1i*x*t')*cftt);
cdf     = cdf + f;
cdf     = max(0,min(1,cdf));
pdf     = 2 * dt(1) * dt(2) * real(exp(-1i*x*t')*cft) / (2*pi)^2;
pdf     = max(0,pdf);

if ~isMeshed
    Zcdf    = reshape(cdf,n2,n1);
    Zpdf    = reshape(pdf,n2,n1);
else
    Zcdf    = cdf;
    Zpdf    = pdf;
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
result.prob                = prob;
result.cf                  = cfOld;
result.isCompound          = options.isCompound;
result.isCircular          = options.isCircular;
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