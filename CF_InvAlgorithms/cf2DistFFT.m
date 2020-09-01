function [result,cdf,pdf,qf] = cf2DistFFT(cf,x,prob,options)
%cf2DistFFT(cf,x,prob,options) Evaluates the CDF/PDF/QF (quantiles) from
% the characteristic function CF of a (continuous) distribution F by using
% the Fast Fourier Transform (FFT) algorithm. 
%
% The algorithm cf2DistFFT evaluates the approximate values CDF(x), PDF(x),
% and/or the quantiles QF(prob) for given x and prob, by interpolation from
% the PDF-estimate computed by the numerical inversion of the given
% characteristic function CF by using the FFT algorithm.
%
% SYNTAX:
%  result = cf2DistFFT(cf,x,prob,options)
%  [result,cdf,pdf,qf] = cf2DistFFT(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle for the characteristic function CF,
%  x       - vector of values from domain of the distribution F, if x = [],
%            cf2DistFFT automatically selects vector x from the domain.
%  prob    - vector of values from [0,1] for which the quantiles will be
%            estimated, if prob = [], cf2DistFFT automatically selects
%            vector prob = [0.9,0.95,0.975,0.99,0.995,0.999].
%  options - structure with the following default parameters:   
%            options.isCompound = false  % treat the compound distributions
%                                        % of the RV Y = X_1 + ... + X_N,
%                                        % where N is discrete RV and X>=0 
%                                        % are iid RVs from nonnegative
%                                        % continuous distribution.
%             options.isInterp   = false % create and use the interpolant
%                                          functions for PDF/CDF/QF/RND
%             options.N = 2^10         % N points used by FFT
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
%                 options.DIST is created automatically after first call: 
%             options.DIST.xMin  = xMin   % the lower limit of X
%             options.DIST.xMax  = xMax   % the upper limit of X
%             options.DIST.xMean = xMean  % the MEAN value of X, 
%             options.DIST.cft   = cft    % CF evaluated at t_j : cf(t_j).   
%
% OUTPUTS:
%  result  - structure with the following results values:  
%             result.x                  = x;
%             result.cdf                = cdf;
%             result.pdf                = pdf;
%             result.prob               = prob;
%             result.qf                 = qf;
%             result.xFTT               = xFFT;
%             result.pdfFFT             = pdfFFT;
%             result.cdfFFT             = cdfFFT;
%             result.SixSigmaRule       = options.SixSigmaRule;
%             result.N                  = N;
%             result.dt                 = dt;
%             result.T                  = t(end);
%             result.PrecisionCrit      = PrecisionCrit;
%             result.myPrecisionCrit    = options.crit;
%             result.isPrecisionOK      = isPrecisionOK;
%             result.xMean              = xMean;
%             result.xStd               = xStd;
%             result.xMin               = xMin;
%             result.xMax               = xMax;
%             result.cf                 = cf;
%             result.options            = options;
%             result.tictoc             = toc;
%
% EXAMPLE 1:
%  % DISTRIBUTION OF A LINEAR COMBINATION OF THE INDEPENDENT RVs
%  % (Normal, Student's t, Rectangular, Triangular & Arcsine distribution)
%  % Y = X_{N} + X_{t} + 5*X_{R} + X_{T} + 10*X_{U}
%  % CFs: Normal, Student's t, Rectangular, Triangular, and Arcsine
%  cf_N  = @(t) exp(-t.^2/2);                                      
%  cf_t  = @(t,nu) min(1,besselk(nu/2, abs(t).*sqrt(nu),1) .* ...
%          exp(-abs(t).*sqrt(nu)) .* (sqrt(nu).*abs(t)).^(nu/2) / ...
%          2^(nu/2-1)/gamma(nu/2));     
%  cf_R  = @(t) min(1,sin(t)./t);   
%  cf_T  = @(t) min(1,(2-2*cos(t))./t.^2);
%  cf_U  = @(t) besselj(0,t);  
%  % Characteristic function of the linear combination Y
%  c    = [1 1 5 1 10]; nu = 1;
%  cf_Y = @(t) cf_N(c(1)*t) .* cf_t(c(2)*t,nu) .* cf_R(c(3)*t) .* ...
%         cf_T(c(4)*t) .* cf_U(c(5)*t);
%  clear options
%  options.N    = 2^10;
%  options.xMin = -50;
%  options.xMax = 50;
%  result = cf2DistFFT(cf_Y,[],[],options);
%  title('CDF of Y = X_{N} + X_{t} + 5*X_{R} + X_{T} + 10*X_{U}')
%
% EXAMPLE 2:
%  % DISTRIBUTION OF A LINEAR COMBINATION OF THE INDEPENDENT CHI2 RVs
%  % (Chi-squared RVs with 1 and 10 degrees of freedom)
%  % Y = 10*X_{Chi2_1} + X_{Chi2_10}
%  % Characteristic functions of X_{Chi2_1} and X_{Chi2_10}
%  %
%  df1       = 1;
%  df2       = 10;
%  cfChi2_1  = @(t) (1-2i*t).^(-df1/2);
%  cfChi2_10 = @(t) (1-2i*t).^(-df2/2);
%  cf_Y      = @(t) cfChi2_1(10*t) .* cfChi2_10(t);
%  clear options
%  options.xMin = 0;
%  result = cf2DistFFT(cf_Y,[],[],options);
%  title('CDF of Y = 10*X_{\chi^2_{1}} + X_{\chi^2_{10}}')
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Exponential distribution)
%  lambda1 = 10;
%  lambda2 = 5;
%  cfX  = @(t) cfX_Exponential(t,lambda2);
%  cf   = @(t) cfN_Poisson(t,lambda1,cfX);
%  x    = linspace(0,8,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = 1;
%  result = cf2DistFFT(cf,x,prob,options)
%
% REMARKS:
%  The outputs of the algorithm cf2DistFFT are approximate values! The
%  precission of the presented results depends on several different
%  factors: 
%  - application of the FFT algorithm for numerical inversion of the CF
%  - selected number of points used by the FFT algorithm (by default
%    options.N = 2^10),
%  - estimated/calculated domain [A,B] of the distribution F, used with the
%    FFT algorithm. Optimally, [A,B] covers large part of the
%    distribution domain, say more than 99%. However, the default
%    automatic procedure for selection of the domain [A,B] may fail. It
%    is based on the 'SixSigmaRule': A = MEAN - SixSigmaRule * STD, and
%    B = MEAN + SixSigmaRule * STD. Alternatively, change the
%    options.SixSigmaRule to different value, say 12, or use the
%    options.xMin and options.xMax to set manually the values of A and B.
%
% REFERENCES:
% [1] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [2] WITKOVSKY, V.: Matlab algorithm TDIST: The distribution of a
%     linear combination of Student’s t random variables. In COMPSTAT
%     2004 Symposium (2004), J. Antoch, Ed., Physica-Verlag/Springer
%     2004, Heidelberg, Germany, pp. 1995–2002.
% [3] WITKOVSKY, V.: WIMMER,G., DUBY, T. Logarithmic Lambert W x F
%     random variables for the family of chi-squared distributions
%     and their applications. Statistics & Probability Letters 96
%     (2015), 223–231.  
% [4] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
%
% SEE ALSO: cf2Dist, cf2DistGP, cf2DistGPT, cf2DistGPA, cf2DistFFT,
%           cf2DistBV, cf2CDF, cf2PDF, cf2QF 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 13:25:21
%
% Revision history:
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%[result,cdf,pdf,qf] = cf2DistFFT(cf,x,prob,options);

%% CHECK THE INPUT PARAMETERS
tic;
narginchk(1, 4);

if nargin < 4, options = []; end
if nargin < 3, prob = []; end
if nargin < 2, x = []; end

if ~isfield(options, 'isCompound')
    options.isCompound = false;
end


if ~isfield(options, 'N')
    if options.isCompound
        options.N = 2^14;
    else
        options.N = 2^10;
    end
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

if ~isfield(options, 'DIST')
    options.DIST = [];
end

% Other options parameters
if ~isfield(options, 'isPlotFFT')
    options.isPlotFFT = false;
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
%% GET/SET the DEFAULT parameters and the OPTIONS
% First, set a special treatment if the real value of CF at infinity (large
% value) is positive, i.e. const = real(cf(Inf)) > 0. In this case the
% compound CDF has jump at 0 of size equal to this value, i.e. cdf(0) =
% const, and pdf(0) = Inf. In order to simplify the calculations, here we
% calculate PDF and CDF of a distribution given by transformed CF, i.e.
% cf_new(t) = (cf(t)-const) / (1-const); which is converging to 0 at Inf,
% i.e. cf_new(Inf) = 0. Using the transformed CF requires subsequent
% recalculation of the computed CDF and PDF, in order to get the originaly
% required values: Set pdf_original(0) =  Inf & pdf_original(x) =
% pdf_new(x) * (1-const), for x > 0. Set cdf_original(x) =  const +
% cdf_new(x) * (1-const).
% 
const = real(cf(1e30));
if options.isCompound
    cfOld = cf;
    if const > 1e-13
        cf    = @(t) (cf(t) - const) / (1-const);
    end
end

if ~isempty(options.DIST)
    xMean              = options.DIST.xMean;
    cft                = options.DIST.cft;
    xMin               = options.DIST.xMin;
    xMax               = options.DIST.xMax;
    N                  = length(cft);
    k                  = (0:(N-1))';
    xRange              = xMax - xMin;
    dt                 = 2*pi / xRange;
    t                  = (k - N/2 + 0.5) * dt;
    xStd               = [];
else
    N                  = 2*options.N;
    dt                 = options.dt;
    T                  = options.T;
    xMin               = options.xMin;
    xMax               = options.xMax;
    xMean              = options.xMean;
    xStd               = options.xStd;
    SixSigmaRule       = options.SixSigmaRule;
    tolDiff            = options.tolDiff;
    cft                = cf(tolDiff*(1:4));
    if isempty(xMean)
        xMean = real((-cft(2) ...
            + 8*cft(1)-8*conj(cft(1)) ...
            + conj(cft(2)))/(1i*12*tolDiff));
    end
    if isempty(xStd)
        xM2 = real(-(conj(cft(4)) ...
            - 16*conj(cft(3)) ...
            + 64*conj(cft(2)) ...
            + 16*conj(cft(1)) ...
            - 130 + 16*cft(1) ...
            + 64*cft(2) ...
            - 16*cft(3)+cft(4))/(144*tolDiff^2));
        xStd  = sqrt(xM2 - xMean^2);
    end
    if (isfinite(xMin) && isfinite(xMax))
        xRange          = xMax - xMin;
    elseif ~isempty(T)
        xRange = 2*pi / (2 * T / N);
        if isfinite(xMax)
            xMin = xMax - xRange;
        elseif isfinite(xMin)
            xMax = xMin + xRange;
        else
            xMin = xMean - xRange/2;
            xMax = xMean + xRange/2;
        end
    elseif ~isempty(dt)
        xRange = 2*pi / dt;
        if isfinite(xMax)
            xMin = xMax - xRange;
        elseif isfinite(xMin)
            xMax = xMin + xRange;
        else
            xMin = xMean - xRange/2;
            xMax = xMean + xRange/2;
        end
    else  
        if isfinite(xMin)
            xMax       = xMean + SixSigmaRule * xStd;
        elseif isfinite(xMax)
            xMin       = xMean - SixSigmaRule * xStd;
        else
            xMin       = xMean - SixSigmaRule * xStd;
            xMax       = xMean + SixSigmaRule * xStd;
        end
        xRange          = xMax - xMin;        
    end
    dt                 = 2*pi / xRange;
    k                  = (0:(N-1))';
    t                  = (k - N/2 + 0.5) * dt;
    cft                = cf(t(N/2+1:end));
    cft                = [conj(cft(end:-1:1));cft];
    %cft(1)             = cft(1)/2;
    %cft(N)             = cft(N)/2;
    options.DIST.xMin  = xMin;
    options.DIST.xMax  = xMax;
    options.DIST.xMean = xMean;
    options.DIST.cft   = cft;
end

%% ALGORITHM

A      = xMin;
B      = xMax;
dx     = (B-A)/N;
c      = (-1).^(A*(N-1)/(B-A))/(B-A);
C      = c * (-1).^((1-1/N)*k);
D      = (-1).^(-2*(A/(B-A))*k);
pdfFFT = max(0,real(C.*fft(D.*cft)));
cdfFFT = min(1,max(0,0.5 + real(1i*C.*fft(D.*cft./t))));
xFFT   = A + k * dx;

% Reset the transformed CF, PDF, and CDF to the original values
if options.isCompound
    cf  = cfOld;
    cdfFFT = const + cdfFFT * (1-const);
    pdfFFT = pdfFFT * (1-const);
    pdfFFT(x==0) = inf;
end

% Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol, 
% PrecisionCrit should be small for t > T, smaller than tolerance
% options.crit 
PrecisionCrit = abs(cft(end)/t(end));
isPrecisionOK = (PrecisionCrit<=options.crit);

%% INTERPOLATE QUANTILE FUNCTION for required prob values: QF(prob)
if isempty(prob)
    prob = [0.9,0.95,0.975,0.99,0.995,0.999];
end
[cdfU,id] = unique(cdfFFT);
xxU   = xFFT(id);
szp   = size(prob);
qfFun = @(prob) interp1([-eps;cdfU],[-eps;xxU+dx/2],prob,'pchip');
qf    = reshape(qfFun(prob),szp);

% INTERPOLATE CDF required values: CDF(x)
% Default values if x = [];
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

% WARNING: Out-of-range
if any(x < xMin) || any(x > xMax)
    warning('VW:CharFunTool:cf2DistFFT',['x out-of-range: ', ...
        '[xMin, xMax] = [',num2str(xMin),...
        ', ',num2str(xMax),'] !']);
end

szx    = size(x);
%cdfFun = @(x) interp1([-eps;xxU+dx/2],[-eps;cdfU],x(:));
cdfFun = @(x) interp1(xxU,cdfU,x(:));
cdf    = reshape(cdfFun(x),szx);

% TRY INTERPOLATE PDF required values: PDF(x)
try
    pdfFun = @(x) interp1(xFFT,pdfFFT,x(:));
    pdf    = reshape(max(0,pdfFun(x)),szx);
catch
    warning('cf2DistFFT: Unable to interpolate the required PDF values')
    pdf = NaN*x;
end
x       = reshape(x,szx);

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
        warning('VW:CharFunTool:cf2DistGPT', ...
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
result.inversionMethod     = 'Discrete version of standard inversion formula';
result.quadratureMethod    = 'Fast Fourier Transform (FFT) algorithm';
result.x                  = x;
result.cdf                = cdf;
result.pdf                = pdf;
result.prob               = prob;
result.qf                 = qf;
if options.isInterp
    result.PDF             = PDF;
    result.CDF             = CDF;
    result.QF              = QF;
    result.RND             = RND;
end
result.xFTT               = xFFT;
result.pdfFFT             = pdfFFT;
result.cdfFFT             = cdfFFT;
result.SixSigmaRule       = options.SixSigmaRule;
result.N                  = options.N;
result.dt                 = dt;
result.T                  = t(end);
%result.PrecisionCrit      = PrecisionCrit;
%result.myPrecisionCrit    = options.crit;
result.isPrecisionOK      = isPrecisionOK;
result.xMean              = xMean;
result.xStd               = xStd;
result.xMin               = xMin;
result.xMax               = xMax;
result.cf                 = cf;
result.const              = const;
result.isCompound         = options.isCompound;
result.options            = options;
result.tictoc             = toc;

%% PLOT THE PDF/CDF, if required
if length(x)==1
    options.isPlot = false;
end
if options.isPlotFFT
    x = xFFT;
    pdf = pdfFFT;
    cdf = cdfFFT;
end
if options.isPlot
    figure
    plot(x,pdf,'-','LineWidth',1)
    grid
    title('PDF Specified by the Characteristic Function CF')
    xlabel('x')
    ylabel('pdf')
    
    figure
    plot(x,cdf,'-','LineWidth',1)
    grid
    title('CDF Specified by the Characteristic Function CF')
    xlabel('x')
    ylabel('cdf')
end

end