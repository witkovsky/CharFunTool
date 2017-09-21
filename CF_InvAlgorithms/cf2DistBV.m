function [result,cdf,pdf,qf] = cf2DistBV(cf,x,prob,options)
%cf2DistBV Calculates the CDF/PDF/QF from the characteristic function CF
%  by using the Gil-Pelaez inversion formulae and the BAKHVALOV-VASILEVA
%  method:  
%    cdf(x) = 1/2 + (1/pi) * Integral_0^inf imag(exp(-1i*t*x)*cf(t)/t)*dt.
%    pdf(x) = (1/pi) * Integral_0^inf real(exp(-1i*t*x)*cf(t))*dt.
% 
%  The FOURIER INTEGRALs are calculated by using the BAKHVALOV-VASILEVA
%  method suggested for computing the oscillatory Fourier integrals based
%  on approximation of the integrand function by the FOURIER-LEGENDRE
%  SERIES EXPANSION, and observation that Fourier transform of the Legendre
%  polynomials is related to the Belssel J functions. For more details see
%  Bakhlanov and Vasileva (1968) and Evans and Webster (1999).
%
% WORK UNDER DEVELOPMENT:
%  CURRENTLY, cf2DistBV WORKS ONLY FOR NON-NEGADIVE DISTRIBUTIONS, i.e with
%  support x >= 0, and setting options.Nonnegative = 'true' ! 
%
% SYNTAX:
%  result = cf2DistBV(cf,x)
%  [result,cdf,pdf,qf] = cf2DistBV(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), 
%  x       - vector of x values where the CDF/PDF is computed, 
%  prob    - vector of values from [0,1] for which the quantiles
%            function is evaluated,
%  options - structure with the following default parameters:
%             options.isNonnegative = true % Gil-Pelaez formulae for
%                                          % non-negative distributions
%             options.isCompound = false % treat the compound distributions
%                                        % of the RV Y = X_1 + ... + X_N,
%                                        % where N is discrete RV and X>=0 
%                                        % are iid RVs from nonnegative
%                                        % continuous distribution.
%             options.isInterp   = false % create and use the interpolant
%                                          functions for PDF/CDF/QF/RND
%             options.xMin = -Inf or 0 % set the lower limit of X
%             options.xMax = Inf       % set the lower limit of X
%             options.SixSigmaRule = 6 % set the rule for computing domain
%             options.tolDiff = 1e-4   % tol for numerical differentiation
%             options.isPlot = true    % plot the graphs of PDF/CDF
%             options.tolCoefs = 1e-12 % tol for Legendre coefficients
%             options.nMax = 100       % nMax of Legendre coefficients
%             options.nLimits = 5      % nLimits+1 integration subintervals
%             options.Limits           %  = [0 10^(-1:nLimits) Inf]
%             options.DIST             % structure with precomputed Info
%             options.qf0              % starting quantile for iterations
%             options.crit = 1e-12;    % stopping criterion for quantiles 
%             options.maxiter = 1000   % max numner of quantile iterations
%             options.xN = 101         % length of the default x vector 
%             options.chebyPts = 2^9   % number of Chebyshev points used
%                                      % for Barycentric Interpolation
%
% OUTPUT:
%  result   - structure with CDF/PDF/QF and further details,
%  cdf      - vector of CDF values evaluated at x,
%  pdf      - vector of PDF values evaluated at x,
%  qf       - vector of QF values evaluated at prob.
%
% EXAMPLE1 (Calculate CDF/PDF of Exponential(1) by inverting its CF)
%  lambda = 1;
%  cf = @(t) cfX_Exponential(t,lambda);
%  result = cf2DistBV(cf)
%
% EXAMPLE2 (PDF/CDF of the compound Binomial-Exponential distribution)
%  n = 25;
%  p = 0.3;
%  lambda = 5;
%  cfX  = @(t) cfX_Exponential(t,lambda);
%  cf   = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,5,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = 1;
%  result = cf2DistBV(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Exponential distribution)
%  lambda1 = 10;
%  lambda2 = 5;
%  cfX  = @(t) cfX_Exponential(t,lambda2);
%  cf   = @(t) cfN_Poisson(t,lambda1,cfX);
%  x    = linspace(0,8,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  options.isInterp   = true;
%  result = cf2DistBV(cf,x,prob,options)
%
% REFERENCES:
% [1] BAKHVALOV, N.S., VASILEVA, L.G. Evaluation of the integrals of
%     oscillating functions by interpolation at nodes of Gaussian
%     quadratures. USSR Computational Mathematics and Mathematical Physics,
%     1968, 8(1): 241-249.
% [2] EVANS, G.A.; WEBSTER, J.R. A comparison of some methods for the
%     evaluation of highly oscillatory integrals. Journal of Computational
%     and Applied Mathematics, 1999, 112(1): 55-69.
% [3] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference, 2001, 94,
%     1–13. 
% [4] WITKOVSKY, V.: Matlab algorithm TDIST: The distribution of a
%     linear combination of Student’s t random variables. In COMPSTAT
%     2004 Symposium (2004), J. Antoch, Ed., Physica-Verlag/Springer
%     2004, Heidelberg, Germany, pp. 1995–2002.
% [5] WITKOVSKY, V., WIMMER, G., DUBY, T. Logarithmic Lambert W x F
%     random variables for the family of chi-squared distributions
%     and their applications. Statistics & Probability Letters, 2015, 96,
%     223–231. 
% [6] WITKOVSKY, V.: Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 2016, 5(3),
%     32-44. 
% [7] WITKOVSKY, V., WIMMER, G., DUBY, T. Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. ArXiv preprint,
%     2017, arXiv:1701.08299.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 21-Sep-2017 09:33:32

%% ALGORITHM
%[result,cdf,pdf,qf] = cf2DistBV(cf,x,prob,options);

%% CHECK THE INPUT PARAMETERS
timeVal = tic;
narginchk(1, 4);

if nargin < 4, options = []; end
if nargin < 3, prob = []; end
if nargin < 2, x = []; end

if ~isfield(options, 'isNonnegative')
    options.isNonnegative = true;
end

if ~isfield(options, 'isCompound')
    options.isCompound = false;
end

if ~isfield(options, 'isInterp')
    options.isInterp = false;
end

if ~isfield(options, 'xMin')
    if options.isNonnegative
        options.xMin = 0;
    else
        options.xMin = -Inf;
    end
end

if ~isfield(options, 'xMax')
    options.xMax = Inf;
end

if ~isfield(options, 'SixSigmaRule')
    options.SixSigmaRule = 6;
end

if ~isfield(options, 'tolDiff')
    options.tolDiff = 1e-4;
end

if ~isfield(options, 'tolCoefs')
    options.tolCoefs = 1e-15;
end

if ~isfield(options, 'nLimits')
    options.nLimits = 5;
end

if ~isfield(options, 'Limits')
    options.Limits = [1e-15 10.^(-1:options.nLimits) 1e+300];
end

if ~isfield(options, 'nMax')
    options.nMax = 100;
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

if ~isfield(options, 'DIST')
    options.DIST = [];
end

if ~isfield(options, 'qf0')
    options.qf0 = real((cf(1e-4)-cf(-1e-4))/(2e-4*1i));
end

if ~isfield(options, 'crit')
    options.crit = 1e-10;
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

%% SET the DEFAULT parameters

% Special treatment for compound distributions. If the real value of CF at
% infinity (large value) is positive constant, i.e. const = real(cf(Inf)) >
% 0. In this case the compound CDF has jump at 0 of size equal to this
% value, i.e. cdf(0) = const, and pdf(0) = Inf. In order to simplify the
% calculations, here we calculate PDF and CDF of a distribution given by
% transformed CF, i.e. cf_new(t) = (cf(t)-const) / (1-const); which is
% converging to 0 at Inf, i.e. cf_new(Inf) = 0. Using the transformed CF
% requires subsequent recalculation of the computed CDF and PDF, in order
% to get the originaly required values: Set pdf_original(0) =  Inf &
% pdf_original(x) = pdf_new(x) * (1-const), for x > 0. Set cdf_original(x)
% =  const + cdf_new(x) * (1-const).

const = real(cf(1e300));
if options.isCompound
    cfOld = cf;
    if const > 1e-13
        cf    = @(t) (cf(t) - const) / (1-const);
    end
    options.isNonnegative = true;
end

if ~isempty(options.DIST)
    xMean        = options.DIST.xMean;
    xMin         = options.DIST.xMin;
    xMax         = options.DIST.xMax;
    xStd         = options.DIST.xStd;
    cdfCoef      = options.DIST.cdfCoef;
    pdfCoef      = options.DIST.pdfCoef;
    scale        = options.DIST.scale;
    shift        = options.DIST.shift;
    nLimits      = options.DIST.nLimits;
    nTrue        = options.DIST.nTrue;
else
    xMin         = options.xMin;
    xMax         = options.xMax;
    SixSigmaRule = options.SixSigmaRule;
    tolDiff      = options.tolDiff;
    cft          = cf(tolDiff*(1:4));
    cftRe        = real(cft);
    cftIm        = imag(cft);
    xMean        = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 ...
                   - 2*cftIm(4)/280) / tolDiff;
    xM2          = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 ...
                   - 16*cftRe(3)/315 + 2*cftRe(4)/560) / tolDiff^2;
    xStd         = sqrt(xM2 - xMean^2);
    if isfinite(xMin) && ~isfinite(xMax)
        xMax     = xMean + SixSigmaRule * xStd;
    elseif isfinite(xMax) && ~isfinite(xMin)
        xMin     = xMean - SixSigmaRule * xStd;
    elseif ~isfinite(xMax) && ~isfinite(xMin)
        xMin     = xMean - SixSigmaRule * xStd;
        xMax     = xMean + SixSigmaRule * xStd;
    end
    pdfFun   = @(t) real(cf(t));
    cdfFun   = @(t) imag(cf(t))./t;
    nMax     = options.nMax;
    tolCoefs = options.tolCoefs;
    Limits   = options.Limits;
    nLimits  = length(Limits);
    cdfCoef  = cell(1,nLimits-1);
    pdfCoef  = cell(1,nLimits-1);
    scale    = cell(1,nLimits-1);
    shift    = cell(1,nLimits-1);
    % Get the coefficients of Legendre series expansion of the required
    % integrand functions: cdfFun and pdfFun
    A = Limits(1);
    B = Limits(2);
    nTrue = 0;
    [cdfCoef{1},scale{1},shift{1},xx,ww,PP] = ...
        LegendreSeries(cdfFun,A,B,nMax,tolCoefs);
    nTrue = max(length(cdfCoef{1}),nTrue);
    pdfCoef{1} = LegendreSeries(pdfFun,A,B,nMax,tolCoefs,xx,ww,PP);
    nTrue = max(length(pdfCoef{1}),nTrue);
    for i = 1:(nLimits-2)
        A = B;
        B = Limits(i+2);
        [cdfCoef{i+1},scale{i+1},shift{i+1}] = ...
            LegendreSeries(cdfFun,A,B,nMax,tolCoefs,xx,ww,PP);
        nTrue = max(length(cdfCoef{i+1}),nTrue);
        pdfCoef{i+1} = LegendreSeries(pdfFun,A,B,nMax,tolCoefs,xx,ww,PP);
        nTrue = max(length(pdfCoef{i+1}),nTrue);
    end
    % Save options.DIST
    options.DIST.xMin    = xMin;
    options.DIST.xMax    = xMax;
    options.DIST.xMean   = xMean;
    options.DIST.xStd    = xStd;
    options.DIST.cdfCoef = cdfCoef;
    options.DIST.pdfCoef = pdfCoef;
    options.DIST.scale   = scale;
    options.DIST.shift   = shift;
    options.DIST.nLimits = nLimits;
    options.DIST.nTrue   = nTrue;
end

if isempty(x)
    x    = linspace(xMin,xMax,options.xN);
end

if options.isInterp
    x0   = x;
    xMin = min(min(x),xMin);
    xMax = max(max(x),xMax);
    x    = ChebPoints(options.chebyPts,[xMin,xMax]);
else
    x0   = [];
end

szx      = size(x);
x        = x(:);

%% ALGORITHM

% Integrate and sum-up over all subintervals
cdf = 0;
pdf = 0;
for i = 1:(nLimits-1)
    [I1,I2] = LegendreSeriesFourierIntegral(x,cdfCoef{i}, ...
        scale{i},shift{i},pdfCoef{i});  
    cdf = cdf + I1;
    pdf = pdf + I2;
end

% Use the Gil-Pelaez formulae for PDF and CDF
if options.isNonnegative
    pdf = 2 * real(pdf) / pi;
    cdf = 1 - 2 * real(cdf) / pi;
else
    warning('VW:CharFunTool:cf2DistBV', ...
        'NOW, cf2DistBV WORKS ONLY FOR NON-NEGADIVE DISTRIBUTIONS !');
    pdf = [];
    cdf = [];
end

% Reset the transformed CF, PDF, and CDF to the original values
if options.isCompound
    cf  = cfOld;
    cdf = const + cdf * (1-const);
    pdf = pdf * (1-const);
    pdf(x==0) = 0;
    pdf(x==xMax) = NaN;
end

% Reshape the results to the original size of x
x    = reshape(x,szx);
pdf  = reshape(max(0,pdf),szx);
cdf  = reshape(min(1,max(0,cdf)),szx);

%% Quantile function: QF evaluated by the Newton-Raphson iterative scheme
if ~isempty(prob)
    isPlot = options.isPlot;
    options.isPlot = false;
    isInterp = options.isInterp;
    options.isInterp = false;
    [n,m]     = size(prob);
    prob      = prob(:);
    maxiter   = options.maxiter;
    crit      = options.crit;
    qf        = options.qf0;
    criterion = true;
    count     = 0;
    [res,cdfQ,pdfQ] = cf2DistBV(cf,qf,[],options);
    options = res.options;
    while criterion
        count  = count + 1;
        correction  = (cdfQ - prob) ./ pdfQ;
        qf = max(xMin,min(xMax,qf - correction));
        [~,cdfQ,pdfQ] = cf2DistBV(cf,qf,[],options);
        criterion = any(abs(correction) ...
            > crit * abs(qf)) ...
            && max(abs(correction)) ...
            > crit && count < maxiter;
    end
    qf   = reshape(qf,n,m);
    prob = reshape(prob,n,m);
    options.isPlot = isPlot;
    options.isInterp = isInterp;
else
    qf = [];
    count = [];
    correction =[];
end

if options.isInterp
    try
        id   = isfinite(pdf);
        pdfFunction = @(xNew) max(0, ...
            InterpBarycentric(x(id),pdf(id),xNew));
        PDF  = @(x) pdfFunction(x);
        id   = isfinite(cdf);
        cdfFunction = @(xNew) max(0,min(1, ...
            InterpBarycentric(x(id),cdf(id),xNew)));
        CDF  = @(x) cdfFunction(x);
        qfFunction = @(prob) InterpBarycentric(cdf(id),x(id),prob);
        QF   = @(prob) qfFunction(prob);
        rndFunction = @(n) qfFunction(rand(n,1));
        RND  = @(n) rndFunction(n);
    catch
        warning('VW:CharFunTool:cf2DistBV', ...
            'Problem using the interpolant function');
        PDF  = [];
        CDF  = [];
        QF   = [];
        RND  = [];
    end
end

if ~isempty(x0)
    x   = x0;
    cdf = CDF(x);
    pdf = PDF(x);
end

% Reset the correct value for compound PDF at 0
if options.isCompound
    pdf(x==0) = Inf;
end

%% RESULT
result.Description        = 'CDF/PDF/QF from the characteristic function CF';
result.x                  = x;
result.cdf                = cdf;
result.pdf                = pdf;
result.prob               = prob;
result.qf                 = qf;
if options.isInterp
    result.PDF            = PDF;
    result.CDF            = CDF;
    result.QF             = QF;
    result.RND            = RND;
end
result.cf                 = cf;
result.isNonnegative      = options.isNonnegative;
result.isCompound         = options.isCompound;
result.isInterp           = options.isInterp;
result.SixSigmaRule       = options.SixSigmaRule;
result.xMean              = xMean;
result.xStd               = xStd;
result.xMin               = xMin;
result.xMax               = xMax;
result.const              = const;
result.details.count      = count;
result.details.correction = correction;
result.details.nTrue      = nTrue;
result.options            = options;
result.isOK               = count < options.maxiter & nTrue < options.nMax;
result.tictoc             = toc(timeVal);

%% PLOT the PDF / CDF
if length(x)==1
    options.isPlot = false;
end
if options.isPlot
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