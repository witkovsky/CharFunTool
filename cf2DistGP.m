function [result,cdf,pdf,qf] = cf2DistGP(cf,x,prob,options)
%cf2DistGP Calculates the CDF/PDF/QF from the characteristic function CF
% by using the Gil-Pelaez inversion formulae:
%   cdf(x) = 1/2 + (1/pi) * Integral_0^inf imag(exp(-1i*t*x)*cf(t)/t)*dt.
%   pdf(x) = (1/pi) * Integral_0^inf real(exp(-1i*t*x)*cf(t))*dt.
%
% SYNTAX:
%  result = cf2DistGP(cf,x)
%  [result,cdf,pdf,qf] = cf2DistGP(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), x
%          - vector of x values where the CDF/PDF is computed, prob
%          - vector of values from [0,1] for which the quantiles
%            function is evaluated,
%  options - structure with the following default parameters:
%             options.isCompound = false  % treat the compound distributions
%                                         % of the RV Y = X_1 + ... + X_N,
%                                         % where N is discrete RV and X>=0 
%                                         % are iid RVs from nonnegative
%                                         % continuous distribution.
%             options.isCircular = false  % treat the circular
%                                         % distributions on (-pi,pi)
%             options.isInterp   = false  % create and use the interpolant
%                                           functions for PDF/CDF/QF/RND
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
%             options.DIST is created automatically after first call:
%             options.DIST.xMin  = xMin   % the lower limit of X
%             options.DIST.xMax  = xMax   % the upper limit of X
%             options.DIST.xMean = xMean  % the MEAN value of X,
%             options.DIST.cft   = cft    % CF evaluated at t_j : cf(t_j).
%
% REMARKS:
% If options.DIST is provided, then cf2DistGP evaluates CDF/PDF based on
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
%  qf       - vector of QF values evaluated at prob.
%
% REMARKS:
% The required integrals are evaluated approximately by using the simple
% trapezoidal rule on the interval(0,T), where T = N * dt is a sufficienly
% large integration upper limit in the frequency domain.
%
% If the optimum values of N and T are unknown, we suggest, as a simple
% rule of thumb, to start with the application of the six-sigma-rule for
% determining the value of dt = (2*pi)/(xMax-xMin), where xMax = xMean +
% 6*xStd, and xMin = xMean - 6*xStd, see [1].
%
% Please note that THIS (TRAPEZOIDAL) METHOD IS AN APPROXIMATE METHOD:
% Frequently, with relatively low numerical precision of the results of the
% calculated PDF/CDF/QF, but frequently more efficient and more precise
% than comparable Monte Carlo methods.
%
% However, the numerical error (truncation error and/or the integration
% error) could be and should be properly controled!
%
% CONTROLING THE PRECISION:
% Simple criterion for controling numerical precision is as follows: Set N
% and T = N*dt such that the value of the integrand function
% Imag(e^(-1i*t*x) * cf(t)/t) is sufficiently small for all t > T, i.e.
%   PrecisionCrit = abs(cf(t)/t) <= tol,
% for pre-selected small tolerance value, say tol = 10^-8. If this
% criterion is not valid, the numerical precission of the result is
% violated, and the method should be improved (e.g. by selecting larger N
% or considering other more sofisticated algorithm - not considered here).
% For more details consult the references below.
%
% EXAMPLE1 (Calculate CDF/PDF of N(0,1) by inverting its CF)
%  cf = @(t) exp(-t.^2/2);
%  result = cf2DistGP(cf)
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
%  result = cf2DistGP(cf,x,prob,options)
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
%  result = cf2DistGP(cf,x,prob,options)
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
% [5] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Preprint submitted
%     to Insurance: Mathematics and Economics.
% [6] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models. Preprint submitted
%     to Journal of Statistical Software.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%[result,cdf,pdf,qf] = cf2DistGP(cf,x,prob,options);

%% CHECK THE INPUT PARAMETERS
tic;
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
if ~isfield(options, 'qf0')
    options.qf0 = (cf(1e-4)-cf(-1e-4))/(2e-4*1i);
end

if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end

if ~isfield(options, 'xN')
    options.xN = 101;
end

if ~isfield(options, 'chebyPts')
    options.chebyPts = 2^8;
end

if ~isfield(options, 'correctedCDF')
    options.correctedCDF = false;
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
        %    prob  = max(0,(prob - const) / (1-const));
    end
end

if ~isempty(options.DIST)
    xMean              = options.DIST.xMean;
    cft                = options.DIST.cft;
    xMin               = options.DIST.xMin;
    xMax               = options.DIST.xMax;
    range              = xMax - xMin;
    dt                 = 2*pi / range;
    N                  = length(cft);
    t                  = (1:N)' * dt;
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
    cft                = cf(tolDiff*(1:4));
    cftRe              = real(cft);
    cftIm              = imag(cft);
    if isempty(xMean)
        if options.isCircular
            % see https://en.wikipedia.org/wiki/Directional_statistics
            xMean = angle(cf(1));
        else
            %             xMean = real((-cft(2) ...
            %                 + 8*cft(1)-8*conj(cft(1)) ...
            %                 + conj(cft(2)))/(1i*12*tolDiff));
            % Alternative
            %             xMean = real(((1/280)*conj(cft(4)) ...
            %                 - (4/105)*conj(cft(3)) ...
            %                 + (1/5)*conj(cft(2)) ...
            %                 - (4/5)*conj(cft(1)) ...
            %                 + (4/5)*cft(1) ...
            %                 - (1/5)*cft(2) ...
            %                 + (4/105)*cft(3) ...
            %                 - (1/280)*cft(4))/(1i*tolDiff));
            xMean = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 ...
                - 2*cftIm(4)/280) / tolDiff;
        end
    end
    if isempty(xStd)
        if options.isCircular
            % see https://en.wikipedia.org/wiki/Directional_statistics
            xStd  = sqrt(-2*log(abs(cf(1))));
        else
            %             xM2 = real(-(conj(cft(4)) ...
            %                 - 16*conj(cft(3)) ...
            %                 + 64*conj(cft(2)) ...
            %                 + 16*conj(cft(1)) ...
            %                 - 130 + 16*cft(1) ...
            %                 + 64*cft(2) ...
            %                 - 16*cft(3)+cft(4))/(144*tolDiff^2));
            % Alternative
            %             xM2 = real(-(-(1/560)*conj(cft(4)) ...
            %                 + (8/315)*conj(cft(3)) ...
            %                 - (1/5)*conj(cft(2)) ...
            %                 + (8/5)*conj(cft(1)) ...
            %                 - (205/72) ...
            %                 + (8/5)*cft(1) ...
            %                 - (1/5)*cft(2) ...
            %                 + (8/315)*cft(3) ...
            %                 - (1/560)*cft(4))/tolDiff^2);
            xM2   = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 ...
                - 16*cftRe(3)/315 + 2*cftRe(4)/560) / tolDiff^2;
            xStd  = sqrt(xM2 - xMean^2);
        end
    end
    if (isfinite(xMin) && isfinite(xMax))
        range          = xMax - xMin;
    elseif ~isempty(dt)
        range          = 2*pi / dt;
    elseif ~isempty(T)
        range          = 2*pi / (T / N);
    else
        if options.isCircular
            xMin       = -pi;
            xMax       = pi;
        else
            if isfinite(xMin)
                xMax       = xMean + SixSigmaRule * xStd;
            elseif isfinite(xMax)
                xMin       = xMean - SixSigmaRule * xStd;
            else
                xMin       = xMean - SixSigmaRule * xStd;
                xMax       = xMean + SixSigmaRule * xStd;
            end
        end
        range          = xMax - xMin;
    end
    dt                 = 2*pi / range;
    t                  = (1:N)' * dt;
    cft                = cf(t);
    cft(N)             = cft(N)/2;
    options.DIST.xMin  = xMin;
    options.DIST.xMax  = xMax;
    options.DIST.xMean = xMean;
    options.DIST.cft   = cft;
end

%% ALGORITHM
% Default values if x = [];
if isempty(x)
    x = linspace(xMax,xMin,options.xN);
end

if options.isInterp
    xOrg = x;
    % Chebyshev points
    x = (xMax-xMin) * (-cos(pi*(0:options.chebyPts) / ...
        options.chebyPts) + 1) / 2 + xMin;
else
    xOrg = [];
end

% WARNING: Out-of-range
if any(x < xMin) || any(x > xMax)
    warning('VW:CFTOOLBOX:cf2DistGP',['x out-of-range (the used support): ', ...
        '[xMin, xMax] = [',num2str(xMin),...
        ', ',num2str(xMax),'] !']);
end

% Evaluate the required functions
[n,m]   = size(x);
x       = x(:);
E       = exp(-1i*x*t');

% CDF estimate computed by using the simple trapezoidal
% quadrature rule
cdf     = (xMean - x)/2 + imag(E * (cft ./ t));
cdf     = 0.5 - (cdf * dt) / pi;

% Correct the CDF (if the computed result is out of (0,1))
% This is useful for circular distributions over intervals of length 2*pi,
% as e.g. the von Mises distribution
corrCDF = 0;
if options.correctedCDF
    if min(cdf) < 0
        corrCDF = min(cdf);
        cdf     = cdf - corrCDF;
    end
    if max(cdf) > 1
        corrCDF = max(cdf)-1;
        cdf     = cdf - corrCDF;
    end
end
cdf     = reshape(max(0,min(1,cdf)),n,m);

% PDF estimate computed by using the simple trapezoidal
% quadrature rule
pdf     = 0.5 + real(E * cft);
pdf     = (pdf * dt) / pi;
pdf     = reshape(max(0,pdf),n,m);
x       = reshape(x,n,m);

% REMARK:
% Note that, exp(-1i*x_i*0) = cos(x_i*0) + 1i*sin(x_i*0) = 1. Moreover,
% cf(0) = 1 and lim_{t -> 0} cf(t)/t = E(X) - x. Hence, the leading term of
% the trapezoidal rule for computing the CDF integral is CDFfun_1 = (xMean
% - x)/2, and PDFfun_1 = 1/2 for the PDF integral, respectively.

% Reset the transformed CF, PDF, and CDF to the original values
if options.isCompound
    cf  = cfOld;
    cdf = const + cdf * (1-const);
    pdf = pdf * (1-const);
    pdf(x==0) = inf;
    pdf(x==xMax) = NaN;
end

% Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol,
% PrecisionCrit should be small for t > T, smaller than tolerance
% options.crit
PrecisionCrit = abs(cft(end)/t(end));
isPrecisionOK = (PrecisionCrit<=options.crit);

%% QF evaluated by the Newton-Raphson iterative scheme
if ~isempty(prob)
    isPlot = options.isPlot;
    options.isPlot = false;
    [n,m]     = size(prob);
    prob      = prob(:);
    maxiter   = options.maxiter;
    crit      = options.crit;
    qf        = options.qf0;
    criterion = true;
    count     = 0;
    [res,cdfQ,pdfQ] = cf2DistGP(cf,qf,[],options);
    options = res.options;
    while criterion
        count  = count + 1;
        correction  = ((cdfQ-corrCDF) - prob) ./ pdfQ;
        qf = max(xMin,min(xMax,qf - correction));
        [~,cdfQ,pdfQ] = cf2DistGP(cf,qf,[],options);
        criterion = any(abs(correction) ...
            > crit * abs(qf)) ...
            && max(abs(correction)) ...
            > crit && count < maxiter;
    end
    qf   = reshape(qf,n,m);
    prob = reshape(prob,n,m);
    options.isPlot = isPlot;
else
    qf = [];
    count = [];
    correction =[];
end

if options.isInterp
    try
        id   = isfinite(pdf);
        pdfFunction = @(xNew) max(0, ...
            interpolateBarycentric(x(id),pdf(id),xNew));
        PDF  = @(x) pdfFunction(x);
        id   = isfinite(cdf);
        cdfFunction = @(xNew) max(0,min(1, ...
            interpolateBarycentric(x(id),cdf(id),xNew)));
        CDF  = @(x) cdfFunction(x);
        qfFunction = @(prob) interpolateBarycentric(cdf(id),x(id),prob);
        QF   = @(prob) qfFunction(prob);
        rndFunction = @(n) qfFunction(rand(n,1));
        RND  = @(n) rndFunction(n);
    catch
        warning('VW:CFTOOLBOX:cf2DistGP', ...
            'Problem using the interpolant function');
        PDF  = [];
        CDF  = [];
        QF   = [];
        RND  = [];
    end
end

if ~isempty(xOrg)
    x   = xOrg;
    cdf = CDF(x);
    pdf = PDF(x);
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
result.SixSigmaRule       = options.SixSigmaRule;
result.N                  = N;
result.dt                 = dt;
result.T                  = t(end);
result.PrecisionCrit      = PrecisionCrit;
result.myPrecisionCrit    = options.crit;
result.isPrecisionOK      = isPrecisionOK;
result.xMean              = xMean;
result.xStd               = xStd;
result.xMin               = xMin;
result.xMax               = xMax;
result.cf                 = cf;
result.const              = const;
result.isCompound         = options.isCompound;
result.isCircular         = options.isCircular;
result.isInterp           = options.isInterp;
result.details.count      = count;
result.details.correction = correction;
result.details.corrCDF    = corrCDF;
result.options            = options;
result.tictoc             = toc;

%% PLOT the PDF / CDF / QF
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
    
    %     % QF
    %     if ~isempty(QF)
    %         figure
    %         prob = (-cos(pi*(0:options.xN)/options.xN)+1)/2;
    %         plot(prob,QF(prob),'Linewidth',1.6)
    %         grid
    %         title('QF Specified by the CF')
    %         xlabel('prob')
    %         ylabel('qf')
    %     end
    
end

end
%% Function interpolateBarycentric
function [funNew,xNew] = interpolateBarycentric(x,fun,xNew,options)
%   interpolateBarycentric Evaluates (interpolates) function values f_new at
%   given points x_new by barycentric interpolation from function values f
%   given at chebpoints x.
%   For more details see:
%   https://en.wikipedia.org/wiki/Barycentric_coordinate_system
%
% SYNTAX:
%   funNew = interpolateBarycentric(x,fun,xNew)
%   [funNew,xNew] = interpolateBarycentric(x,fun,xNew,options)
%
% EXAMPLE1 (Barycentric interpolant of the Sine function on (-pi,pi))
%   x    = chebpts(32,[-pi,pi])';
%   f    = sin(x);
%   xNew = linspace(-pi,pi,201)';
%   fNew = interpolateBarycentric(x,f,xNew,options);
%   plot(x,f,xNew,fNew,'.')
%   disp([xNew fNew sin(xNew)])
%
% EXAMPLE2 (Barycentric interpolant of the Normal CDF)
%   x    = linspace(-8,8,101)';
%   f    = normcdf(x);
%   xNew = linspace(-5,5,201)';
%   fNew = interpolateBarycentric(x,f,xNew);
%   plot(x,f,xNew,fNew,'.')
%   disp([xNew fNew normcdf(xNew)])
%
% EXAMPLE3 (Barycentric interpolant of the Normal quantile function)
%   x    = linspace(-8,8,2^7)';
%   cdf  = normcdf(x);
%   prob = linspace(1e-4,1-1e-4,101)';
%   qf   = interpolateBarycentric(cdf,x,prob);
%   plot(cdf,x,prob,qf,'.')
%   disp([prob qf norminv(prob)])
%
% EXAMPLE4 (Barycentric interpolant of the Compound distribution)
%   cfX  = @(t) cfX_Exponential(t,5);
%   cf   = @(t) cfN_Poisson(t,10,cfX);
%   x    = chebpts(2^9,[0,10]);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.isCompound = 1;
%   result = cf2DistGP(cf,x,prob,options)
%   CDF  = @(x) interpolateBarycentric(result.x,result.cdf,x);
%   PDF  = @(x) interpolateBarycentric(result.x,result.pdf,x);
%   QF   = @(p) interpolateBarycentric(result.cdf,result.x,p);
%   prob = linspace(1e-4,1-1e-4,201);
%   figure;plot(result.cdf,result.x,prob,QF(prob),'.');

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Jan-2017 17:28:22

%% ALGORITHM
%cf = cfX_LogGamma(t,alpha,beta,tol);

%% CHECK THE INPUT PARAMETERS

if nargin < 2
    error(message('VW:interpolateBarycentric:TooFewInputs'));
end

if nargin < 3, xNew = []; end
if nargin < 4, options = []; end

if isempty(xNew)
    xNew = linspace(min(x),max(x))';
end

if ~isfield(options, 'isChebPts')
    options.isChebPts = true;
end

%% Set the parameters

x      = x(:);
fun    = fun(:);
[n,m]  = size(xNew);
xNew   = xNew(:);
nx     = length(x);
nxNew  = length(xNew);
funNew = zeros(nxNew,1);

w     = (-1).^(0:nx-1).';
w(1)  = w(1)/2;
w(nx) = w(nx)/2;
for i = 1:nxNew
    A = 0;
    B = 0;
    for j = 1:nx
        if xNew(i) == x(j)
            exactVal  = true;
            funNew(i) = fun(j);
        else
            exactVal = false;
            weight   = w(j) ./ (xNew(i)-x(j));
            A = A + weight * fun(j);
            B = B + weight;
        end
        if exactVal == true
            break;
        else
            funNew(i) = A / B;
        end
    end
end

xNew   = reshape(xNew,n,m);
funNew = reshape(funNew,n,m);
end