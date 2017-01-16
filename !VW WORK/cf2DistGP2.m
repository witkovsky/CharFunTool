function [result,cdf,pdf,qf] = cf2DistGP2(cf,x,prob,options)
%cf2DistGP2 Calculates the CDF/PDF/QF from the characteristic function CF 
% by using the Gil-Pelaez inversion formulae: 
%   cdf(x) = 1/2 + (1/pi) * Integral_0^inf imag(exp(-1i*t*x)*cf(t)/t)*dt. 
%   pdf(x) = (1/pi) * Integral_0^inf real(exp(-1i*t*x)*cf(t))*dt.
%
% SYNTAX:
%  result = cf2DistGP2(cf,x)
%  [result,cdf,pdf,qf] = cf2DistGP2(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), x            
%          - vector of x values where the CDF/PDF is computed, prob         
%          - vector of values from [0,1] for which the quantiles
%            function is evaluated,
%  options - structure with the following default parameters:   
%             options.isCompound = false  % treat the compond distributions
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
% If options.DIST is provided, then cf2DistGP2 evaluates CDF/PDF based on
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
%  result = cf2DistGP2(cf)
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
%  result = cf2DistGP2(cf,x,prob,options)
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
%  result = cf2DistGP2(cf,x,prob,options)
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
%[result,cdf,pdf,qf] = cf2DistGP2(cf,x,prob,options);

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
if ~isfield(options, 'qf0')
    options.qf0 = (cf(1e-4)-cf(-1e-4))/(2e-4*1i);
end

if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
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
    if const > 1e-13
        cfOld = cf;
        cf    = @(t) (cf(t) - const) / (1-const);
    %    prob  = max(0,(prob - const) / (1-const));
    end
end

if ~isempty(options.DIST)
    xMean              = options.DIST.xMean;
    cft1                = options.DIST.cft1;
    cft2                = options.DIST.cft2;
    xMin               = options.DIST.xMin;
    xMax               = options.DIST.xMax;
    range              = xMax - xMin;
    N                  = 2*length(cft1);
    dt                 = 2*pi / range;
    dt1                = options.DIST.dt1;
    dt2                = options.DIST.dt2;
    t1                 = (1:N/2)' * dt1;
    t2                 = t1(end) + (1:N/2)' * dt2;
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
        range          = xMax - xMin;
    elseif ~isempty(dt)
        range          = 2*pi / dt;
    elseif ~isempty(T)
        range          = 2*pi / (T / N);
    else
        if isfinite(xMin)
            xMax       = xMean + SixSigmaRule * xStd;
        elseif isfinite(xMax)
            xMin       = xMean - SixSigmaRule * xStd;
        else
            xMin       = xMean - SixSigmaRule * xStd;
            xMax       = xMean + SixSigmaRule * xStd;
        end
        range          = xMax - xMin;        
    end
    dt                 = 2*pi / range;
    dt1                = dt/2;
    dt2                = 10*dt;
    %t                  = (1:N)' * dt;
    t1                 = (1:N/2)' * dt1;
    t2                 = t1(end) + (1:N/2)' * dt2;
    cft1               = cf(t1);
    cft1(N/2)          = cft1(N/2)/2;
    cft2               = cf(t2);
    cft2(1)            = cft2(1)/2;
    cft2(N/2)          = cft2(N/2)/2;
    options.DIST.xMin  = xMin;
    options.DIST.xMax  = xMax;
    options.DIST.xMean = xMean;
    options.DIST.cft1   = cft1;
    options.DIST.cft2   = cft2;
    options.DIST.dt1   = dt1;
    options.DIST.dt2   = dt2;
end

%% ALGORITHM
% Default values if x = [];
if isempty(x)
    x = linspace(xMin,xMax,101)';    
end

% WARNING: Out-of-range
if any(x < xMin) || any(x > xMax)
warning(['x out-of-range (the used support): ', ...
        '[xMin, xMax] = [',num2str(xMin),...
        ', ',num2str(xMax),'] !']);
end

% Evaluate the required functions
[n,m]   = size(x);
x       = x(:);
E1       = exp(-1i*x*t1');
E2       = exp(-1i*x*t2');

% CDF estimate computed by using the simple trapezoidal
% quadrature rule 
cdf     = ((xMean - x)/2 + imag(E1 * (cft1 ./ t1)))*dt1 + imag(E2 * (cft2 ./ t2))*dt2;
cdf     = 0.5 - cdf / pi;
cdf     = reshape(max(0,min(1,cdf)),n,m);

% PDF estimate computed by using the simple trapezoidal
% quadrature rule 
pdf     = (0.5 + real(E1 * cft1))*dt1 + real(E2 * cft2)*dt2;
pdf     = pdf / pi;
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
end

% Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol,
% PrecisionCrit should be small for t > T, smaller than tolerance
% options.crit 
PrecisionCrit = abs(cft2(end)/t2(end));
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
    [res,cdfQ,pdfQ] = cf2DistGP2(cf,qf,[],options);
    options = res.options;
    while criterion
        count  = count + 1;
        correction  = (cdfQ - prob) ./ pdfQ;
        qf = max(xMin,qf - correction);
        [~,cdfQ,pdfQ] = cf2DistGP2(cf,qf,[],options);
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

%% RESULT
result.x                  = x;
result.cdf                = cdf;
result.pdf                = pdf;
result.prob               = prob;
result.qf                 = qf;
result.SixSigmaRule       = options.SixSigmaRule;
result.N                  = N;
result.dt                 = dt;
result.T                  = t2(end);
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
result.details.count      = count;
result.details.correction = correction;
result.options            = options;
result.tictoc             = toc;

%% PLOT the PDF / CDF
if length(x)==1
    options.isPlot = false;
end
if options.isPlot
    figure
    plot(x,pdf,'.-')
    grid
    title('PDF Specified by the CF')
    xlabel('x')
    ylabel('pdf')
    %
    figure
    plot(x,cdf,'.-')
    grid
    title('CDF Specified by the CF')
    xlabel('x')
    ylabel('cdf')
end

end