function result = cf2DIST(cf,options)
%CF2PDF  Evaluates the PDF/CDF from the characteristic function CF
%        by the Fast Fourier Transform (FFT).
%
% SYNTAX:
% result = cf2DIST(cf,options)
%
% INPUT:
% cf       - function handle for the characteristic
%            function.
% options  - structure with the following parameters:
%   xmin   - minimum (support value) x for the evaluated PDF,
%            in order to get good result xmin should be such that
%            pdf(xmin) = 0, or pdf(xmin) ~ 0.
%   xmax   - maximum (support value) x for the evaluated PDF,
%            in order to get good result xmax should be such that
%            PDF(xmax) = 0, or PDF(xmax) ~ 0. PDF is evaluated
%            for n equidistant points in [xmin,xmax].
%   n      - number of points, a power of 2, default
%            value is options.n = 2^10.
%   isPlot - logical indicator for plotting the PDF, default
%            value is options.isPlot = true.
%
% OUTPUT:
% result   - structure with the following parameters:
%   pdf    - n dimensional vector of PDF values
%            evaluated at x.
%   x      - n dimensional vector of points where the
%            PDF was evaluated.
%   cf     - n dimensional vector of CHF values evaluated at t.
%   t      - n dimensional vector of points where the CHF was
%            evaluated.
%
% EXAMPLE 1:
% % PDF OF A LINEAR COMBINATION OF THE INDEPENDENT RVs
% % (Normal, Uniform, Triangular, and Student's t distribution)
% % Y = 1*X_N + 10*X_U + 2*X_T + 0.5*X_t
% % Characteristic functions of X_N, X_U, X_T, X_t
% chf_N  = @(t) exp(-t.^2/2);              % Normal(0,1)
% chf_U  = @(t) min(1,sin(t)./t);          % Uniform(-1,1)
% chf_T  = @(t) min(1,(2-2*cos(t))./t.^2); % Triangular(-1,1)
% nu     = 1;                              % degrees of freedom
% chf_t  = @(t) min(1,besselk(nu/2, abs(t).*sqrt(nu),1) ...
%          .* exp(-abs(t).*sqrt(nu)) ...
%          .* (sqrt(nu).*abs(t)).^(nu/2) ...
%          / 2^(nu/2-1)/gamma(nu/2));      % Student's t
% % Characteristic function of the linear combination Y
% chf_Conv = @(t) chf_N(1*t) .* chf_U(10*t) ...
%            .* chf_T(2*t) .* chf_t(0.5*t);
% options.xmin = -50;
% options.xmax = 50;
% options.n    = 2^10;
% options.isZeroSymmetric = true;
% result = cf2DIST(chf_Conv,options);
% title('PDF of Y = 1*X_{N} + 10*X_{U} + 2*X_{T} + 0.5*X_{t}')
%
% EXAMPLE 2:
% % PDF OF A LINEAR COMBINATION OF THE INDEPENDENT RVs
% % (Chi-squared with 1 and 10 degrees of freedom)
% % Y = 10*X_Chi2_1 + 1*X_Chi2_10
% % Characteristic functions of X_Chi2_1, X_Chi2_10
% nu1 = 1;
% chf_Chi2_1  = @(t) (1-2i*t).^(-nu1/2);
% nu2 = 10;
% chf_Chi2_10 = @(t) (1-2i*t).^(-nu2/2);
% % Characteristic function of the linear combination Y
% chf_Conv = @(t) chf_Chi2_1(10*t) .* chf_Chi2_10(t);
% options.isForcedSymmetric = true;
% result = cf2DIST(chf_Conv,options);
% title('PDF of Y = 10*X_{\chi^2_{1}} + 1*X_{\chi^2_{10}}')
%
% REFERENCES:
%
% WITKOVSKÝ , V. On the exact computation of the density and of
% the quantiles of linear combinations of t and F random
% variables. Journal of Statistical Planning and Inference 94
% (2001), 1–13.
%
% WITKOVSKÝ , V. Matlab algorithm TDIST: The distribution of a
% linear combination of Student’s t random variables. In COMPSTAT
% 2004 Symposium (2004), J. Antoch, Ed., Physica-Verlag/Springer
% 2004, Heidelberg, Germany, pp. 1995–2002.
%
% WITKOVSKÝ , V., WIMMER,G., DUBY, T. Logarithmic Lambert W x F
% random variables for the family of chi-squared distributions
% and their applications. Statistics & Probability Letters 96
% (2015), 223–231.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 12-Mar-2015 10:27:22

%% CHECK THE INPUT PARAMETERS
tic;
narginchk(1, 2);

if nargin < 2, options = []; end

if ~isfield(options, 'n')
    options.n = 2^10;
end

if ~isfield(options, 'xmin')
    options.xmin = [];
end

if ~isfield(options, 'xmax')
    options.xmax = [];
end

if ~isfield(options, 'QuadratureRule')
    options.QuadratureRule = 'MPR';
end

if ~isfield(options, 'rho')
    options.rho = 0.15;
end

if ~isfield(options, 'SixSigmaRule')
    options.SixSigmaRule = 6;
end

if ~isfield(options, 'pdf0')
    options.pdf0 = NaN;
end

if ~isfield(options, 'isForcedSymmetric')
    options.isForcedSymmetric = false;
end

if ~isfield(options, 'isZeroSymmetric')    
    if options.isForcedSymmetric
        options.isZeroSymmetric = true;
    else
        options.isZeroSymmetric = false;
    end
end

if ~isfield(options, 'isChebfun')
    options.isChebfun = false;
end

if ~isfield(options, 'nChebPts')
    options.nChebPts = 100;
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

%% SET THE INTEGRATION RULE
rule = options.QuadratureRule;
switch lower(rule)
    case 'lpr'
        Qrule = [0,   1];
    case 'mpr'
        Qrule = [1/2, 1];
    case 'rpr'
        Qrule = [1,   1];
    case 'simpson'
        Qrule = [0,   1/6;
                 1/2, 4/6;
                 1,   1/6];
    case 'o3'
        Qrule = [0,   1/4;
                 1/2, 2/4;
                 1,   1/4];
    case 'trapezoid'
        Qrule = [0,   1/2;
                 1,   1/2];
    case 'gauss3'
        Qrule = [0.112701665379258,   0.277777777777778;
                 0.500000000000000,   0.444444444444444;
                 0.887298334620742,   0.277777777777778];
    case 'gauss5'
        Qrule = [0.046910077030668,   0.118463442528095;
                 0.230765344947158,   0.239314335249683;
                 0.500000000000000,   0.284444444444444;
                 0.769234655052842,   0.239314335249683;
                 0.953089922969332,   0.118463442528095];
    case 'gauss7'
        Qrule = [0.025446043828621,   0.064742483084435;
                 0.129234407200303,   0.139852695744638;
                 0.297077424311301,   0.190915025252560;
                 0.500000000000000,   0.208979591836735;
                 0.702922575688699,   0.190915025252560;
                 0.870765592799697,   0.139852695744638;
                 0.974553956171379,   0.064742483084435];
    otherwise
        Qrule = [1/2, 1];
        warning('Unknown Quadrature Rule - USED rule: Middle Point Rule')
end

%% MOMENTS AND SUPPORT (the approximate distribution range)
cfd = @(t) numDiff(cf,t,1e-4);
mu1 = real(numDiff(cf)/1i);
mu2 = real(numDiff(cfd)/1i^2);

if options.isForcedSymmetric
    Mean  = 0;
    Var   = mu2;
    Std   = sqrt(Var);
    Upper = Mean + options.SixSigmaRule*Std;
    Lower = max(0,Mean - options.SixSigmaRule*Std);
else
    Mean  = mu1;
    Var   = mu2 - mu1^2;
    Std   = sqrt(Var);
    Upper = Mean + options.SixSigmaRule*Std;
    Lower = Mean - options.SixSigmaRule*Std;
end

if isempty(options.xmax) && isempty(options.xmax)
    options.xmax = Upper;
    options.xmin = Lower;
end

%% ALGORITHM
cf_org = cf;
if options.isForcedSymmetric
    cf  = @(t) real(cf(t));
    N   = 2*options.n;
    B   = options.xmax;     % Approximate upper limit of sym(X)
    A   = -B;               % Approximate lower limit of sym(X)
else
    N   = options.n;        % Number of sub-intervals
    A   = options.xmin;     % Approximate lower limit of X
    B   = options.xmax;     % Approximate upper limit of X
end

k    = 0:(N-1);              % k indices 0 : N-1
j    = 0:(N-1);              % j indices 0 : N-1
bp   = Qrule(:,1);           % Quadrature base points
w    = Qrule(:,2);           % Quadrature weights
dx   = (B-A)/N;              % dx is width of the subintervals [x_i,x_{i+1}]
x    = A + k * dx;           % left points of the subintervals [x_i,x_{i+1}]
sft  = dx*(w'*bp);
xCDF = x + sft;           % shift the left points x to the "middle"

nBP = numel(bp);
pdf = 0;
% PDF from CF by FFT (see e.g. Hürlimann 2013, Witkovsky 2015)
for id  = 1:nBP
    C   = ((-1).^((A/(B-A) + k/N)*(N-2*bp(id))))/(B-A);
    u   = ((j + bp(id)) - N/2)/(B-A);
    chf = cf(2*pi*u);
    phi = (-1).^(-(2*A/(B-A))*j) .* chf;
    pdf = pdf + w(id) * real(C .* fft(phi));
end
pdf(isnan(pdf)) = 0;

% PDF / CDF
pdf = max(0,pdf);
cdf = cumsum(pdf*dx);
cdf = max(0,cdf);
cdf = min(1,cdf);

if options.isZeroSymmetric
    cdf = cdf + 1/2-(cdf(N/2+1)+cdf(N/2))/2;
end

if options.isForcedSymmetric
    x2PDF   = x(N/2+1:end);
    %x2PDF(1)= 0;
    x2CDF   = xCDF(N/2:end);
    x2CDF(1)= 0;
    pdf2    = 2*pdf(N/2+1:end);
    pdf2(1) = options.pdf0;
    % pdf2(2) = pdf2(2)/2;
    cdf2    = 2*cdf(N/2:end)-1;
    cdf2(1) = 0;
    A       = 0;
else
    x2PDF   = x;
    x2CDF   = xCDF;
    pdf2    = pdf;
    cdf2    = cdf;
end

%% CHEBFUN PDF/CDF (create chebfun at n chebyshev points over [A,B])
PDF = [];
CDF = [];
QF  = [];
if options.isChebfun
    pts = chebpts(options.nChebPts,[A,B]);
    f = interp1(x2PDF,pdf2,pts(2:end-1));
    pdfVal = reshape([options.pdf0;f(:);0],size(pts));
    pdfVal(isnan(pdfVal)) = 0;
    pdfVal = max(0,pdfVal);
    PDF  = chebfun(pdfVal,[A,B]);
    CDF  = cumsum(PDF);
end

%% QUANTILES
% if options.isChebfun
%     pts       = chebpts(options.nChebPts,[1e-4,1-1e-4]);
% else
%     pts = linspace(1e-4,1-1e-4,101);
% end
% 
% prob      = pts;
% maxiter   = 100;
% count     = 0;
% crit      = 1e-12;
% criterion = true;
% qf  = mu1 * ones(size(prob));
% while criterion
%     count  = count + 1;
%     correction  = (interp1(x2CDF,cdf2,qf) - prob) ./ interp1(x2PDF,pdf2,qf);
%     qf = qf - correction;
%     criterion = any(abs(correction) > crit * abs(qf)) ...
%         && max(abs(correction)) > crit && count < maxiter; 
% end
% qf(end) = Inf;
% 
% if options.isChebfun
%     QF  = chebfun(qf,[1e-4,1-1e-4]);
% end

%% RESULT
mytime         = toc;
result.pdf     = pdf;
result.cdf     = cdf;
result.xPDF       = x;
result.pdf2    = pdf2;
result.cdf2    = cdf2;
result.x2PDF   = x2PDF;
result.x2CDF   = x2CDF;
%result.prob    = prob;
%result.qf      = qf;
result.dx      = dx;
result.N       = N;
result.chf     = chf;
result.t       = 2*pi*u;
result.A       = A;
result.B       = B;
result.tmin    = 2*pi*u(1);
result.tmax    = 2*pi*u(end);
result.cf      = cf_org;
result.cf_used = cf;
result.PDF     = PDF;
result.CDF     = CDF;
result.QF      = QF;
result.Mean    = mu1;
result.Var     = mu2-mu1^2;
result.Std     = sqrt(result.Var);

result.QRule   = options.QuadratureRule;
result.options = options;
result.tictoc  = mytime;

%% PLOT
if options.isPlot
    figure
    plot(x2PDF,pdf2,'.-')
    grid
    title('PDF Specified by the Characteristic Function CHF')
    xlabel('x')
    ylabel('pdf')
    
    figure
    plot(x2CDF,cdf2,'.-')
    grid
    title('CDF Specified by the Characteristic Function CHF')
    xlabel('x')
    ylabel('cdf')
end

end
function d1 = numDiff(fun,x,h)
% NUMDIFF  Auxiliary function to calculate numerical (first) derivative.
%          Given below is the five point method for the first derivative
%          (five-point stencil in one dimension). See  Abramowitz & Stegun,
%          Table 25.2.
%          See also http://en.wikipedia.org/wiki/Numerical_differentiation.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 12-Mar-2015 08:36:01

if nargin < 3
    h = 1e-3;
end

if nargin < 3
    x = 0;
end

d1 = (-fun(x+2*h)+8*fun(x+h)-8*fun(x-h)+fun(x-2*h))/(12*h);
end
%% END

