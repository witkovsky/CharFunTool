function result = PolyCal(x,y,xFit,options)
% PolylCal - Comparative polynomial calibration. PolyCal estimates the
% patrameters of the polynomial calibration function together with 
% their uncertainties. 
%
% In metrology, the straight-line calibration and the polynomial
% calibration problem is addressed in the ISO Technical Specifications
% ISO/TS 28037:2010 [1] and ISO/TS 28038:2018 [2] which are based on the
% GUM uncertainty framework [3] using the law of propagation of
% uncertainty. 
%
% PolyCal is an alternative approach for estimating the calibration
% parameters together with their associated uncertainties and, moreover,
% with their derived state-of-knowledge distributions, based on the EIV
% (Errors-In-Variables) modeling approach [4] and the characteristic
% function approach [5-6]. For detailed description of the suggested
% concept for polynomial comparative calibration see [7-8].  
%
% Here, in the considered calibration experiment, the direct measurements
% taken by two independent measurement devices are represented by
% m-dimensional vectors X = mu + (XA + XB + XB0) and  Y = nu + (YA + YB +
% YB0), where mu = (mu_1,...,mu_m)' and nu = (nu_1,...,nu_m)' represent the
% true unknown values of the measurands expressed in units of the measuring
% devices. These parameters are functionally related by the polynomial
% calibration function, nu = beta_0 + beta_1 mu^1 + ... + beta_p mu^p. The
% components of the measurement errors are assumed to be mutually
% independent zero-mean random variables, fully specified by their
% probability distributions (here by their characteristic functions)derived
% by Type A and Type B methods of evaluation. 
% 
% PolyCal estimates the patrameters of the polynomial calibration function
% beta = (beta_0,...,beta_p) together with their uncertainties (covariance
% matrix) and their marginal state-of-knowledge distributions. Moreover,
% the algorithm evaluates the fitted values yFit of the calibration
% function for given xFit, yFit =  beta_0 + beta_1 xFit^1 + ... + beta_p
% xFit^p, together with their uncertainties and the state-of-knowledge
% distributions. The output from the PolyCal algorithm can be further used
% to evaluate the results of the future measurements by the calibrated
% instrument, together with their uncertainties and the state-of-knowledge
% distributions.
% 
% SYNTAX:
% result = PolyCal(x,y,xFit,options)
% 
% INPUTS:
%  x      - m-dimensional vector of measurements taken by the instrument X.
%  y      - m-dimensional vector of measurements taken by the instrument Y.
%  xFit   - n-dimensional vector where the polynomial calibration function
%           will be evaluated, yFit =  beta_0 + beta_1 xFit^1 + ... + beta_p
%           xFit^p. If options.xFitWeights is specified as non-empty
%           n-dimensional vector of weights, the resulted yFit is calculated
%           as a weighted mean of the particular yFit(j) for j = 1,...,n.
% options - structure with the optional parameters:
%           options.uXA - m-dimensional vector with the uncertainties
%                   of the independent components of the vector XA =
%                   (XA_1,...,XA_m). If empty, uXA is calculated from cfXA.
%           options.uXB - m-dimensional vector with the uncertainties
%                   of the independent components of the vector XB =
%                   (XB_1,...,XB_m). If empty, uXB is calculated from cfXB.
%           options.uXB0 - uncertainty of the independent random variable
%                   XB0, common for all components of X = (X_1,...,X_m). If
%                   empty, uXB0 is calculated from cfXB0. 
%           options.uYA - m-dimensional vector with the uncertainties
%                   of the independent components of the vector YA =
%                   (YA_1,...,YA_m). If empty, uYA is calculated from cfYA.
%           options.uYB - m-dimensional vector with the uncertainties
%                   of the independent components of the vector YB =
%                   (YB_1,...,YB_m). If empty, uYB is calculated from cfYB.
%           options.uYB0 - uncertainty of the independent random variable
%                   YB0, common for all components of Y = (Y_1,...,Y_m). If
%                   empty, uYB0 is calculated from cfYB0. 
%           options.cfXA - m-dimensional cell with function handles of the
%                   characteristic functions of independent components of
%                   XA = (XA_1,...,XA_m). 
%           options.cfXB - m-dimensional cell with function handles of the
%                   characteristic functions of independent components of
%                   XB = (XB_1,...,XB_m). 
%           options.cfXB0 - function handle of the characteristic function
%                   of the random variable XB0.
%           options.cfYA - m-dimensional cell with function handles of the
%                   characteristic functions of independent components of
%                   YA = (YA_1,...,YA_m). 
%           options.cfYB - m-dimensional cell with function handles of the
%                   characteristic functions of independent components of
%                   YB = (YB_1,...,YB_m). 
%           options.cfYB0 - function handle of the characteristic function
%                   of the random variable YB0.
%           options.deltaXB - m-dimensional vector with the known Type B
%                   nonstochastic corrections for the components of X =
%                   (X_1,...,X_m). 
%           options.deltaXB0  - the known Type B nonstochastic correction
%                   common for all components of X = (X_1,...,X_m). 
%           options.deltaYB- m-dimensional vector with the known Type B
%                   nonstochastic corrections for the components of Y =
%                   (Y_1,...,Y_m). 
%           options.deltaYB0 - the known Type B nonstochastic correction
%                   common for all components of Y = (Y_1,...,Y_m).
%           options.order - order of the polynomial calibration function.
%                   Hence, the dimension of the parameter vector beta0 is
%                   (order+1). 
%           options.Aijk - the (m,order+1,q)-dimensional array of
%                   coefficients Aijk, which specifies the generalized
%                   polynomial calibration model. 
%           options.q - third dimension of the (m,order+1,q)-dimensional
%                   array of coefficients Aijk, which specifies the
%                   generalized polynomial calibration model.
%           options.alpha - nominal significance level. Default value is
%                   options.alpha = 0.05.
%           options.probLow - lower probability used to evaluate confidence
%                   intervals for the fitted values. Im empty,
%                   options.probLow = alpha/2.
%           options.probUpp - upper probability used to evaluate confidence
%                   intervals for the fitted values. Im empty,
%                   options.probUpp = 1 - alpha/2.
%           options.yLow - n-dimensional vector of lower bounds of the
%                   confidence intevals used to be inverted by the
%                   Bonferroni method.   
%           options.yUpp - n-dimensional vector of upper bounds of the
%                   confidence intevals used to be inverted by the
%                   Bonferroni method. 
%           options.isplot - flag indicator for ploting the default plots.
%           options.mu0 - m-dimensional vector of initial values of the
%                   vector parameter mu (mean value of the measurement
%                   vector X) used for linrarization of the calibration
%                   model. 
%           options.nu0 - m-dimensional vector of initial values of the
%                   vector parameter nu (mean of the measurement vector Y)
%                   used for linrarization of the calibration model.
%           options.beta0 - (order+1)-dimensional vector of initial values
%                   of the vector parameter beta (parameters of the
%                   polynomial calibration function of given order) used
%                   for linrarization of the calibration model.
%           options.isIntFit - flag indicator for computing the coverage
%                   intervals of the fitted values yFit calculated from
%                   their state-of-knowledge distributions.  
%           options.xFitWeights - n-dimensional vector of weights (the size
%                   should be equal to xFit). If non-empty, the resulted
%                   yFit is calculated as a weighted mean of the particular
%                   yFit(j) for j = 1,...,n.
%           options.tolDiff - selected small difference used for numerical
%                   differenciation. If empty we set here options.tolDiff =
%                   1e-2.  
%
% OUTPUTS:
%  result - structure with the following parameters:
%  % Estimated parameters
%          result.Description - Description of the method.
%          result.muEstimate - m-dimensional vector of the estimates of the
%                 parameter mu. 
%          result.nuEstimate - m-dimensional vector of the estimates of the
%                 parameter nu. 
%          result.betaEstimate - (order+1)-dimensional vector of the
%                 estimates of the calibration function parameter beta. 
%          result.uBeta - (order+1)-dimensional vector of the
%                 uncertainties of the calibration function parameter beta. 
%          result.SigmaBeta - (order+1) x (order+1)-dimensional covariance
%                 matrix of the calibration function parameter beta. 
%          result.cfBeta - (order+1)-dimensional cell with function handles
%                 of the characteristic functions of the state-of-knowledge
%                 distrfibutions of the components of the estimator
%                 betaEstimate.
%  % Other useful inputs and outputs
%          result.order - order of the polynomial calibration function.
%          result.m - dimension of the measurement vectors x and y.
%          result.LX - [(order+1) x m]-dimensional coefficient matrix used
%                 for estimation of betaEstimate.
%          result.KY - [(order+1) x m]-dimensional coefficient matrix used
%                 for estimation of betaEstimate.
%          result.uXA   - see options.uXA.
%          result.uXB   - see options.uXB.
%          result.uXB0  - see options.uXB0.
%          result.uYA   - see options.uYA.
%          result.uYB   - see options.uYB.
%          result.uYB0  - see options.uYB0.
%          result.cfXA  - see options.cfXA.
%          result.cfXB  - see options.cfXB.
%          result.cfXB0 - see options.cfXB0.
%          result.cfYA  - see options.cfYA.
%          result.cfYB  - see options.cfYB.
%          result.cfYB0 - see options.cfYB0.
%  % Fitted values
%          result.xFit - n-dimensional vector of values where the
%                 calibration function is evaluated. If empty, the default
%                 value is xFit = x. 
%          result.yFit - n-dimensional vector of fitted values. 
%          result.uFit - n-dimensional vector of the uncertainties of the
%                 fitted values. 
%          result.cfFit - n-dimensional cell with function handles
%                 of the characteristic functions of the state-of-knowledge
%                 distrfibutions of the components of the fitted yFit.
%          result.IntFit - (n x 2)-dimensional matrix with coverage
%                 intervals of the components of the fitted yFit.
%  % Bonferroni confidence intervals
%          result.yLow - see options.yLow
%          result.uUpp - see options.yUPP
%          result.xLow -  n-dimensional vector of lower bounds of the
%                   inverted confidence intevals by the Bonferroni method.
%          result.xUpp -  n-dimensional vector of upper bounds of the
%                   inverted confidence intevals by the Bonferroni method. 
%  % OTHER inputs and outputs
%          result.x - m-dimensional vector of measurements by the
%                 instrument X.
%          result.y - m-dimensional vector of measurements by the
%                 instrument Y.
%          result.SigmaX - (m x m)-dimensional covariance matrix of the
%                 random vector of the measurements X. 
%          result.SigmaY- (m x m)-dimensional covariance matrix of the
%                 random vector of the measurements Y. 
%          result.ucX - m-dimensional vector of uncertainties of the
%                 measurements X. 
%          result.ucY - m-dimensional vector of uncertainties of the
%                 measurements Y. 
%          result.q    - see options.q.
%          result.Aijk - see optionsAijk.
%          result.A - matrix A used for final linearization of polynomial
%                 comparative calibration model of given order, locally at  
%                 the specified mu0, nu0, and beta0.
%          result.B - matrix B used for final linearization of polynomial
%                 comparative calibration model of given order, locally at  
%                 the specified mu0, nu0, and beta0.
%          result.c - vector c used for final linearization of polynomial
%                 comparative calibration model of given order, locally at  
%                 the specified mu0, nu0, and beta0.
%          result.D - matrix D used for final linearization of polynomial
%                 comparative calibration model of given order, locally at  
%                 the specified mu0, nu0, and beta0.
%          result.W - evaluated matrix W = D*SigmaX*D' + SigmaY.
%          result.Q - (m +(order+1)) x (m +(order+1)) -dimensional matrix
%                 calculated in the final iteration. The blocks of the
%                 matrix Q are usefull, e.g. to express covariance matrices
%                 of the estimators.
%          result.crit - final value of the criterion crit.
%          result.tolerance - used tolerance. 
%          result.loops - number of iteration loops.
%          result.options - copy of the option structure.
%
% EXAMPLE 1 (simple)
%  % Comparative Polynomial Calibration - Artificial Data
%  x = [4.7, 5.5, 6.5, 7.5, 8.4]';
%  y = [5.8, 7.9, 10.4, 12.7, 15.2]';
%  xFit = linspace(5,8,11)';
%  clear options
%  options.order = 1;
%  options.cfXA = {@(t)cf_Normal(0.01*t), ...
%          @(t)cf_Normal(0.1*t), ...
%          @(t)cf_Normal(0.1*t), ...
%          @(t)cf_Student(0.1*t,5), ...
%          @(t)cf_RectangularSymmetric(0.2*t)};
%  options.cfXB = {@(t)cf_RectangularSymmetric(0.1*t), ...
%          @(t)cf_RectangularSymmetric(0.1*t), ...
%          @(t)cf_TriangularSymmetric(0.1*t), ...
%          @(t)cf_TriangularSymmetric(0.1*t), ...
%          @(t)cf_RectangularSymmetric(0.2*t)};
%  options.cfXB0 = [];
%  options.cfYA = {@(t)cf_Normal(0.15*t), ...
%          @(t)cf_Normal(0.2*t), ...
%          @(t)cf_Normal(0.2*t), ...
%          @(t)cf_Normal(0.3*t), ...
%          @(t)cf_Normal(0.3*t)};
%  options.cfYB = {@(t)cf_RectangularSymmetric(0.2*t), ...
%          @(t)cf_RectangularSymmetric(0.2*t), ...
%          @(t)cf_TriangularSymmetric(0.3*t), ...
%          @(t)cf_TriangularSymmetric(0.3*t), ...
%          @(t)cf_RectangularSymmetric(0.3*t)};
%  options.cfYB0 = @(t)cf_ArcsineSymmetric(0.1*t);
%  options.tolDiff = 1e-2;
%  result = PolyCal(x,y,xFit,options)
%
% EXAMPLE 2 (continued EXAMPLE 1)
%  % Plot the PDF/CDF of the parameters beta0 and beta1
%  cfBeta = result.cfBeta;
%  resultBeta0 = cf2DistGP(cfBeta{1},[],[],options);
%  figure;
%  resultBeta1 = cf2DistGP(cfBeta{2},[],[],options);
%  
% EXAMPLE 3 (advanced)
%  % CF of a weighted linear combination (a weighted mean)
%  % of the yFit(j) fitted from xFit(j) for j = 1,...,n, where the weights
%  % are given  by options.xFitWeights.
%  x = [4.7, 5.5, 6.5, 7.5, 8.4]';
%  y = [5.8, 7.9, 10.4, 12.7, 15.2]';
%  xNewBestEstimate = 7;
%  cfxFit = @(t) exp(1i*t*xNewBestEstimate) .* cf_Student(0.1*t,10);
%  [prob,weights] = LegendrePoints(21,0,1);
%  clear options
%  options.isPlot = false;
%  options.SixSigmaRule = 10;
%  [~,~,~,xFit] = cf2DistGP(cfxFit,[],prob,options);
%  options.xFitWeights = weights;
%  options.order = 1;
%  options.cfXA = {@(t)cf_Normal(0.01*t), ...
%          @(t)cf_Normal(0.1*t), ...
%          @(t)cf_Normal(0.1*t), ...
%          @(t)cf_Student(0.1*t,5), ...
%          @(t)cf_RectangularSymmetric(0.2*t)};
%  options.cfXB = {@(t)cf_RectangularSymmetric(0.1*t), ...
%          @(t)cf_RectangularSymmetric(0.1*t), ...
%          @(t)cf_TriangularSymmetric(0.1*t), ...
%          @(t)cf_TriangularSymmetric(0.1*t), ...
%          @(t)cf_RectangularSymmetric(0.2*t)};
%  options.cfXB0 = [];
%  options.cfYA = {@(t)cf_Normal(0.15*t), ...
%          @(t)cf_Normal(0.2*t), ...
%          @(t)cf_Normal(0.2*t), ...
%          @(t)cf_Normal(0.3*t), ...
%          @(t)cf_Normal(0.3*t)};
%  options.cfYB = {@(t)cf_RectangularSymmetric(0.2*t), ...
%          @(t)cf_RectangularSymmetric(0.2*t), ...
%          @(t)cf_TriangularSymmetric(0.3*t), ...
%          @(t)cf_TriangularSymmetric(0.3*t), ...
%          @(t)cf_RectangularSymmetric(0.3*t)};
%  options.cfYB0 = @(t)cf_ArcsineSymmetric(0.1*t);
%  options.tolDiff = 1e-2;
%  result = PolyCal(x,y,xFit,options);
%  cfyFit = result.cfyFit{1};
%  % Plot the Characteristic function of yFit  
%  figure
%  t = linspace(-10,10,2^9);
%  plot(t,real(cfyFit(t)),t,imag(cfyFit(t)))
%  title('Characteristic Function of the Combined yFit')
%  % Plot the PDF/CDF of yFit  
%  figure
%  options.isPlot = true;
%  resultyFit = cf2DistGP(cfyFit,[],[],options);
%
% REFERENCES
%
% [1] ISO/TS 28037:2010. Determination and Use of Straight-Line Calibration
%     Functions. International StandardsOrganization, Geneva, September
%     (2010).   
% [2] ISO/TS 28038:2018. Determination and Use of Polynomial Calibration
%     Functions. International Standards Organization, Geneva, December
%     (2018). 
% [3] JCGM 100:2008. Evaluation of measurement data – Guide to the
%     expression of uncertainty in measurement (GUM 1995 with minor
%     corrections), ISO, BIPM, IEC, IFCC, ILAC, IUPAC, IUPAP and OIML,
%     (2008).  
% [4] KUBÁČEK L. Foundations of Estimation Theory, Elsevier, Amsterdam
%     (1988). 
% [5] WITKOVSKÝ V. Numerical inversion of a characteristic function: An
%     alternative tool to form the probability distribution of output
%     quantity in linear measurement models. ACTA IMEKO 5 (3), (2016)
%     32–44.   
% [6] WITKOVSKÝ V. CharFunTool: The Characteristic Functions Toolbox
%     (MATLAB). https://github.com/witkovsky/CharFunTool, (2020).
% [7] WITKOVSKÝ V. and WIMMER G. Generalized polynomial comparative
%     calibration: Parameter estimation and applications. In: Advances in
%     Measurements and Instrumentation: Reviews. S.Y. Yurish (Ed.).
%     Barcelona, Spain, IFSA (International Frequency Sensor Association
%     Publishing) Publishing S.L., (2018) 15-52. 
% [8] WITKOVSKÝ V. and WIMMER G. PolyCal - MATLAB algorithm for comparative
%     polynomial calibration and its applications. AMCTM 2020.  

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Sep-2020 16:46:23

%% Check the Inputs
narginchk(2, 4);

if nargin < 4, options = []; end
if nargin < 3, xFit = []; end

if ~isfield(options, 'uXA')
    options.uXA = [];
end

if ~isfield(options, 'uXB')
    options.uXB = [];
end

if ~isfield(options, 'uXB0')
    options.uXB0 = [];
end

if ~isfield(options, 'uYA')
    options.uYA = [];
end

if ~isfield(options, 'uYB')
    options.uYB = [];
end

if ~isfield(options, 'uYB0')
    options.uYB0 = [];
end

if ~isfield(options, 'cfXA')
    options.cfXA = [];
end

if ~isfield(options, 'cfXB')
    options.cfXB = [];
end

if ~isfield(options, 'cfXB0')
    options.cfXB0 = [];
end

if ~isfield(options, 'cfYA')
    options.cfYA = [];
end

if ~isfield(options, 'cfYB')
    options.cfYB = [];
end

if ~isfield(options, 'cfYB0')
    options.cfYB0 = [];
end

if ~isfield(options, 'deltaXB')
    options.deltaXB = [];
end

if ~isfield(options, 'deltaXB0')
    options.deltaXB0 = [];
end

if ~isfield(options, 'deltaYB')
    options.deltaYB = [];
end

if ~isfield(options, 'deltaYB0')
    options.deltaYB0 = [];
end

if ~isfield(options, 'alpha')
    options.alpha = 0.05;
end

if ~isfield(options, 'isplot')
    options.isplot = true;
end

if ~isfield(options, 'maxiter')
    options.maxiter = 100;
end

if ~isfield(options, 'tolerance')
    options.tolerance = 1e-12;
end

if ~isfield(options, 'order')
    options.order = [];
end

if ~isfield(options, 'Aijk')
    options.Aijk = [];
end

if ~isfield(options, 'q')
    options.q = [];
end

if ~isfield(options, 'mu0')
    options.mu0 = [];
end

if ~isfield(options, 'nu0')
    options.nu0 = [];
end

if ~isfield(options, 'beta0')
    options.beta0 = [];
end

if ~isfield(options, 'isIntFit')
    options.isIntFit = true;
end

if ~isfield(options, 'probLow')
    options.probLow = [];
end

if ~isfield(options, 'probUpp')
    options.probUpp = [];
end

if ~isfield(options, 'yLow')
    options.yLow = [];
end

if ~isfield(options, 'yUpp')
    options.yUpp = [];
end

if ~isfield(options, 'isXY')
    options.isXY = true; 
end

if ~isfield(options, 'xFitWeights')
    options.xFitWeights = [];
end

if ~isfield(options, 'tolDiff')
    options.tolDiff = [];
end

isXY = options.isXY;
if isXY
    uXA      = options.uXA;
    uXB      = options.uXB;
    uXB0     = options.uXB0;
    uYA      = options.uYA;
    uYB      = options.uYB;
    uYB0     = options.uYB0;
    cfXA     = options.cfXA;
    cfXB     = options.cfXB;
    cfXB0    = options.cfXB0;
    cfYA     = options.cfYA;
    cfYB     = options.cfYB;
    cfYB0    = options.cfYB0;
    deltaXB  = options.deltaXB;
    deltaXB0 = options.deltaXB0; 
    deltaYB  = options.deltaYB;
    deltaYB0 = options.deltaYB0;  
    x        = x(:);
    y        = y(:);
else
    uXA      = options.uYA;
    uXB      = options.uYB;
    uXB0     = options.uY0;
    uYA      = options.uXA;
    uYB      = options.uXB;
    uYB0     = options.uXB0;
    cfXA     = options.cfYA;
    cfXB     = options.cfYB;
    cfXB0    = options.cfYB0;
    cfYA     = options.cfXA;
    cfYB     = options.cfXB;
    cfYB0    = options.cfXB0;
    deltaXB  = options.deltaYB;
    deltaXB0 = options.deltaYB0; 
    deltaYB  = options.deltaXB;
    deltaYB0 = options.deltaXB0; 
    x        = y(:);
    y        = x(:);
end

%% % Set the Inputs
order       = options.order;
Aijk        = options.Aijk;
q           = options.q;
alpha       = options.alpha;
probLow     = options.probLow;
probUpp     = options.probUpp;
yLow        = options.yLow;
yUpp        = options.yUpp;
isplot      = options.isplot;
mu0         = options.mu0;
nu0         = options.nu0;
beta0       = options.beta0;
isIntFit    = options.isIntFit;
xFitWeights = options.xFitWeights;
tolDiff     = options.tolDiff;

m       = length(x);
one1    = ones(m,1);
zero0   = zeros(m,1);

if isempty(tolDiff)
    tolDiff = 1e-2;
end

% Uncertainties of the input variables
if isempty(uXA) && ~isempty(cfXA)
    uXA = PolyCalCF2Std(cfXA,tolDiff);
elseif isempty(uXA)
    uXA = zero0; 
end

if isempty(uXB) && ~isempty(cfXB)
    uXB = PolyCalCF2Std(cfXB,tolDiff);
elseif isempty(uXB)
    uXB = zero0; 
end

if isempty(uXB0) && ~isempty(cfXB0)
    uXB0 = PolyCalCF2Std(cfXB0,tolDiff);
elseif isempty(uXB0)
    uXB0 = 0; 
end

if isempty(uYA) && ~isempty(cfYA)
    uYA = PolyCalCF2Std(cfYA,tolDiff);
elseif isempty(uYA)
    uYA = one1; 
end

if isempty(uYB) && ~isempty(cfYB)
    uYB = PolyCalCF2Std(cfYB,tolDiff);
elseif isempty(uYB)
    uYB = zero0; 
end

if isempty(uYB0) && ~isempty(cfYB0)
    uYB0 = PolyCalCF2Std(cfYB0,tolDiff);
elseif isempty(uYB0)
    uYB0 = 0; 
end

% Systematic shifts
if isempty(deltaXB), deltaXB = zero0; end
if isempty(deltaXB0), deltaXB0 = 0; end
if isempty(deltaYB), deltaYB = zero0; end
if isempty(deltaYB0), deltaYB0 = 0; end
x    = (x(:) - deltaXB(:)) - deltaXB0;
y    = (y(:) - deltaYB(:)) - deltaYB0;
xMin = min(x);
xMax = max(x);

% Characteristic functions of the input variables
if isempty(cfXA)
    cfXA  = cell(m,1); 
    for i = 1:m
        cfXA{i} = @(t) exp(-(uXA(i)*t).^2/2); % ~ N(0,uXA^2)
    end
end

if isempty(cfXB)
    cfXB  = cell(m,1); 
    for i = 1:m
        cfXB{i} = @(t) exp(-(uXB(i)*t).^2/2); % ~ N(0,uXB^2)
    end
end

if isempty(cfXB0)
    cfXB0 = @(t) exp(-(uXB0*t).^2/2); % ~ N(0,uXB0^2)
end

if isempty(cfYA)
    cfYA  = cell(m,1); 
    for i = 1:m
        cfXA{i} = @(t) exp(-(uYA(i)*t).^2/2); % ~ N(0,uYA^2)
    end
end

if isempty(cfYB)
    cfYB  = cell(m,1); 
    for i = 1:m
        cfYB{i} = @(t) exp(-(uYB(i)*t).^2/2); % ~ N(0,uYB^2)
    end
end

if isempty(cfYB0)
    cfYB0 = @(t) exp(-(uYB0*t).^2/2); % ~ N(0,uYB0^2)
end

% Covariance matrices of the measurement vectors X and Y
SigmaX    = diag(uXA.^2) + diag(uXB.^2) + uXB0^2 * ones(m);
SigmaY    = diag(uYA.^2) + diag(uYB.^2) + uYB0^2 * ones(m);

% Uncertainties of the the measurements X and Y
ucX     = sqrt(diag(SigmaX));
ucY     = sqrt(diag(SigmaY));

if isempty(order)
    order = 1;
end

if isempty(xFit)
    xFit = linspace(xMin,xMax)';
end

if isempty(mu0)
    mu0 = x;
end

if isempty(nu0)
    nu0 = y;
end

if isempty(q)
    if ~isempty(Aijk)
        [~,~,q] = size(Aijk);
    end
end

if isempty(beta0)
    if isempty(Aijk)
        beta0 = flip(polyfit(mu0,nu0,order))';
    else
        B0 = update(mu0,nu0,beta0,m,order,q,Aijk);
        beta0 = B0\nu0;
    end
end

if isempty(probUpp)
    probUpp = 1-alpha/2;
end

if isempty(probLow)
    probLow = alpha/2;
end

%% Algorithm

% Initialization
maxiter      = options.maxiter;
tolerance    = options.tolerance;
crit         = 1;
loops        = 1;
muEstimate   = mu0;
nuEstimate   = nu0;
betaEstimate = beta0;
id1          = 1:m;
id2          = (m+1):(m+order+1);

% Iteration
while (crit > tolerance) && (loops <= maxiter)
    mu0          = muEstimate;
    nu0          = nuEstimate;
    beta0        = betaEstimate;
    muOld        = muEstimate;
    nuOld        = nuEstimate;
    betaOld      = betaEstimate;
    [B,D]        = PolyCalUpdate(mu0,nu0,beta0,m,order,q,Aijk);
    W            = D*SigmaX*D + SigmaY;
    Q            = [W B; B' zeros(order+1)]\eye(m+order+1);
    Q11          = Q(id1,id1);
    Q21          = Q(id2,id1);
    yDx          = y - D*(x-mu0);
    muEstimate   = x + SigmaX*D*Q11*yDx;
    nuEstimate   = y - SigmaY*Q11*yDx;
    betaEstimate = Q21*yDx; % = (B'*Winv*B)\(B'*Winv*yDx)
    crit = norm([muEstimate;nuEstimate;betaEstimate] - ... 
        [muOld;nuOld;betaOld])^2/(2*m+order+1); 
    loops = loops+1;
end

% Final result
mu0          = muEstimate;
nu0          = nuEstimate;
beta0        = betaEstimate;
[B,D,A,c]    = PolyCalUpdate(mu0,nu0,beta0,m,order,q,Aijk);
W            = D*SigmaX*D + SigmaY;
Q            = [W B; B' zeros(order+1)]\eye(m+order+1);
Q11          = Q(id1,id1);
Q21          = Q(id2,id1);
Q22          = Q(id2,id2);
KY           = Q21;
LX           = -KY*D;
SigmaBeta    = -Q22;
uBeta        = sqrt(diag(SigmaBeta));
yDx          = y - D*(x-mu0);
muEstimate   = x + SigmaX*D*Q11*yDx;
nuEstimate   = y - SigmaY*Q11*yDx;
betaEstimate = KY*yDx; 

% Calculate the CFs of the polynomial parameters
cfBeta = cell(order+1,1);
for i = 1:(order+1)
    wX = LX(i,:);
    wY = KY(i,:);
    cfBeta{i} = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,...
        betaEstimate(i));
end

% Calculate the fitted values
xFit  = xFit(:);
nxFit = length(xFit);
w     = ones(nxFit,order+1);

for i = 2:(order+1)
    w(:,i) = w(:,i-1).*xFit;
end

if isempty(xFitWeights) 
    % If xFitWeights is empty calculate cfFit{j} for each xFit(j) 
    yFit    = zeros(nxFit,1);
    uyFit   = zeros(nxFit,1);
    cfyFit  = cell(nxFit,1);
    cfyFit0 = cell(nxFit,1);
    for j = 1:nxFit
        yFit(j)    = w(j,:)*betaEstimate;
        uyFit(j)   = sqrt(w(j,:)*SigmaBeta*w(j,:)');
        wX         = w(j,:)*LX;
        wY         = w(j,:)*KY;
        cfyFit{j}  = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,yFit(j));
        cfyFit0{j} = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0);
    end
else
    % If xFitWeights is specified, resulted CF is weighted mean from cfFit{j}
    yFit   = 0;
    uyFit  = 0;
    cfyFit{1} = @(t) 0;
    cfyFit0{1} = @(t) 0;
    for j = 1:nxFit
        yFit  = w(j,:)*betaEstimate;
        uyFit = uyFit + (w(j,:)*SigmaBeta*w(j,:)')*xFitWeights(j)^2;
        wX    = w(j,:)*LX;
        wY    = w(j,:)*KY;
        cfPolyCal = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,yFit);
        cfyFit{1} = @(t) cfyFit{1}(t) + xFitWeights(j) * cfPolyCal(t);
        cfPolyCal = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0);
        cfyFit0{1} = @(t) cfyFit0{1}(t) + xFitWeights(j) * cfPolyCal(t);
    end
    uyFit = sqrt(uyFit);
end

if isIntFit
    nx = length(cfyFit);
    IntFit = zeros(nx,2);
    prob  = [alpha/2,1-alpha/2];
    options.isPlot = false;
    for j = 1:nx
        res = cf2DistGP(cfyFit{j},[],prob,options);
        IntFit(j,:) = res.qf;
    end
else
    IntFit = [];
end

% Bonferroni-type confidence intervals for x given intervals [yLow,yUpp]
xLow = [];
xUpp = [];
if ~isempty(yLow) && ~isempty(yUpp)
    nyLow = length(yLow);
    nyUpp = length(yUpp);
    xLow  = zeros(nyLow,1);
    xUpp  = zeros(nyUpp,1);
    if nyLow == nyUpp
        for i = 1:nyLow
            [xL,xU] = PolyCalInterval(yLow,yUpp,betaEstimate,LX,KY,...
                xMin,xMax,probLow,probUpp,order,cfXA,cfXB,cfXB0,...
                cfYA,cfYB,cfYB0,options);
            xLow(i) = xL;
            xUpp(i) = xU;
        end
    else
        warning('The size of yLow does not fit with the size of yUpp')
    end
end

%% Results
% Estimated parameters
result.Description = 'PolyFit - Polynomial Comparative Calibration';
result.muEstimate   = muEstimate;
result.nuEstimate   = nuEstimate;
result.betaEstimate = betaEstimate;
result.uBeta        = uBeta;
result.SigmaBeta    = SigmaBeta;
result.cfBeta       = cfBeta;
% Other useful inputs and outputs
result.order        = order;
result.m            = m;
result.LX           = LX;
result.KY           = KY;
result.uXA          = uXA;
result.uXB          = uXB;
result.uXB0         = uXB0;
result.uYA          = uYA;
result.uYB          = uYB;
result.uYB0         = uYB0;
result.cfXA         = cfXA;
result.cfXB         = cfXB;
result.cfXB0        = cfXB0;
result.cfYA         = cfYA;
result.cfYB         = cfYB;
result.cfYB0        = cfYB0;
% Fitted values
result.xFit         = xFit;
result.yFit         = yFit;
result.uyFit        = uyFit;
result.cfyFit       = cfyFit;
result.cfyFit0      = cfyFit0;
result.IntFit       = IntFit;
% Bonferroni confidence intervals
result.yLow         = yLow;
result.uUpp         = yUpp;
result.xLow         = xLow;
result.xUpp         = xUpp;
% OTHER inputs and outputs
result.x            = x;
result.y            = y;
result.SigmaX       = SigmaX;
result.SigmaY       = SigmaY;
result.Q            = Q;
result.ucX          = ucX;
result.ucY          = ucY;
result.q            = q;
result.Aijk         = Aijk;
result.A            = A;
result.B            = B;
result.c            = c;
result.D            = D;
result.W            = W;
result.crit         = crit;
result.tolerance    = tolerance;
result.loops        = loops;
result.options      = options;

%% Plot
if isplot
    figure
    plot(x,y,'o')
    hold on
    for i = 1:m
        plot([x(i)-2*ucX(i),x(i)+2*ucX(i)],[y(i),y(i)],'r-')
        plot([x(i),x(i)],[y(i)-2*ucY(i),y(i)+2*ucY(i)],'r-')
    end
    plot(muEstimate,nuEstimate,'x-')
    plot(xFit,yFit,'b+-')
    if isempty(IntFit)
        plot(xFit,yFit+2*uyFit,'r-.')
        plot(xFit,yFit-2*uyFit,'r-.')
    else
        plot(xFit,IntFit(:,1),'b-.')
        plot(xFit,IntFit(:,2),'b-.')
        %
        plot(xFit,yFit+1.96*uyFit,'r-.')
        plot(xFit,yFit-1.96*uyFit,'r-.')
    end
    hold off
    grid on
end

end