%function [result,cdf,pdf,qf] = cf2DistFFT(cf,x,prob,options)
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
%             options.isCompound = false % treat the compond distributions
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
%[result,cdf,pdf,qf] = cf2DistFFT(cf,x,prob,options);