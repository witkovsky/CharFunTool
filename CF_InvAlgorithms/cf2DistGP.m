function [result,cdf,pdf,qf] = cf2DistGP(cf,x,prob,options)
%cf2DistGP Calculates the CDF/PDF/QF from the characteristic function CF
% by using the Gil-Pelaez inversion formulae:
%   cdf(x) = 1/2 + (1/pi) * Integral_0^inf imag(exp(-1i*t*x)*cf(t)/t)*dt.
%   pdf(x) = (1/pi) * Integral_0^inf real(exp(-1i*t*x)*cf(t))*dt.
%
%  The FOURIER INTEGRALs are calculated by using the simple TRAPEZOIDAL
%  QUADRATURE method, see below.
% 
% ALIAS NAME:
%  cf2DistGP is an alias name for the algorithm cf2DistGPT.
%
% SYNTAX:
%  result = cf2DistGP(cf,x)
%  [result,cdf,pdf,qf] = cf2DistGP(cf,x,prob,options)
%
% INPUT:
%  cf      - function handle of the characteristic function (CF), 
%  x       - vector of x values where the CDF/PDF is computed, 
%  prob    - vector of values from [0,1] for which the quantiles
%            function is evaluated,
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
%  options.isCompound = true;
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
%  options.isCompound = true;
%  options.isInterp = true;
%  result = cf2DistGP(cf,x,prob,options)
%  PDF = result.PDF
%  CDF = result.CDF
%  QF  = result.QF
%  RND = result.RND
%
% REFERENCES:
% [1] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference, 2001, 94,
%     1–13. 
% [2] WITKOVSKY, V.: Matlab algorithm TDIST: The distribution of a
%     linear combination of Student’s t random variables. In COMPSTAT
%     2004 Symposium (2004), J. Antoch, Ed., Physica-Verlag/Springer
%     2004, Heidelberg, Germany, pp. 1995–2002.
% [3] WITKOVSKY, V., WIMMER, G., DUBY, T. Logarithmic Lambert W x F
%     random variables for the family of chi-squared distributions
%     and their applications. Statistics & Probability Letters, 2015, 96,
%     223–231. 
% [4] WITKOVSKY, V.: Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 2016, 5(3),
%     32-44. 
% [5] WITKOVSKY, V., WIMMER, G., DUBY, T. Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. ArXiv preprint,
%     2017, arXiv:1701.08299.
%
% SEE ALSO: cf2Dist, cf2DistGP, cf2DistGPT, cf2DistGPA, cf2DistFFT,
%           cf2DistBV, cf2CDF, cf2PDF, cf2QF

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Dec-2018 20:23:52

%% ALGORITHM

narginchk(1, 4);
if nargin < 4, options = []; end
if nargin < 3, prob = []; end
if nargin < 2, x = []; end

[result,cdf,pdf,qf] = cf2DistGPT(cf,x,prob,options);

end