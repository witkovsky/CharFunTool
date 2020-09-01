function [result,cdf,pdf,qf] = cf2Dist(cf,x,prob,options,algorithm)
%cf2Dist Calculates the CDF/PDF/QF from the characteristic function CF
% by using the inversion algorithm specified by algorithm (string).
% 
% ALIAS NAME:
%  cf2Dist is an alias name for the algorithm specified by ALGORITHM name. 
%
% SYNTAX:
%  result = cf2Dist(cf,x,prob,options)
%  [result,cdf,pdf,qf] = cf2Dist(cf,x,prob,options,algorithm)
%
% INPUT:
%  cf        - function handle of the characteristic function (CF), 
%  x         - vector of x values where the CDF/PDF is computed, 
%  prob      - vector of values from [0,1] for which the quantiles
%              function is evaluated,
%  options   - options structure of the selected inversion algorithm
%  algorithm - selected algorithm for numerical inversion of the
%              characteristic function. Currently available inversion
%              algorithms include: 
%              'cf2DistBTAV'  (Bromwich-Talbot-Abate-Valko method), 
%              'cf2DistBV',   (Bakhvalov-Vasilieva method)
%              'cf2DistFFT',  (Fast Fourier Transform algorithm method)
%              'cf2DistGP',   (Gil-Pelaez with Trapezoidal quadrature)
%              'cf2DistGPT',  (Gil-Pelaez with Trapezoidal quadrature)
%              'cf2DistGPA'   (Gil-Pelaez with adaptive Gauss-Kronrod
%                              quadrature and convergence acceleration 
%                              techniques).  
%              If empty, default value is algorithm = 'cf2DistGP'.
% OUTPUT:
%  result   - structure with CDF/PDF/QF and further details,
%  cdf      - vector of CDF values evaluated at x,
%  pdf      - vector of PDF values evaluated at x,
%  qf       - vector of QF values evaluated at prob.
%
% EXAMPLE1 (Calculate CDF/PDF of N(0,1) by inverting its CF)
%  cf = @(t) exp(-t.^2/2);
%  result = cf2Dist(cf)
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
%  algorithm = 'cf2DistGPT';
%  result = cf2Dist(cf,x,prob,options,algorithm)
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
%  algorithm = 'cf2DistGPA';
%  result = cf2Dist(cf,x,prob,options,algorithm)
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
% SEE ALSO: cf2Dist, cf2DistGP, cf2DistGPT, cf2DistGPA, cf2CDF, cf2PDF,
%           cf2QF 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-Sep-2020 13:25:21
%
% Revision history:
% Ver.: 02-Dec-2018 20:23:52

%% ALGORITHM

narginchk(1, 5);
if nargin < 5, algorithm = []; end
if nargin < 4, options = []; end
if nargin < 3, prob = []; end
if nargin < 2, x = []; end

if isempty(algorithm)
    [result,cdf,pdf,qf] = cf2DistGP(cf,x,prob,options);
else
    [result,cdf,pdf,qf] = eval([algorithm,'(cf,x,prob,options)']);
end

end