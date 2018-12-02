function [qf,result] = cf2QF(cf,prob,options)
%cf2QF Evaluates the quantiles (QF) from the distribution with 
%  given characteristic function CF by the Gil-Pelaez inversion formula
%  and application of the Newton-Raphson iterative algorithm. The algorithm
%  uses the adaptive Gauss-Kronod quadrature rule for numerical integration
%  of the oscillatory integrand function divided into sub-intervals (found
%  by a fast root-finding algorithm) and subsequent application of the
%  convergence acceleration techniques for computing the limit of the
%  resulted alternating series, for more details see e.g. Cohen et al.
%  (2000), Sidi (2011), Sidi (2017.)
%
% NOTE
%  The algorithms cf2CDF/cf2PDF/cf2Qf are suggested as alternative to the
%  algorithm cf2DistGP based on using the simple trapezoidal for
%  integration, which can be inefficient for integration of highly
%  oscillatory integrand functions (which can appear if the
%  characteristic function decays very slowly to 0 for large t -> Infty,
%  e.g., for chi-square distribution with 1 degree of freedom). This is an
%  experimental (working) version of the algorithm cf2QF.    
%
% SYNTAX:
% [qf,result] = cf2QF(cf,prob,options)
% INPUT:
%  cf       - function handle of the characteristic function
%  prob     - vector of probabilities where the QF is computed
%  options  - structure with the following parameters:
%             maxiter : indicator of maximum number of Newtob-Raphson
%             iterations. Default value is options.maxiter = 1000.
%             crit : value of the criterion limit for stopping rule.
%             Default value is options.crit = 1e-12.
%             verbose : indicator for producing result with a more detailed
%             output. Default value is options.verbose = false.
%             isPlot : logical indicator for plotting the integrand
%             function and calculation of their zeros. Default value is
%             options.isPlot = false.
%
% OUTPUT:
%  qf       - vector of the quantile values evaluated at prob.
%  result   - structure with QF further details:
%
% EXAMPLE:
% % QF of the standard Normal distribution
%  cf   = @(t) exp(-t.^2/2);
%  prob = [0.9 0.95 0.99];
%  qf   = cf2QF(cf,prob)
%
% EXAMPLE:
% % QF of the Chi-squared distribution with DF = 1
%  df = 1;
%  cf = @(t) (1 - 2i*t).^(-df/2);
%  prob = [0.9 0.95 0.99];
%  [qf,result] = cf2QF(cf,prob)
%
% EXAMPLE:
% % QF of the linear combination (convolution) of central Chi-squared RVs
%  df   = [1 2 3];
%  cf = @(t) cf_ChiSquare(t,df) ;
%  prob = [0.9 0.95 0.99];
%  [qf,result] = cf2QF(cf,prob)
%
% REFERENCES:
% [1] Gil-Pelaez, J., 1951. Note on the inversion theorem. Biometrika,
%     38(3-4), pp.481-482.
% [2] Imhof, J.: Computing the distribution of quadratic forms in normal
%     variables. Biometrika 48, 419–426 (1961).
% [3] Cohen, H., Villegas, F.R. and Zagier, D., 2000. Convergence
%     acceleration of alternating series. Experimental mathematics, 9(1),
%     pp.3-12.
% [4] Weniger, E.J., 1989. Nonlinear sequence transformations for the
%     acceleration of convergence and the summation of divergent series.
%     Computer Physics Reports, 10(5-6), pp.189-371. 
% [5] Sidi, A., 2011. A user-friendly extrapolation method for computing
%     infinite range integrals of products of oscillatory functions. IMA
%     Journal of Numerical Analysis, 32(2), pp.602-631.
% [6] Sidi, A., 2017. Acceleration of Convergence of Some Infinite
%     Sequences $\boldsymbol {\{A_n\}} $ Whose Asymptotic Expansions
%     Involve Fractional Powers of $\boldsymbol {n} $. arXiv preprint
%     arXiv:1703.06495.
% [7] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3),
%     32-44.
%
% SEE ALSO: cf2Dist, cf2DistGP, cf2DistGPT, cf2DistGPA, cf2DistFFT,
%           cf2DistBV, cf2CDF, cf2PDF, cf2QF

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Dec-2018 13:36:33

%% CHECK/SET THE INPUT PARAMETERS
StartTime = cputime;
narginchk(2, 3);
if nargin < 3, options = []; end

if ~isfield(options, 'verbose')
    options.verbose = false;
end

if ~isfield(options, 'qf0')
    options.qf0 = (cf(1e-4)-cf(-1e-4))/(2e-4*1i);
end

if ~isfield(options, 'crit')
    options.crit = 1e-12;
end

if ~isfield(options, 'maxiter')
    options.maxiter = 1000;
end

if ~isfield(options, 'isPlot')
    options.isPlot = false;
end

options.isPlot = false;

%% ALGORITHM
szp       = size(prob);
prob      = prob(:);
maxiter   = options.maxiter;
crit      = options.crit;
qf        = options.qf0;
criterion = true;
count     = 0;
while criterion
    count = count + 1;
    CDFqf = cf2CDF(cf,qf,options);
    PDFqf = cf2PDF(cf,qf,options);
    correction  = (CDFqf - prob) ./ PDFqf;
    qf = qf - correction;
    criterion = any(abs(correction) > crit * abs(qf)) ...
        && max(abs(correction)) > crit && count < maxiter;
end
qf       = reshape(qf,szp);
prob     = reshape(prob,szp); 
tictoc   = cputime - StartTime;

%% RESULTS
result.quantile = qf;
result.prob = prob;
result.CDF = CDFqf;
result.PDF = PDFqf;
result.correction = correction;
result.count = count;
result.cf = cf;
result.options = options;
result.tictoc = tictoc;

end