function [cdf,result] = cf2CDF_GPA(cf,x,options)
%cf2CDF_GPA  Evaluates the CDF from the characteristic function CF by the
%  Gil-Pelaez inversion formula. The algorithm uses the adaptive
%  Gauss-Kronrod quadrature rule for numerical integration of the
%  oscillatory integrand function divided into sub-intervals (found by a
%  fast root-finding algorithm) and subsequent application of the
%  convergence acceleration techniques for computing the limit of the
%  resulted alternating series, for more details see, e.g., Cohen et al.
%  (2000), Sidi (2011).
%
% NOTE
%  The algorithms cf2CDF_GPA/cf2PDF/cf2Qf are suggested as alternative to the
%  algorithm cf2DistGP based on using the simple trapezoidal for
%  integration, which can be inefficient for integration of highly
%  oscillatory integrand functions (which can appear if the
%  characteristic function decays very slowly to 0 for large t -> Infty,
%  e.g., for chi-square distribution with 1 degree of freedom). This is an
%  experimental (working) version of the algorithm cf2CDF_GPA.
%
% SYNTAX:
%  [cdf,result] = cf2CDF_GPA(cf,x,options)
%
% INPUT:
%  cf       - function handle of the characteristic function
%  x        - vector of x values where the CDF is computed
%  options  - structure with the following parameters:
%             options.isAccelerated    % indicator of the activated
%                                      % acceleration. 
%             options.nPeriods = 25    % the upper integration limit: UPPER
%                                      % = nPeriods * pi / x. The the basic
%                                      % integration interval [0,UPPER] is
%                                      % devided into two subintervals [0 A]
%                                      % and [A UPPER], where A is first
%                                      % zero of the integrand function.  
%             options.isPlot = false   % logical indicator for plotting the
%                                      % integrand function and calculation
%                                      % of their zeros.
%             options.tol = 1e-10      % absolute tolerance for numerical
%                                      % integration
%             options.tolFindRoots     % tolerance for the root-finding
%                                      % procedure. Default value is
%                                      % options.tolFindRoots = 1e-32.
%                                      % Alternatively, set lower levels up
%                                      % to options.tolFindRoots = 1e-321.
%             options.maxiterFindRoots % maximum number of iterations for
%                           = 100      % the root-finding procedure. 
%             options.nPoly = 2^5      % order of the chebyshev polynomials
%                                      % used for the root-finding
%                                      % procedure. 
% OUTPUT:
%  cdf      - vector of CDF values evaluated at x.
%  result   - structure with further details on computing the CDF.
%
% EXAMPLE:
% % CDF of the standard Normal distribution
%  cf  = @(t) exp(-t.^2/2);
%  x   = [1.281551565544601 1.644853626951472 2.326347874040841];
%  cdf = cf2CDF_GPA(cf,x)
%
% EXAMPLE:
% % CDF of the Chi-squared distribution with DF = 1
%  df = 1;
%  cf = @(t) (1 - 2i*t).^(-df/2);
%  x  = [2.705543454095416 3.841458820694126 6.634896601021214];
%  [cdf,result] = cf2CDF_GPA(cf,x)
%
% EXAMPLE:
% % CDF of the linear combination (convolution) of central Chi-squared RVs
%  df   = [1 2 3];
%  cf = @(t) cf_ChiSquare(t,df) ;
%  x = [10.644640675668420  12.591587243743977  16.811893829770927];
%  clear options
%  options.isPlot = true;
%  [cdf,result] = cf2CDF_GPA(cf,x,options)
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
% [6] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3),
%     32-44.
%
% SEE ALSO: cf2Dist, cf2DistGP, cf2DistGPT, cf2DistGPA, cf2DistFFT,
%           cf2DistBV, cf2CDF_GPA, cf2PDF, cf2QF

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 19-Dec-2018 21:59:07

%% CHECK/SET THE INPUT PARAMETERS
tic;
narginchk(2, 3);

if nargin < 3, options = []; end

if ~isfield(options, 'nPeriods')
    options.nPeriods = 25;
end

if ~isfield(options, 'isPlot')
    options.isPlot = false;
end

if ~isfield(options, 'division')
    options.division = [5 8 11];
end

if ~isfield(options, 'isAccelerated')
    options.isAccelerated = true;
end

if ~isfield(options, 'firstRootID')
    options.firstRootID = 1;
end

if ~isfield(options, 'nPoly')
    options.nPoly = 2^5;
end

if ~isfield(options, 'tol')
    options.tol = 1e-10;
end

if ~isfield(options, 'tolFindRoots')
    options.tolFindRoots = 1e-32;
end

if ~isfield(options, 'maxiterFindRoots')
    options.maxiterFindRoots = 100;
end

if ~isfield(options, 'verbose')
    options.verbose = false;
end

if ~isfield(options, 'xMean')
    options.xMean = [];
end

%% ALGORITHM
szx      = size(x);
x        = x(:);
nx       = length(x);
cdf      = zeros(nx,1);
error    = zeros(nx,1);
isPlot   = options.isPlot;
division = options.division;
xMean    = options.xMean;
for id = 1:nx
    xid  = x(id);
    fCDF = @(t) funCDF(t,cf,xid,xMean);
    div  = max(0.5,abs(xid));
    A    = (pi/2)/div;
    B    = pi*options.nPeriods/div;
    if options.isAccelerated
        try
            [roots,warnings] = FindRoots(fCDF,A,B,options.nPoly, ...
                isPlot,options.tolFindRoots,options.maxiterFindRoots);
        catch
            roots = [];
            warnings = 1;
        end
        if isempty(roots)
            R1 = A;
        elseif numel(roots) == 1
            R1 = roots(1);
        else
            if options.firstRootID < numel(roots)
                R1 = roots(options.firstRootID);
                Intervals = [roots(options.firstRootID:end-1)'; ...
                    roots(options.firstRootID+1:end)'];
            else
                R1 = roots(1);
                Intervals = [roots(1:end-1)';roots(2:end)'];
            end
        end
        nRoots = length(roots);
    else
        roots = [];
        nRoots = length(roots);
        warnings = 1;
    end
    if options.isAccelerated  && warnings == 0 && nRoots > 10
        isAcceleration = true;
        ERR1 = options.tol;
        INT1 = integral(fCDF,0,R1,'RelTol',0,'AbsTol',ERR1);
        [~,ERR2,~,Q] = IntegralGK(fCDF,Intervals,division,isPlot);
        INT  = INT1 + AcceleratedSum(sign(Q(1))*abs(Q));
        ERR = ERR1 + ERR2;
    else
        isAcceleration = false;
        ERR = options.tol;
        INT = integral(fCDF,0,Inf,'RelTol',0,'AbsTol',ERR);
    end
    cdf(id)   = 1/2 - INT/pi;
    error(id) = ERR;
end
cdf    = reshape(max(0,min(1,cdf)),szx);
x      = reshape(x,szx);
error  = reshape(error,szx);
tictoc = toc;

%% RESULTS
if nargout > 1
    result.Description         = 'CDF from the characteristic function CF';
    result.inversionMethod     = 'Gil-Pelaez';
    result.quadratureMethod    = 'adaptive Gauss-Kronrod with acceleration';
    result.CDF = cdf;
    result.x = x;
    if options.verbose
        result.ErrBound = error;
        result.ERR_xLast = ERR;
        result.warnings_xLast = warnings;
        result.isAcceleration_xLast = isAcceleration;
        result.fun_xLast = fCDF;
        result.roots_xLast = roots;
    end
    result.cf = cf;
    result.funCDF  = fCDF;
    result.options = options;
    result.tictoc  = tictoc;
end

end
%% Function funCDF
function fun = funCDF(t,cf,x,xMean)
%funCDF  Auxiliary function. Evaluates the integrand function for CDF

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 30-Nov-2018 11:14:06

fun = imag(cf(t).*exp(-1i*t*x))./t;

if ~isempty(xMean)
    fun(t==0) = xMean - x;
end

end