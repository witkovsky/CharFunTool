function [qf,result] = cf2QF_BTAV(cf,prob,options)
%cf2QF_BTAV
%  Evaluates the quantiles (QF) from the NON-NEGATIVE DISTRIBUTION
%  specified by the given characteristic function CF by using the
%  Bromwich-Talbot-Abate-Valko (BTAV) inversion method (originally
%  suggested as numerical inversion for the Laplace transform function).
%  Here we assume that the specified CF is a characteristic function of
%  nonnegative distribution which is well defined for complex valued
%  arguments.
% 
% SYNTAX:
% [qf,result] = cf2QF_BTAV(cf,prob,options)
%
% INPUT:
%  cf       - function handle of the characteristic function
%  prob     - vector of probabilities where the QF is computed
%  options  - structure with the following parameters:
%             options.quadrature       % quadrature method. Default method
%                  = 'trapezoidal'     % is quadrature = 'trapezoidal'.
%                                      % Alternatively use quadrature =
%                                      % 'matlab' (for the MATLAB built in
%                                      % adaptive Gauss-Kronod quadrature) 
%                                      % Gauss-Kronod integral. 
%             options.tol = 1e-12      % absolute tolerance for the MATLAB
%                                      % Gauss-Kronod integral. 
%             options.crit = 1e-12     % value of the criterion limit for
%                                      % stopping rule.
%             options.nTerms = 50      % number of terms used in the
%                                      % trapezoidal quadrature
%             options.maxiter = 100    % indicator of the maximum number of
%                                      % Newton-Raphson iterations.
%             options.qf0              % starting values of the quantiles.
%                                      % By default, the algorithm starts
%                                      % from the mean of the distribution
%                                      % estimated from the specified CF
%             options.Mpar_BTAV = 10   % parameter M for the deformed
%                                      % Bromwich curve.
%
% OUTPUT:
%  qf       - vector of the quantile values evaluated at prob.
%  result   - structure with further details.
%
% EXAMPLE 1
% % QF of the Chi-squared distribution with DF = 1
%  df = 1;
%  cf = @(t) (1 - 2i*t).^(-df/2);
%  prob = [0.9 0.95 0.99];
%  [qf,result] = cf2QF_BTAV(cf,prob)
%
% EXAMPLE 2 
% % QF of the linear combination (convolution) of central Chi-squared RVs
%  df   = [1 2 3];
%  cf = @(t) cf_ChiSquare(t,df) ;
%  prob = [0.9 0.95 0.99];
%  [qf,result] = cf2QF_BTAV(cf,prob)
%
% EXAMPLE 3:
% % Quantiles of Bartlett distribution computed by the algorithm cf2QF_BTAV
%     k     = 5;
%     df    = [1 2 3 4 5 6 7 8 9 10];
%     prob  = [0.9 0.95 0.99];
%     Table = zeros(10,3);
%     for i = 1:10
%         nu = df(i)*ones(k,1);
%         cf = @(t) cfTest_Bartlett(t,nu);
%         qf = cf2QF_BTAV(cf,prob);
%         Table(i,:) = qf;
%     end
%     disp(Table)
%
% NOTE OF CAUTION
%  The method was suggested for inverting proper Laplace tranform
%  functions. Here we use available characteristic functions for
%  creatig Laplace transform function, by using M(s) = cf(1i*s). In
%  general, the implemented algorithms for computing CFs assume that the
%  argument t is real. In specific situations CF is well defined also for
%  complex arguments. However, numerical issues could appear in any step
%  during the calculations. The result and the inner calculations should be
%  chcecked and properly controlled. 
%
% REFERENCES:
% [1] Talbot, A., 1979. The accurate numerical inversion of Laplace
%     transforms. IMA Journal of Applied Mathematics, 23(1), pp.97-120.
% [2] Abate, J. and Valkó, P.P., 2004. Multi-precision Laplace transform
%     inversion. International Journal for Numerical Methods in
%     Engineering, 60(5), pp.979-993.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 22-Dec-2018 23:56:13

%% CHECK/SET THE INPUT PARAMETERS
StartTime = cputime;
narginchk(2, 3);
if nargin < 3, options = []; end

if ~isfield(options, 'quadrature')
    options.quadrature = 'trapezoidal';
end

if ~isfield(options, 'tol')
    options.tol =1e-12;
end

if ~isfield(options, 'crit')
    options.crit = 1e-12;
end

if ~isfield(options, 'nTerms')
    options.nTerms = 50;
end

if ~isfield(options, 'maxiter')
    options.maxiter = 100;
end

if ~isfield(options, 'qf0')
    options.qf0 = (cf(1e-4)-cf(-1e-4))/(2e-4*1i);
end

if ~isfield(options, 'Mpar_BTAV')
    options.Mpar_BTAV = 10;
end

%% ALGORITHM
szp        = size(prob);
prob       = prob(:)';
maxiter    = options.maxiter;
crit       = options.crit;
qf         = options.qf0;
criterion  = true;
count      = 0;
while criterion
    count      = count + 1;
    [cdf,pdf]  = cf2CDFPDF_BTAV(cf,qf,options);
    correction = (cdf - prob) ./ pdf;
    qf = qf - correction;
    criterion = max(abs(cdf - prob)) > 30*eps && ...
        max(abs(correction)) > crit && count < maxiter;
end
qf       = reshape(qf,szp);
prob     = reshape(prob,szp);
tictoc   = cputime - StartTime;

%% RESULTS
result.inversionMethod = 'Bromwich-Talbot-Abate-Valkó';
result.quadratureMethod = options.quadrature;
result.quantile = qf;
result.prob = prob;
result.cdf  = cdf;
result.pdf  = pdf;
result.nIterations = count;
result.lastCorrection = correction;
result.cf = cf;
result.options = options;
result.tictoc = tictoc;

end