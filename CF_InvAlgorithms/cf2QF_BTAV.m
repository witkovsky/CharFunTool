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
% INPUT:
%  cf       - function handle of the characteristic function
%  prob     - vector of probabilities where the QF is computed
%  options  - structure with the following parameters:
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
%             options.fun_BTAV = 10    % parameter M for the deformed
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
%         disp(['df = ',num2str(df(i))])
%         nu = df(i)*ones(k,1);
%         cf = @(t) cfTest_Bartlett(t,nu);
%         qf = cf2QF_BTAV(cf,prob,options);
%         disp(qf(:)')
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

if ~isfield(options, 'fun_BTAV')
    options.fun_BTAV = 10;
end
%% ALGORITHM
szp        = size(prob);
prob       = prob(:)';
maxiter    = options.maxiter;
crit       = options.crit;
qf         = options.qf0;
n          = options.nTerms;
M          = options.fun_BTAV;
t          = linspace(0,pi,n+1)';
criterion  = true;
count      = 0;
while criterion
    count  = count + 1;
    [~,cdfFun,pdfFun] = fun_BTAV(t,qf,cf,[],M);
    CDF  = (2*pi/n)*real(cdfFun(1,:)/2 + sum(cdfFun(2:n,:)));
    PDF  = (2*pi/n)*real(pdfFun(1,:)/2 + sum(pdfFun(2:n,:)));
    correction  = (CDF - prob) ./ PDF;
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
result.CDF = CDF;
result.PDF = PDF;
result.nQuadratureTerms = n;
result.iterationCount = count;
result.lastCorrection = correction;
result.lastTermCDF = abs(cdfFun(n,:));
result.lastTermPDF = abs(pdfFun(n,:));
result.cdfFun = cdfFun;
result.pdfFun = pdfFun;
result.cf = cf;
result.options = options;
result.tictoc = tictoc;

end
%% function fun_BTAV
function [fun,cdfFun,pdfFun] = fun_BTAV(phi,x,cf,funtype,M)
%fun_BTAV 
%  Auxiliary function calculates the integrand function for computing the
%  CDF/PDF of the NON-NEGATIVE DISTRIBUTION specified by its characteristic
%  function, by using the Bromwich-Talbot-Abate-Valko (BTAV) inversion
%  method (originally suggested as numerical inversion for the Laplace
%  transform function). In particular, the  Bromwich contour is properly
%  deformed and the integrand functions are:
%    pdfFun = exp(s0.*x) .* cf(1i*s0) .* s1;
%    cdfFun = pdf ./ s0;
%  where cf is an anonymous characteristic function of nonnegative
%  distribution which is well defined for complex valued arguments, and
%    s0 = const * phi .* (cot(phi)+1i) ./ x
%    s1 = const * 1i * (1 + 1i*(phi + (phi.*cot(phi)-1).*cot(phi))) ./ x
%  for -pi <= phi <= pi.
%  Then the CDF/PDF can be numerically evaluated at specified values x by
%  integrating the integrand functions over the interval (-pi,pi) - by
%  using any quadrature rule (e.g. the simple trapezoidal or the more
%  advanced adaptive Gauss-Kronod quadrature rule):
%    CDF = real(integral(cdfFun,-pi,pi,'ArrayValued',true))
%    PDF = real(integral(pdfFun,-pi,pi,'ArrayValued',true))
%  For more details see Talbo (1979) and Abate & Valko (2004).
%
% SYNTAX:
%  [fun,cdfFun,pdfFun] = fun_BTAV(phi,x,cf,funtype,M)
%
% EXAMPLE 1
% % CDF/PDF of chi-square distribution with 1 degree of freedom, df = 1
% % By using integral - Matlab adaptive Gauss-Kronod quadrature
%  x      = [2.705543454095416 6.634896601021214 28.373987362798132];
%  cf     = @(t)(1-2i*t).^(-1/2);
%  cdfFun = @(phi) fun_BTAV(phi,x,cf,'cdf');
%  pdfFun = @(phi) fun_BTAV(phi,x,cf,'pdf');
%  CDF    = real(integral(cdfFun,-pi,pi,'ArrayValued',true));
%  PDF    = real(integral(pdfFun,-pi,pi,'ArrayValued',true));
%
% EXAMPLE 2 
% % CDF/PDF of chi-square distribution with 1 degree of freedom, df = 1
% % By using simple trapezoidal quadrature
%  x      = [2.705543454095416 6.634896601021214 28.373987362798132];
%  cf     = @(t)(1-2i*t).^(-1/2);
%  N      = 100;
%  phi    = linspace(-pi,pi,N+1)';
%  [~,cdfFun,pdfFun] = fun_BTAV(phi(2:end-1),x,cf);
%  CDF    = (2*pi/N)*real(sum(cdfFun));
%  PDF    = (2*pi/N)*real(sum(pdfFun));
%
% EXAMPLE 3 
% % CDF/PDF of chi-square distribution with 1 degree of freedom, df = 1
% % By using efficient trapezoidal quadrature on half interval
%  x      = [2.705543454095416 6.634896601021214 28.373987362798132];
%  cf     = @(t)(1-2i*t).^(-1/2);
%  N      = 50;
%  phi    = linspace(0,pi,N+1)';
%  [~,cdfFun,pdfFun] = fun_BTAV(phi,x,cf);
%  CDF    = (2*pi/N)*real(cdfFun(1,:)/2 + sum(cdfFun(2:end-1,:)));
%  PDF    = (2*pi/N)*real(pdfFun(1,:)/2 + sum(pdfFun(2:end-1,:)));
%
% REFERENCES:
% [1] Talbot, A., 1979. The accurate numerical inversion of Laplace
%     transforms. IMA Journal of Applied Mathematics, 23(1), pp.97-120.
% [2] Abate, J. and Valkó, P.P., 2004. Multi-precision Laplace transform
%     inversion. International Journal for Numerical Methods in
%     Engineering, 60(5), pp.979-993.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 22-Dec-2018 12:37:57

%% ALGORITHM
%[fun,cdfFun,pdfFun] = fun_BTAV(phi,x,cf,funtype,M)

%% CHECK THE INPUT PARAMETERS
narginchk(3, 5);

if nargin < 5, M = []; end
if nargin < 4, funtype = []; end

if isempty(M)
    M = 10;
end

if isempty(funtype)
    funtype = 'cdf';
end

Mp  = 2*M/5;
phi = phi(:);
x   = x(:)';
ct  = cot(phi);
s0  = phi .* (ct+1i);
s0(phi==0) = 1;
s0  = Mp * s0 ./ x;
s1  = 1 + 1i*(phi+(phi.*ct-1).*ct);
s1(phi==0) = 1;
s1  = s1 ./ x;

const  = Mp/(2*pi);
pdfFun = const * exp(s0.*x) .* cf(1i*s0) .* s1;
cdfFun = pdfFun ./ s0;

switch lower(funtype)
    case 'cdf'
        fun = cdfFun;
    case 'pdf'
        fun = pdfFun;
end
end