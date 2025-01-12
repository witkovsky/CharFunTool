function [testStat,result] = TestStat_MedianExponential(medians, options)
% TestStat_MedianExponential evaluates the GOF test statistic based on the
%  difference of the theoretical characteristic functions, cfM(t), and the
%  empirical characteristic functions, cfE(t), calculated from the
%  n observed medians from the N-dimensional random sample from the
%  exponential distribution EXP(lambda). Here, 
%
%    testStat = n*integral(abs(cfM(t)-cfE(t))^2*exp(-r*abs(t)^p),-Inf,Inf)
%
%  were r and p are the specified control parameters (as default values we
%  set r = 2, p = 2).
% 
% SYNTAX:
% [testStat,result] = TestStat_MedianExponential(medians, options)
%
% INPUT:
%  medians  - n-dimensional vector of observed iid medians
%  options  - structure with the following parameters:
%              options.lambda = 1;     % parameter of the exponential
%                                      % distribution
%              options.r = 2;          % control parameter r 
%              options.p = 2;          % control parameter p
%              options.N = 10;         % size of the random sample from
%                                      % exponential distribution used to
%                                      % calculate the medians
%              options.Upp = 10;       % upper integration limit
%              options.cf_MedianExponential = [];  % function handle of the
%                                      % characteristic function of the
%                                      % sample median from the sample of
%                                      % size N from the exponential
%                                      % distribution EXP(lambda)
%
% OUTPUT:
%  testStat - evaluated test statistic
%  result   - structure with further details.
%
% %EXAMPLE 1: Observed 1000 medians from sample of size N = 20 from EXP(1)
%  n = 1000;
%  N = 20;
%  lambda = 1;
%  R = exprnd(1/lambda,n,N); % Generated random values from Exponential 
%  M = median(R,2);
%  clear options
%  options.lambda = lambda; 
%  options.N = N;
%  options.r = 2;
%  options.p = 4;
%  options.Upp = 5;
%  [testStat,result] = TestStat_MedianExponential(M,options)
%  integrandFun = result.integrandFun;
%  cfM = result.cfM;
%  cfE = result.cfE;
%  t = linspace(-3,3,101);
%  figure
%  plot(t,integrandFun(t))
%  xlabel('t')
%  ylabel('function')
%  title('Integrand Function')
%  figure
%  t = linspace(-20,20,201);
%  plot(t, real(cfM(t)), t, imag(cfM(t)))
%  hold on
%  plot(t, real(cfE(t)), t, imag(cfE(t)))
%  grid
%  xlabel('t')
%  ylabel('CF')
%  title('Theoretical and Empirical Characteristic Function')
%  legend({'cfM real' 'cfM imag' 'cfE real' 'cfE imag'},'Location','ne')
%
% %EXAMPLE 2: Observed 15 medians from sample of size N = 7 from EXP(5)
%  n = 30;
%  N = 20;
%  lambda = 0.25;
%  R = exprnd(1/lambda,n,N); % Generated random values from Exponential 
%  M = median(R,2);
%  clear options
%  options.lambda = lambda; 
%  options.N = N;
%  options.r = 1;
%  options.p = 2;
%  options.Upp = 5;
%  [testStat,result] = TestStat_MedianExponential(M, options)
%  integrandFun = result.integrandFun;
%  cfM = result.cfM;
%  cfE = result.cfE;
%  t = linspace(-5,5,101);
%  figure
%  plot(t,integrandFun(t))
%  xlabel('t')
%  ylabel('function')
%  title('Integrand Function')
%  figure
%  t = linspace(-30,30,201);
%  plot(t, real(cfM(t)), t, imag(cfM(t)))
%  hold on
%  plot(t, real(cfE(t)), t, imag(cfE(t)))
%  grid
%  xlabel('t')
%  ylabel('CF')
%  title('Theoretical and Empirical Characteristic Function')
%  legend({'cfM real' 'cfM imag' 'cfE real' 'cfE imag'},'Location','ne')

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Jan-2024 13:28:02

%% CHECK/SET THE INPUT PARAMETERS

narginchk(1, 2);
if nargin < 2, options = []; end

% Validate the input arguments
validateattributes(medians,{'numeric'},{'vector','nonnegative'});
validateattributes(options,{'struct'},{});

if ~isfield(options, 'lambda')
    options.lambda = 1;
end

if ~isfield(options, 'r')
    options.r = 2;
end

if ~isfield(options, 'p')
    options.p = 2;
end

if ~isfield(options, 'N')
    options.N = 10;
end

if ~isfield(options, 'Upp')
    options.Upp = 10;
end

if ~isfield(options, 'cf_MedianExponential')
    options.cf_MedianExponential = [];
end

%% Algorithm

% Extract the input parameters
medians = medians(:);
n = length(medians);
N = options.N;
r = options.r;
p = options.p;
Upp = options.Upp;

% Calculate the theoretical characteristic function
if isempty(options.cf_MedianExponential)
    lambda = options.lambda;
    cfM = @(t) cf_MedianExponential(t, lambda, N);
    options.cf_MedianExponential = cfM;
else
    cfM = options.cf_MedianExponential;
end

% Calculate the empirical characteristic function
cfE = @(t) cfE_Empirical(t,medians);

% Set the integrand: fun = n * abs(cfM(t)-cfE(t)).^2 .* exp(-r * abs(t).^p)
integrandFun = @(t) integrand(t, cfM, cfE, n, r, p);

% Evaluate the test statistic
testStat = 2 * integral(integrandFun, 0, Upp);

% Store the results in a structure
if nargout > 1
    result.testStat = testStat;
    result.integrandFun = integrandFun;
    result.cfE = cfE;
    result.cfM = cfM;
    result.ObservedMedians = medians;
    result.N = options.N;
    result.n = n;
    result.lambda = options.lambda;
    result.options = options;
end

end
%% Function integrand
function fun = integrand(t, cfM, cfE, n, r, p)

fun = n * abs(cfM(t)-cfE(t)).^2 .* exp(-r * abs(t).^p);

end