function [testStat,result] = TestStat_MedianUniform(medians, options)
% TestStat_MedianUniform evaluates the GOF test statistic based on the
%  difference of the theoretical characteristic functions, cfM(t), and the
%  empirical characteristic functions, cfE(t), calculated from the
%  n observed medians from the N-dimesnsional random sample from the
%  uniform distribution U(a,b). Here, 
%
%    testStat = n*integral(abs(cfM(t)-cfE(t))^2*exp(-r*abs(t)^p),-Inf,Inf)
%
%  were r and p are the specified control parameters (as default values we
%  set r = 2, p = 2).
% 
% SYNTAX:
% [testStat,result] = TestStat_MedianUniform(medians, options)
%
% INPUT:
%  medians  - n-dimensional vector of observed iid medians
%  options  - structure with the following parameters:
%              options.a = -1;         % lower bound of the uniform
%                                      % distribution
%              options.b = 1;          % upper bound of the uniform
%                                      % distribution
%              options.r = 2;          % control parameter r 
%              options.p = 2;          % control parameter p
%              options.N = 10;         % size of the random sample from
%                                      % uniform distribution used to
%                                      % calculate the medians
%              options.Upp = 10;       % upper integration limit
%              options.cf_MedianUniform = [];  % function handle of the
%                                      % characteristic function of the
%                                      % sample median from the sample of
%                                      % size N from the uniform
%                                      % distribution U(a,b) 
%
% OUTPUT:
%  testStat - evaluated test statistic
%  result   - structure with further details.
%
% %EXAMPLE 1: Observed 1000 medians from sample of size N = 20 from U(-1,1)
%  n = 1000;
%  a = -1;
%  b = 1;
%  N = 20;
%  r = 2;
%  p = 2;
%  Upp = 5;
%  clear options
%  options.a = a;
%  options.b = b;
%  options.N = N;
%  options.r = r;
%  options.p = p;
%  options.Upp = Upp;
%  R = rand(n,N)*(b-a)+a; % Generated random values from U(a,b) 
%  M = median(R,2);
%  [testStat,result] = TestStat_MedianUniform(M,options)
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
% %EXAMPLE 2: Observed 50 medians from sample of size N = 20 from U(0,1)
%  n = 50;
%  a = 0;
%  b = 1;
%  N = 20;
%  r = 2;
%  p = 2;
%  Upp = 5;
%  clear options
%  options.a = a;
%  options.b = b;
%  options.N = N;
%  options.r = r;
%  options.p = p;
%  options.Upp = Upp;
%  R = rand(n,N)*(b-a)+a; % Generated random values from U(a,b) 
%  M = median(R,2);
%  [testStat,result] = TestStat_MedianUniform(M, options)
%  integrandFun = result.integrandFun;
%  cfM = result.cfM;
%  cfE = result.cfE;
%  t = linspace(-5,5);
%  figure
%  plot(t,integrandFun(t))
%  xlabel('t')
%  ylabel('function')
%  title('Integrand Function')
%  figure
%  t = linspace(-30,30);
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
validateattributes(medians,{'numeric'},{'vector'});
validateattributes(options,{'struct'},{});

if ~isfield(options, 'a')
    options.a = -1;
end

if ~isfield(options, 'b')
    options.b = 1;
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

if ~isfield(options, 'cf_MedianUniform')
    options.cf_MedianUniform = [];
end

%% Algorithm

medians = medians(:);
n = length(medians);

% Extract the input parameters    
a = options.a;
b = options.b;
N = options.N;
r = options.r;
p = options.p;
Upp = options.Upp;


% Calculate the theoretical characteristic function
if isempty(options.cf_MedianUniform)
    cfM = @(t) cf_MedianUniform(t, a, b, N);
    options.cf_MedianUniform = cfM;
else
    cfM = options.cf_MedianUniform;
end

% Calculate the theoretical characteristic function
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
    result.a = options.a;
    result.b = options.b;
    result.options = options;
end

end
%% Function integrand
function fun = integrand(t, cfM, cfE, n, r, p)

fun = n * abs(cfM(t)-cfE(t)).^2 .* exp(-r * abs(t).^p);

end