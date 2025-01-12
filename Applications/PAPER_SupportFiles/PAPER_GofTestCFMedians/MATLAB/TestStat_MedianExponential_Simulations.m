function testStat = TestStat_MedianExponential_Simulations(M, N, n, options)
%TestStat_MedianExponential_Simulations Generates M independent
%  realizations of the test statistic testStat =
%  TestStat_MedianExponential_Simulations(medians, options), the GOF test
%  statistic based on the difference of the theoretical characteristic
%  functions, cfM(t), and the empirical characteristic functions, cfE(t),
%  calculated from the n observed medians from the N-dimesnsional random
%  sample from the exponential distribution EXP(lambda) [default value is
%  lambda = 1]. Here, 
%
%    testStat = n*integral(abs(cfM(t)-cfE(t))^2*exp(-r*abs(t)^p),-Inf,Inf)
%
%  were r and p are the specified control parameters (as default values we
%  set r = 2, p = 2).
% 
% SYNTAX:
%  testStat = TestStat_MedianExponential_Simulations(M, N, n, options)
%
% INPUT:
%  M        - number of generated testStat. If empty, the default value
%             is M = 1000. 
%  N        - size of the random sample from uniform distribution used to
%             calculate the medians. If empty, the default value
%             is N = 10.  
%  n        - sample size of the medians used to defined the empirical CF.
%             If empty, the default value is n = 100. 
%  options  - structure with the following parameters:
%              options.M = 1000;       % number of random realizations of
%                                      % the test statistic 
%              options.N = 10;         % size of the random sample from
%                                      % uniform distribution used to
%                                      % calculate the medians
%              options.n = 100;        % number of observed medians from
%                                      % the N-dimensional random samples
%                                      % from the uniform distribution
%              options.lambda = 1;     % parameter of the exponential
%                                      % distribution
%              options.r = 2;          % control parameter r 
%              options.p = 2;          % control parameter p
%              options.Upp = 10;       % upper integration limit
%              options.cf_MedianExponential = @(t) cf_MedianExponential(t, lambda, N);
%                                      % function handle of the
%                                      % characteristic function of the
%                                      % sample median from the sample of
%                                      % size N from the exponential
%                                      % distribution EXP(lambda) 
%
% %EXAMPLE
%  M = 1000;
%  N = 21;
%  n = 50;
%  testStat = sort(TestStat_MedianExponential_Simulations(M, N, n));
%  prob = [0.9 0.95 0.975 0.99];
%  quantiles  = quantile(testStat,prob)
%  ecdf(testStat)
%  grid
%  xlabel('testStat')
%  ylabel('cdf')
%  title('Empirical CDF of the testStat')
%
% %EXAMPLE 
%  M = [];
%  N = [];
%  n = [];
%  clear options
%  options.M = 10000;
%  options.N = 20;
%  options.n = 50;
%  options.lambda = 2;
%  options.r = 2;
%  options.p = 2;
%  testStat = sort(TestStat_MedianExponential_Simulations(M, N, n, options));
%  prob = [0.9 0.95 0.975 0.99];
%  quantiles  = quantile(testStat,prob)
%  ecdf(testStat)
%  grid
%  xlabel('testStat')
%  ylabel('cdf')
%  title('Empirical CDF of the testStat')

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Jan-2024 13:40:13

%% CHECK/SET THE INPUT PARAMETERS
narginchk(0, 4);
if nargin < 4, options = []; end
if nargin < 3, n = []; end
if nargin < 2, N = []; end
if nargin < 1, M = []; end

if ~isfield(options, 'M')
    options.M = 1000;
end

if ~isempty(M)
    options.M = M;
end

if ~isfield(options, 'N')
    options.N = 10;
end

if ~isempty(N)
    options.N = N;
end

if ~isfield(options, 'n')
    options.n = 100;
end

if ~isempty(n)
    options.n = n;
end

if ~isfield(options, 'lambda')
    options.lambda = 1;
end

if ~isfield(options, 'r')
    options.r = 2;
end

if ~isfield(options, 'p')
    options.p = 2;
end

if ~isfield(options, 'Upp')
    options.Upp = 10;
end

if ~isfield(options, 'cf_MedianExponential')
    options.cf_MedianExponential = @(t) cf_MedianExponential(t,  options.lambda , options.N );
end

%% Algorithm

testStat = zeros(options.M,1);

for i = 1:options.M
    R = exprnd(1/options.lambda,options.n,options.N); % Generated random values from Exponential 
    Medians = median(R,2);
    testStat(i) = TestStat_MedianExponential(Medians, options);
end

end