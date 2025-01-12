function testStat = TestStat_MedianUniform_Simulations(M, N, n, options)
%TestStat_MedianUniform_Simulations Generates M independent realizations of
%  the test statistic testStat =
%  TestStat_MedianUniform_Simulations(medians, options), the GOF test
%  statistic based on the difference of the theoretical characteristic
%  functions, cfM(t), and the empirical characteristic functions, cfE(t),
%  calculated from the n observed medians from the N-dimensional random
%  sample from the uniform distribution U(-1,1). Here,
%
%    testStat = n*integral(abs(cfM(t)-cfE(t))^2*exp(-r*abs(t)^p),-Inf,Inf)
%
%  were r and p are the specified control parameters (as default values we
%  set r = 2, p = 2).
% 
% SYNTAX:
%  testStat = TestStat_MedianUniform_Simulations(M, N, n, options)
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
%              options.a = -1;         % lower bound of the uniform
%                                      % distribution
%              options.b = 1;          % upper bound of the uniform
%                                      % distribution
%              options.r = 2;          % control parameter r 
%              options.p = 2;          % control parameter p
%              options.Upp = 10;       % upper integration limit
%              options.cf_MedianUniform = @(t) cf_MedianUniform(t, -1, 1, N);
%                                      % function handle of the
%                                      % characteristic function of the
%                                      % sample median from the sample of
%                                      % size N from the uniform
%                                      % distribution U(-1,1) 
%
% %EXAMPLE
%  M = 1000;
%  N = 21;
%  n = 50;
%  testStat = sort(TestStat_MedianUniform_Simulations(M, N, n));
%  prob = [0.9 0.95 0.975 0.99];
%  quantiles  = quantile(testStat,prob)
%  ecdf(testStat)
%  grid
%  xlabel('testStat')
%  ylabel('cdf')
%  title('Empirical CDF of the testStat')
%
% %EXAMPLE (Requires LONG TIME CALCULATIONS if N is even and large)
%  M = [];
%  N = [];
%  n = [];
%  clear options
%  options.M = 500;
%  options.N = 4;
%  options.n = 50;
%  testStat = sort(TestStat_MedianUniform_Simulations(M, N, n, options));
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

if ~isfield(options, 'Upp')
    options.Upp = 10;
end

if ~isfield(options, 'cf_MedianUniform')
    options.cf_MedianUniform = @(t) cf_MedianUniform(t,  options.a ,  options.b , options.N );
end

%% Algorithm

testStat = zeros(options.M,1);

for i = 1:options.M
    % Generated random values from U(a,b) 
    R = rand(options.n,options.N)*(options.b - options.a) + options.a; 
    Medians = median(R,2);
    testStat(i) = TestStat_MedianUniform(Medians, options);
end

end