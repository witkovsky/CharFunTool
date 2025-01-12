function testStat = TestStat_Simulations(M,N,n,options)
%TestStat_Simulations Generates M independent realizations of the test
% statistic testStat = TestStat_MedianUniform(medians, options), the GOF
% test statistic based on the difference of the theoretical characteristic
% functions, cfM(t), and the empirical characteristic functions, cfE(t),
% calculated from the n observed medians from the N-dimesnsional random
% sample from the uniform distribution U(-1,1). Here,
%    testStat = 2*integral(n*abs(cfM(t)-cfE(t))^2*exp(-r*abs(t)^p),0,Upp)
%  were r and p are the specified control parameters (as default values we
%  set r = 2, p = 4).
% 
% SYNTAX:
%  testStat = TestStat_Simulations(M,N,n,options)
%
% INPUT:
%  M        - number of generated testStat
%  N        - size of the random sample from uniform distribution used to
%             calculate the medians 
%  n        - sample size of the medians used to defined the empirical CF
%  options  - structure with the following parameters:
%              options.a = -1;         % lower bound of the uniform
%                                      % distribution
%              options.b = 1;          % upper bound of the uniform
%                                      % distribution
%              options.r = 2;          % control parameter r 
%              options.p = 4;          % control parameter p
%              options.N = 10;         % size of the random sample from
%                                      % uniform distribution used to
%                                      % calculate the medians
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
%  N = 11;
%  n = 100;
%  options.r = 2; 
%  options.p = 4; 
%  options.Upp = 3;
%  testStat = sort(TestStat_Simulations(M,N,n,options));
%  prob = [0.9 0.95 0.975 0.99];
%  quantiles  = quantile(testStat,prob)
%  ecdf(testStat)
%  grid
%  xlabel('testStat')
%  ylabel('cdf')
%  title('Empirical CDF of the testStat')

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Dec-2023 10:49:33

%% CHECK/SET THE INPUT PARAMETERS

if isempty(M)
    M = 1000;
end

if isempty(N)
    N = 10;
end

if isempty(n)
    n = 10;
end

narginchk(3, 4);
if nargin < 4, options = []; end

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
    options.p = 4;
end

if ~isfield(options, 'N')
    options.N = N;
end

if ~isfield(options, 'Upp')
    options.Upp = 10;
end

if ~isfield(options, 'cf_MedianUniform')
    options.cf_MedianUniform = @(t) cf_MedianUniform(t,  options.a ,  options.b , options.N );
end

%% Algorithm

a = options.a;
b = options.b;

testStat = zeros(M,1);

for i = 1:M
    M = median(rand(n,N)*(b-a)+a,2);
    testStat(i) = TestStat_MedianUniform(M, options);
end

end