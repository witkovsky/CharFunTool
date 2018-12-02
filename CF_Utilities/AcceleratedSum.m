function [aSum,Sum] = AcceleratedSum(data,method)
%AcceleratedSum Estimates the limit of infinite sum of alternating series,
%  based on methods for fast convergence acceleration of alternating
%  series. If data is a vector of nonnegative real numbers,
%  AcceleratedSum(data) estimates the infinite sum of the alternating
%  series values: aSum = data(1) - data(2) + data(3) - data(4) + ... . 
%
%  AcceleratedSum(data,method) where data is a vector of length N and
%  method is a  string selects one of the following linear acceleration
%  methods: 
%   'euler', Euler's method with error ~ 2^-N
%   'binom', a binomial weighting with error ~ 3^-N
%   'cheby', a weighting based on Chebyshev polynomials, error ~ 5.828^-N
%  Default method is 'cheby'.
%
% SYNTAX:
%  [aSum,Sum] = AcceleratedSum(data,method)
%
% EXAMPLE 1	(Estimate of log(2))
%  % Here we consider estimate of log(2) = 1 - 1/2 + 1/3 - 1/4 + ...
%  N     = 15;
%  data  = 1./(1:N);
%  [aSum,Sum] = AcceleratedSum(data);
%  true  = log(2);
%  error = abs(true-aSum);
%  disp([true;aSum;Sum;error])
%
% EXAMPLE 2	(Estimate of pi)
%  % Here we consider estimate of pi = 4 - 4/3 + 4/5 - 4/7 + ...
%  N     = 15;
%  data  = 4./(1:2:2*N);
%  [aSum,Sum] = AcceleratedSum(data);
%  true  = pi;
%  error = abs(true-aSum);
%  disp([true;aSum;Sum;error])
%
% ACKNOWLEDGEMENS:
%  The algorithm AcceleratedSum is adapted version of the algorithm altsum
%  by Jonas Lundgren (2009), <splinefit@gmail.com>, see Matlab Central File
%  Exchange:
%  https://www.mathworks.com/matlabcentral/fileexchange/25200-altsum
%
% REFERENCES
%  [1] Weniger, E.J., 1989. Nonlinear sequence transformations for the
%      acceleration of convergence and the summation of divergent series.
%      Computer Physics Reports, 10(5-6), pp.189-371.  
%  [2] Cohen, H., Villegas, F.R. and Zagier, D., 2000. Convergence
%      acceleration of alternating series. Experimental mathematics, 9(1),
%      pp.3-12.  

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-Dec-2018 10:36:30

%% FUNCTION
%  s = AcceleratedSum(data,method)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, method = []; end

if isempty(method)
    method = 'cheby';
end

%% ALGORITHM
% Index vector
n = numel(data);
k = (1:n)';

% Select method
switch method
    case 'euler'
        % Euler's method ~ 2^-n
        p = 1/2;
        Q = (n-k+1)./k;
    case 'binom'
        % Binomial weighting ~ 3^-n
        p = 2/3;
        Q = (n-k+1)./k/2;
    case 'cheby'
        % Chebyshev weighting ~ 5.828^-n
        p = 4/(3+sqrt(8));
        Q = (n-k+1)./k.*(n-k+1/2)./(2*n-k);
    otherwise
        p = 4/(3+sqrt(8));
        Q = (n-k+1)./k.*(n-k+1/2)./(2*n-k);
end

% d = cumprod([p^n; Q]);
d = exp(cumsum([n*log(p); log(Q)]));        
w = cumsum(d);

% Alternate signs
data = data(:);
data(2:2:n) = -data(2:2:n);

% Reverse order
data = data(n:-1:1);

% Weighted sum
aSum = sum(w(1:n).*data)/w(n+1);

% Standard sum for comparisons
if nargout > 1
    Sum = sum(data);
end
end