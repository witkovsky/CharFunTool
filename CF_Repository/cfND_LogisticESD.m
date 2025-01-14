function cf = cfND_LogisticESD(t, mu, Sigma)
%% cfND_LogisticESD
%  Computes the characteristic function (CF) of an N-dimensional logistic 
%  random variable with elliptically symmetric distribution (ESD). The 
%  random variable X = (X1, ..., XN) is defined by its location parameter 
%  vector `mu` and covariance matrix `Sigma`.
%
%  The characteristic function of the 'standardized' N-variate logistic
%  distribution is given by: 
%    CF(t) = exp(1i * t * μ') * (T / sinh(T)),
%  where T = sqrt(3 * sum((t * Σ) .* t, 2)).
%
%  Here, we consider the 'standardized' N-variate logistic distribution and
%  its characteristic function as defined in Balakrishnan, Ma, and Wang
%  (2015), see [1], such that E(X) = μ and Var(X) = Σ.
% 
%  In particular, an N-variate random vector X is said to be a logistic
%  random vector if it has the same distribution as U * Z + μ, where:
%    - Z ~ N(0, Σ) is an N-dimensional normal random vector,
%    - U is a Kolmogorov–Smirnov random variable,
%    - Z and U are independent,
%    - μ is an N-dimensional (non-random) vector.
%
%  For such a logistic random vector X, we have:
%    - Mean: E(X) = μ,
%    - Covariance: Cov(X, X) = E(U^2) * Σ = (π^2 / 3) * Σ.
%
%  Hence, the 'standardized' N-variate logistic random vector X is defined 
%  as:
%    X = (1 / sqrt(E(U^2))) * U * Z + μ,
%  where Z is an N-variate normal random vector with mean 0 and 
%  variance–covariance matrix Σ, U is a Kolmogorov–Smirnov random variable, 
%  Z and U are independent, and μ is an N-dimensional (non-random) vector.
%  Here, E(X) = μ and Cov(X, X) = Σ.
%
% SYNTAX:
%  cf = cfND_LogisticESD(t, mu, Sigma)
%
% INPUTS:
%  t      - (M x N)-matrix t = [t1, ..., tN] of real values, where the CF is
%           evaluated.  
%  mu     - N-dimensional vector μ = [μ1, ..., μN]'; of real location
%           parameters. If empty, default value is μ = [0, ..., 0]'.   
%  Sigma  - (N x N) covariance matrix Σ. If empty, default value
%           is Σ = eye(N).   
%
% OUTPUT:
%  cf     - (M x 1) vector of CF values at the specified points in `t`.
%
% REFERENCES:
% [1] Balakrishnan, N., Ma, C., and Wang, R. (2015). Logistic vector random
%     fields with logistic direct and cross covariances. Journal of 
%     Statistical Planning and Inference, 161, pp.109-118.
% [2] https://en.wikipedia.org/wiki/Logistic_distribution
%
% EXAMPLES:
%  % Example 1: CF of the Bivariate Logistic RV
%  mu    = [0 0];
%  Sigma = [1 0.5; 0.5 1];
%  cf    = @(t) cfND_LogisticESD(t, mu, Sigma);
%  tt    = linspace(-4, 4, 51);
%  [t1, t2] = meshgrid(tt, tt);
%  t     = [t1(:), t2(:)];
%  figure; 
%  contour(t1, t2, real(reshape(cf(t), 51, 51)), 'ShowText', 'on');
%  title('Contour Plot of the CF of the Bivariate Logistic Distribution');
%
%  % Example 2: CF of the 3-variate Logistic RV
%  mu    = [1 2 3]';
%  Sigma = [1 0.5 0.1; 0.5 1 0.3; 0.1 0.3 2];
%  cf    = @(t) cfND_LogisticESD(t, mu, Sigma);
%  tt    = linspace(-4, 4, 11);
%  [t1, t2, t3] = meshgrid(tt);
%  t     = [t1(:), t2(:), t3(:)];
%  disp([t cf(t)]);
%
%  % Example 3: CF of convolution of 3 3-variate Logistic Distributions
%  mu1    = [1 2 3]';
%  Sigma1 = [1 0.5 0.1; 0.5 1 0.3; 0.1 0.3 2];
%  Coef1  = [0.03 0.47 0.64; 0.27 0.51 0.24; 0.29 0.65 0.89];
%  mu2    = [-1 0 1]';
%  Sigma2 = [3 0.25 0.0; 0.25 1 0.3; 0.0 0.3 1];
%  Coef2  = [0.86 0.21 0.39; 0.88 0.25 0.96; 0.61 0.16 0.82];
%  mu3    = [0 0 0]';
%  Sigma3 = [1 0 0; 0 1 0; 0 0 1];
%  Coef3  = [1 0 0; 0 1 0; 0 0 1];
%  cf1    = @(t) cfND_LogisticESD(t * Coef1, mu1, Sigma1);
%  cf2    = @(t) cfND_LogisticESD(t * Coef2, mu2, Sigma2);
%  cf3    = @(t) cfND_LogisticESD(t * Coef3, mu3, Sigma3);
%  cf     = @(t) cf1(t) .* cf2(t) .* cf3(t);
%  tt    = linspace(-4, 4, 11);
%  [t1, t2, t3] = meshgrid(tt);
%  t     = [t1(:), t2(:), t3(:)];
%  disp([t cf(t)]);
%
%  % Example 4: CF of convolution of 3 bivariate Logistic Distributions
%  mu1    = [1 2]';
%  Sigma1 = [1 0.5; 0.5 1];
%  Coef1  = [0.03 0.47; 0.27 0.51];
%  mu2    = [-1 0 ]';
%  Sigma2 = [3 0.25; 0.25 1];
%  Coef2  = [0.86 0.21; 0.88 0.25];
%  mu3    = [0 0]';
%  Sigma3 = [1 0; 0 1];
%  Coef3  = [1 0; 0 1];
%  cf1    = @(t) cfND_LogisticESD(t * Coef1, mu1, Sigma1);
%  cf2    = @(t) cfND_LogisticESD(t * Coef2, mu2, Sigma2);
%  cf3    = @(t) cfND_LogisticESD(t * Coef3, mu3, Sigma3);
%  cf     = @(t) cf1(t) .* cf2(t) .* cf3(t);
%  clear options;
%  options.isInterp = true;
%  options.chebyPts = 101;
%  result = cf2Dist2D(cf, [], options);

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-May-2024 16:29:07'
% Updated: '13-Jan-2025 23:27:34'

%% ALGORITHM
% cf = cfND_LogisticESD(t,mu,Sigma)

%% Input Validation and Default Parameters
narginchk(1, 3);
[~, N] = size(t);
if nargin < 3 || isempty(Sigma), Sigma = eye(N); end
if nargin < 2 || isempty(mu), mu = zeros(N, 1); end
if numel(mu) ~= N
    error('Dimension mismatch: `mu` must match the number of columns in `t`.');
end
if ~isequal(size(Sigma), [N, N]) || ~issymmetric(Sigma) || any(eig(Sigma) < 0)
    error('Invalid covariance matrix: `Sigma` must be symmetric positive semi-definite.');
end

%% Characteristic Function Computation
T = sqrt(3 * sum((t * Sigma) .* t, 2));
cf = exp(1i * t * mu(:)) .* (T ./ sinh(T));
cf(T == 0) = 1; % Handle T = 0 to avoid NaN

end