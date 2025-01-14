function cf = cfND_UniformSphere(t, mu, Sigma)
%% cfND_UniformSphere
%  Computes the characteristic function (CF) of an N-dimensional random
%  variable with a uniform sphere distribution. The random variable X, X =
%  (X1, ..., XN), is defined by its location parameter vector `mu` and
%  covariance matrix `Sigma`.
%
%  For a given mean vector `mu` and covariance matrix `Sigma`, the CF of the 
%  uniform sphere distribution is given by:
%    CF(t) = exp(1i * t * mu') * (sin(T) / T),
%  where:
%    - T = sqrt(pi * sum((t * Sigma) .* t, 2)).
%
% SYNTAX:
%  cf = cfND_UniformSphere(t, mu, Sigma)
%
% INPUTS:
%  t      - (M x N)-matrix where the CF is evaluated. Each row represents 
%           an evaluation point in N dimensions.
%  mu     - (N x 1)-vector of location parameters. Default: `mu = zeros(N, 1)`.
%  Sigma  - (N x N)-covariance matrix. Default: `Sigma = eye(N)`.
%
% OUTPUT:
%  cf     - (M x 1)-vector of CF values evaluated at the specified points in `t`.
%
% REFERENCES:
% [1] Balakrishnan, N., Ma, C., and Wang, R. (2015). Logistic vector random 
%     fields with logistic direct and cross covariances. Journal of 
%     Statistical Planning and Inference, 161, pp.109-118.
%
% EXAMPLE 1:
% % CF of the Bivariate UNIFORM SPHERE DISTRIBUTION
%   mu    = [0 0];
%   Sigma = [1 0.5; 0.5 1];
%   cf    = @(t) cfND_UniformSphere(t,mu,Sigma);
%   tt    = linspace(-4,4,51);
%   [t1 t2] = meshgrid(tt,tt);
%   t = [t1(:) t2(:)];
%   figure; 
%   contour(t1,t2,real(reshape(cf(t),51,51)),'ShowText','on')
%   title('Contour Plot of the CF of the Bivariate Logistic Distribution')
%
% EXAMPLE 2:
% % CF of the 3-variate UNIFORM SPHERE DISTRIBUTION
%   mu    = [1 2 3]';
%   Sigma = [1 0.5 0.1; 0.5 1 0.3; 0.1 0.3 2];
%   cf    = @(t) cfND_UniformSphere(t,mu,Sigma);
%   tt    = linspace(-4,4,11);
%   [t1 t2 t3] = meshgrid(tt);
%   t = [t1(:) t2(:) t3(:)];
%   disp([t cf(t)])
%
% EXAMPLE 3:
% % CF of convolution of 3 3-variate UNIFORM SPHERE DISTRIBUTIONS
%   mu1    = [1 2 3]';
%   Sigma1 = [1 0.5 0.1; 0.5 1 0.3; 0.1 0.3 2];
%   Coef1  = [0.03 0.47 0.64; 0.27 0.51 0.24; 0.29 0.65 0.89];
%   mu2    = [-1 0 1]';
%   Sigma2 = [3 0.25 0.0; 0.25 1 0.3; 0.0 0.3 1];
%   Coef2  = [0.86 0.21 0.39; 0.88 0.25 0.96; 0.61 0.16 0.82];
%   mu3    = [0 0 0]';
%   Sigma3 = [1 0 0; 0 1 0; 0 0 1];
%   Coef3  = [1 0 0; 0 1 0; 0 0 1];
%   cf1    = @(t) cfND_UniformSphere(t*Coef1,mu1,Sigma1);
%   cf2    = @(t) cfND_UniformSphere(t*Coef2,mu2,Sigma2);
%   cf3    = @(t) cfND_UniformSphere(t*Coef3,mu3,Sigma3);
%   cf     = @(t) cf1(t).*cf2(t).*cf3(t);
%   tt    = linspace(-4,4,11);
%   [t1 t2 t3] = meshgrid(tt);
%   t = [t1(:) t2(:) t3(:)];
%   disp([t cf(t)])
%
% EXAMPLE 4:
% % PDF/CDF of bivariate UNIFORM SPHERE DISTRIBUTION
%   mu     = [0 0]';
%   Sigma  = [1 0; 0 1];
%   cf     = @(t) cfND_UniformSphere(t,mu,Sigma);
%   clear options;
%   options.isInterp = true;
%   options.chebyPts = 101;
%   result = cf2Dist2D(cf,[],options) 
%
% EXAMPLE 5:
% % CF of convolution of 3 bivariate UNIFORM SPHERE DISTRIBUTION
%   mu1    = [1 2]';
%   Sigma1 = [1 0.5; 0.5 1];
%   Coef1  = [0.03 0.47; 0.27 0.51];
%   mu2    = [-1 0 ]';
%   Sigma2 = [3 0.25; 0.25 1];
%   Coef2  = [0.86 0.21; 0.88 0.25];
%   mu3    = [0 0]';
%   Sigma3 = [1 0; 0 1];
%   Coef3  = [1 0; 0 1];
%   cf1    = @(t) cfND_UniformSphere(t*Coef1,mu1,Sigma1);
%   cf2    = @(t) cfND_UniformSphere(t*Coef2,mu2,Sigma2);
%   cf3    = @(t) cfND_UniformSphere(t*Coef3,mu3,Sigma3);
%   cf     = @(t) cf1(t).*cf2(t).*cf3(t);
%   clear options;
%   options.isInterp = true;
%   options.chebyPts = 101;
%   result = cf2Dist2D(cf,[],options) 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-May-2024 07:09:34
% Updated: '13-Jan-2025 23:27:34'

%% ALGORITHM
% cf = cfND_UniformSphere(t,mu,Sigma)

%% Input Validation and Default Parameters
narginchk(1, 3);

% Get dimensions of `t`
if ~ismatrix(t)
    error('InputSizeMismatch: `t` must be a 2D matrix.');
end
N = size(t,2);

% Set default values for `mu` and `Sigma`
if nargin < 3 || isempty(Sigma), Sigma = eye(N); end
if nargin < 2 || isempty(mu), mu = zeros(N, 1); end

% Validate dimensions of `mu` and `Sigma`
if numel(mu) ~= N
    error('DimensionMismatch: Length of `mu` must match the number of columns in `t`.');
end
if ~isequal(size(Sigma), [N, N])
    error('DimensionMismatch: `Sigma` must be an (N x N) matrix.');
end
if ~issymmetric(Sigma) || any(eig(Sigma) < 0)
    error('InvalidInput: `Sigma` must be symmetric and positive semi-definite.');
end

%% Compute the Characteristic Function
% Compute T = sqrt(pi * sum((t * Sigma) .* t, 2))
T = sqrt(pi * sum((t * Sigma) .* t, 2));

% Compute the CF
cf = exp(1i * t * mu(:)) .* (sin(T) ./ T);

% Handle the case where T == 0 to avoid division by zero
cf(T == 0) = 1;

end