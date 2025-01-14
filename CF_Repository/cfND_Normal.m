function cf = cfND_Normal(t, mu, Sigma)
%% cfND_Normal
%  Computes the characteristic function (CF) of an N-variate normal 
%  (Gaussian) random variable, X = (X1, ..., XN), with mean vector `mu` 
%  and covariance matrix `Sigma`.
%
%  The CF of a multivariate normal distribution is given by:
%    CF(t) = exp(1i * t * mu') * exp(-0.5 * t * Sigma * t')
%
%  Reference:
%  https://en.wikipedia.org/wiki/Multivariate_normal_distribution
%
% SYNTAX:
%  cf = cfND_Normal(t, mu, Sigma)
%
% INPUTS:
%  t     - (M x N) matrix of real values where the CF is evaluated. 
%          Each row corresponds to an evaluation point in N dimensions.
%  mu    - (N x 1) vector of means for the normal distribution. 
%          Default: [0; ...; 0] (N-dimensional zero vector).
%  Sigma - (N x N) covariance matrix. Must be symmetric positive semi-definite.
%          Default: eye(N) (identity matrix).
%
% OUTPUT:
%  cf    - (M x 1) vector of CF values at the specified points in `t`.
%
% EXAMPLE 1:
% % CF of the Bivariate Normal RV
%   mu    = [0 0];
%   Sigma = [1 0.5; 0.5 1];
%   cf    = @(t) cfND_Normal(t,mu,Sigma);
%   tt    = linspace(-4,4,51);
%   [t1 t2] = meshgrid(tt,tt);
%   t = [t1(:) t2(:)];
%   figure; 
%   contour(t1,t2,real(reshape(cf(t),51,51)),'ShowText','on')
%   title('Contour Plot of the CF of the Bivariate Normal Distribution')
%
% EXAMPLE 2:
% % CF of the 3-variate Normal RV
%   mu    = [1 2 3]';
%   Sigma = [1 0.5 0.1; 0.5 1 0.3; 0.1 0.3 2];
%   cf    = @(t) cfND_Normal(t,mu,Sigma);
%   tt    = linspace(-4,4,11);
%   [t1 t2 t3] = meshgrid(tt);
%   t = [t1(:) t2(:) t3(:)];
%   disp([t cf(t)])
%
% EXAMPLE 3:
% % CF of convolution of 3 3-variate Normal Distributions
%   mu1    = [1 2 3]';
%   Sigma1 = [1 0.5 0.1; 0.5 1 0.3; 0.1 0.3 2];
%   Coef1  = [0.03 0.47 0.64; 0.27 0.51 0.24; 0.29 0.65 0.89];
%   mu2    = [-1 0 1]';
%   Sigma2 = [3 0.25 0.0; 0.25 1 0.3; 0.0 0.3 1];
%   Coef2  = [0.86 0.21 0.39; 0.88 0.25 0.96; 0.61 0.16 0.82];
%   mu3    = [0 0 0]';
%   Sigma3 = [1 0 0; 0 1 0; 0 0 1];
%   Coef3  = [1 0 0; 0 1 0; 0 0 1];
%   cf1    = @(t) cfND_Normal(t*Coef1,mu1,Sigma1);
%   cf2    = @(t) cfND_Normal(t*Coef2,mu2,Sigma2);
%   cf3    = @(t) cfND_Normal(t*Coef3,mu3,Sigma3);
%   cf     = @(t) cf1(t).*cf2(t).*cf3(t);
%   tt    = linspace(-4,4,11);
%   [t1 t2 t3] = meshgrid(tt);
%   t = [t1(:) t2(:) t3(:)];
%   disp([t cf(t)])
%
% EXAMPLE 4:
% % CF of convolution of 3 bivariate Normal Distributions
%   mu1    = [1 2]';
%   Sigma1 = [1 0.5; 0.5 1];
%   Coef1  = [0.03 0.47; 0.27 0.51];
%   mu2    = [-1 0 ]';
%   Sigma2 = [3 0.25; 0.25 1];
%   Coef2  = [0.86 0.21; 0.88 0.25];
%   mu3    = [0 0]';
%   Sigma3 = [1 0; 0 1];
%   Coef3  = [1 0; 0 1];
%   cf1    = @(t) cfND_Normal(t*Coef1,mu1,Sigma1);
%   cf2    = @(t) cfND_Normal(t*Coef2,mu2,Sigma2);
%   cf3    = @(t) cfND_Normal(t*Coef3,mu3,Sigma3);
%   cf     = @(t) cf1(t).*cf2(t).*cf3(t);
%   clear options;
%   options.isInterp = true;
%   options.chebyPts = 101;
%   result = cf2Dist2D(cf,[],options) 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: '01-May-2024 16:29:07'
% Updated: '13-Jan-2025 23:27:34'

%% Input Validation and Default Values
narginchk(1, 3);

% Get dimensions of `t`
N = size(t,2);

% Set default values for `mu` and `Sigma`
if nargin < 3 || isempty(Sigma), Sigma = eye(N); end
if nargin < 2 || isempty(mu), mu = zeros(N, 1); end

% Validate input dimensions
if numel(mu) ~= N
    error('Dimension mismatch: `mu` must have length equal to the number of columns in `t`.');
end
if ~isequal(size(Sigma), [N, N])
    error('Dimension mismatch: `Sigma` must be an (N x N) matrix.');
end
if ~issymmetric(Sigma) || any(eig(Sigma) < 0)
    error('Invalid covariance matrix: `Sigma` must be symmetric positive semi-definite.');
end

%% Compute the Characteristic Function
% 1. Linear phase term: exp(1i * t * mu')
linearTerm = 1i * (t * mu(:));

% 2. Quadratic term: exp(-0.5 * t * Sigma * t')
quadraticTerm = -0.5 * sum((t * Sigma) .* t, 2);

% 3. Combine terms
cf = exp(linearTerm + quadraticTerm);

end