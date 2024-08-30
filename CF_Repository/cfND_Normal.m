function cf = cfND_Normal(t,mu,Sigma)
%% cfND_Normal  
%  Characteristic function of N-VARIATE NORMAL (GAUSSIAN) random variable,
%  say X = (X1,...,Xn), with location parameters specified by (N x 1) mean
%  vector mu = (mu1,...,muN)' (real), and scale parameters specified by the
%  (N x N) covariance matrix Sigma.
%
%  Here, we consider the multivariate normal distribution and its
%  characteristic function as specified at WIKIPEDIA
%  https://en.wikipedia.org/wiki/Multivariate_normal_distribution
%
%  
% SYNTAX:
%  cf = cfND_Normal(t,mu,Sigma)
%
% INPUTS:
%  t      - (M x N)-matrix t = [t1,...,tN] of real values, where the CF is
%           evaluated.  
%  mu     - N-dimensional vector mu = [mu1,...,muN]'; of real location
%           parameters. If empty, default value is mu = [0,...,0]'.   
%  Sigma  - (N x N) covariance matrix Sigma. If empty, default value
%           is Sigma = eye(N).   
%
% WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
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
% Ver.: 01-May-2024 16:29:07'

%% ALGORITHM
% cf = cfND_Normal(t,mu,Sigma)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, Sigma = []; end
if nargin < 2, mu = []; end

%% Characteristic function

if ismatrix(t)
    [~,N] = size(t);
else
    error('InputSizeMismatch');
end

if isempty(mu), mu = zeros(N,1); end
if isempty(Sigma), Sigma = eye(N); end

cf = exp(1i*t*mu(:)) .* exp(-sum((t*Sigma).*t,2)/2);

end