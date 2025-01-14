function cf = cfND_Empirical(t, data)
%% cfND_Empirical
%  Computes the characteristic function (CF) of the N-dimensional empirical 
%  distribution specified by the data matrix `data`. The CF is equivalent 
%  to a weighted N-dimensional mixture distribution with equal weights.
%
%  For a given data matrix `data` and evaluation points `t = [t1, ..., tN]`, 
%  the CF is defined as:
%    CF(t) = (1/n) * sum_{j=1}^n exp(1i * (d_j1 * t1 + ... + d_jN * tN)),
%  where `data = [d1; ...; dn]` is the matrix of empirical observations.
%
% SYNTAX:
%  cf = cfND_Empirical(t, data)
%
% INPUTS:
%  t      - (M x N)-matrix where the CF is evaluated. Each row represents 
%           an evaluation point in N dimensions. If `t` is an M-vector, it 
%           is assumed that `t1 = t, ..., tN = t`.
%  data   - (n x N)-matrix of empirical observations `data = [d1; ...; dn]`. 
%           Each row corresponds to an N-dimensional observation.
%
% OUTPUT:
%  cf     - (M x 1)-vector of CF values evaluated at the specified points in `t`.
%
% REFERENCES:
% [1] https://en.wikipedia.org/wiki/Empirical_distribution_function
%
% EXAMPLE 1:
% % CF of the Bivariate Empirical Distributions specified by the data
%   n       = 100;
%   N       = 2;
%   data    = rand(n,N);
%   cf      = @(t) cfND_Empirical(t,data);
%   tt      = linspace(-4,4,51);
%   [t1 t2] = meshgrid(tt,tt);
%   t = [t1(:) t2(:)];
%   figure; 
%   contour(t1,t2,real(reshape(cf(t),51,51)),'ShowText','on')
%   title('Contour Plot of the Bivariate Mixture Distribution')
%
% EXAMPLE 2:
% % CF of the weighted 3-variate Mixture Distribution
%   n       = 100;
%   N       = 3;
%   data    = rand(n,N);
%   cf      = @(t) cfND_Empirical(t,data);
%   tt      = linspace(-4,4,21);
%   [t1 t2 t3] = meshgrid(tt);
%   t       = [t1(:) t2(:) t3(:)];
%   disp([t cf(t)])
%
% EXAMPLE 3:
% % PDF/CDF of weighted Bivariate Mixture Distribution 
% % smoothed by Gaussian kernel 
%   n       = 100;
%   N       = 2;
%   data    = rand(n,N);
%   weight  = 1/n;
%   bandwidth = 0.1;
%   cf = @(t) cfND_DiracMixture(t,data,weight)  .* cf2D_Normal(bandwidth*t);
%   clear options;
%   options.isInterp = true;
%   result = cf2Dist2D(cf,[],options)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-May-2024 16:29:07'
% Updated: '13-Jan-2025 23:27:34'

%% ALGORITHM
% cf = cfND_Empirical(t,data)

%% Input Validation
narginchk(2, 2);

% Check if data is provided and valid
if isempty(data)
    error('InputError: `data` must be a non-empty matrix of empirical observations.');
end

% Dimensions of `data`
[n, N] = size(data);

% Validate `t`
if size(t, 2) ~= N && size(t, 1) ~= N
    error('DimensionMismatch: `t` must have the same number of columns as the dimension of `data`.');
end

%% Compute the Characteristic Function
% Set uniform weights (1/n for each data point)
weight = ones(n, 1) / n;

% Use the `cfND_DiracMixture` function to compute the CF
cf = cfND_DiracMixture(t, data, weight);

end