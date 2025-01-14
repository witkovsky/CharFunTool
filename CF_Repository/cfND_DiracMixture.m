function cf = cfND_DiracMixture(t, data, weight)
%% cfND_DiracMixture
%  Computes the characteristic function (CF) of a weighted N-dimensional 
%  mixture distribution formed by a mixture of `n` independent N-dimensional 
%  Dirac random vectors. Each vector is concentrated at a fixed N-dimensional 
%  point given by rows of the `data` matrix.
%
%  For a given data matrix `data` and evaluation points `t = [t1, ..., tN]`, 
%  the CF is defined as:
%    CF(t) = sum_{j=1}^n weight_j * exp(1i * (d_j1 * t1 + ... + d_jN * tN)),
%  where `data = [d1; ...; dn]` is the matrix of fixed points, `weight` is 
%  the vector of weights for the mixture, and `exp(1i * (...))` is the CF of 
%  the N-dimensional Dirac random variable concentrated at the vector `d_j`.
%
% SYNTAX:
%  cf = cfND_DiracMixture(t, data, weight)
%
% INPUTS:
%  t      - (M x N)-matrix where the CF is evaluated. Each row represents 
%           an evaluation point in N dimensions. If `t` is an M-vector, it 
%           is assumed that `t1 = t, ..., tN = t`.
%  data   - (n x N)-matrix of fixed points `data = [d1; ...; dn]`. Each row 
%           corresponds to the location of a Dirac random vector in N dimensions.
%  weight - Vector of weights of the distribution mixture. If empty, 
%           the default value is `weight = 1/n` (uniform weights).
%
% OUTPUT:
%  cf     - (M x 1)-vector of CF values evaluated at the specified points in `t`.
%
% REFERENCES:
% [1] https://en.wikipedia.org/wiki/Empirical_distribution_function
%
% EXAMPLE 1:
% % CF of the weighted n-dimensional Mixture Distribution
%   n       = 100;
%   N       = 3;
%   data    = rand(n,N);
%   weight  = 1/n;
%   cf      = @(t) cfND_DiracMixture(t,data,weight);
%   tt      = linspace(-4,4,21);
%   [t1 t2 t3] = meshgrid(tt);
%   t       = [t1(:) t2(:) t3(:)];
%   disp([t cf(t)])
%
% EXAMPLE 2:
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
% cf = cfND_DiracMixture(t,data,weight)

%% Input Validation and Default Parameters
narginchk(2, 3);

% Dimensions of data
[n, N] = size(data);

% Default value for weight
if nargin < 3 || isempty(weight)
    weight = ones(n, 1) / n;
end

%% Dimension Check for Consistency
% Extract the first column of data for consistency checks
d1 = data(:, 1);

% Use distchck to verify that `data` and `weight` have compatible sizes
[err, d1, weight] = distchck(2, d1, weight);
if err > 0
    error('InputSizeMismatch: Dimensions of `data` and `weight` are inconsistent.');
end

%% Validate input dimensions
if numel(weight) ~= n
    error('Dimension mismatch: `weight` must have the same length as the number of rows in `data`.');
end
if size(t, 2) ~= N && size(t, 1) ~= N
    error('Dimension mismatch: `t` must have the same number of columns as the dimension of `data`.');
end

%% Reshape and Validate the Argument t

% Get the original size of `t`
size_t = size(t);
size_t_Min = min(size_t(1), size_t(2));
size_t_Max = max(size_t(1), size_t(2));

switch size_t_Min
    case 1 % If `t` is a vector, assume t1 = t, ..., tN = t
        % Create a new `t` matrix where all columns are the same as the vector
        if size_t_Max > N
            % If the vector length exceeds N, broadcast it to match N dimensions
            size_reshape = size(t);
            t = t(:); % Ensure `t` is a column vector
            t = t * ones(1, N); % Replicate the vector across N columns
        elseif size_t_Max == N
            % If the vector length matches N, treat it as a row vector
            t = t(:)'; % Ensure `t` is a row vector
            size_reshape = [1, 1]; % Save the output shape as single-row
        else
            error('InputSizeMismatch: The size of `t` does not match the required dimensions.');
        end
        
    case N % If `t` has N rows, transpose it to have N columns
        if size_t_Max > N && size_t(1) == N
            t = t'; % Transpose `t` to have N columns
            size_reshape = [1, size_t(2)]; % Save the output shape
        else
            size_reshape = [size_t(1), 1]; % Default shape if no transposition is needed
        end
        
    otherwise
        % If `t` does not match expected dimensions, raise an error
        error('InputSizeMismatch: `t` must be a vector or a matrix with N columns or N rows.');
end

%% Compute the Characteristic Function
% Compute the argument of the exponential
arg = t(:, 1) * d1';
for i = 2:N
    arg = arg + t(:, i) * data(:, i)';
end

% Compute the weighted sum of exponential terms
cf = exp(1i * arg) * weight;

% Reshape the result for output
cf = reshape(cf, size_reshape);

end