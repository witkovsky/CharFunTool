function cf = cf2D_Empirical(t,data)
%% cf2D_Empirical  
%  Characteristic function of the Bivariate Empirical Distribution
%  specified by the N x 2 data matrix, D = [d_11,d_12;...;d_N1,d_N2].  
%  
%  For given data and t = [t1,t2] the characteristic function is given as a
%  weighted Bivariate Mixture Distribution ith equal weighs, 
%  cf(t) = (1/N)*exp(1i*(d_11*t1+d_12*t2)) +...+
%  (1/N)*exp(1i*(d_N1*t1+d_N2*t2)), where exp(1i*(d_1*t1+d_2*t2)).
%
%  
% SYNTAX:
%  cf = cf2D_Empirical(t,data)
%
% INPUTS:
%  t      - (n x 2)-matrix t = [t1 t2] or a cell t = {t1 t2} of real
%           values, where the CF is evaluated. If t is cell of two vectors,
%           then cf is a (n1 x n2)-matrix where n1 = length(t1), n2 =
%           length(t2). 
%  data   - (n x 2)-matrix data = [d1 d2] or a cell data = {d1 d2} of real
%           values. If data is cell of two vectors d1 and d2, then the new
%           data matrix is created from the meshgrid, i.e. [d2,d1] =
%           meshgrid(data{2},data{1}) and hence new data = [d1(:),d2(:)]. 
%  weight - vector of weights of the distribution mixture. If empty,
%           default value is weight = 1/length(d1). 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
%
% EXAMPLE 1:
% % CF of the Bivariate Empirical Distributions specified by the data
%   N       = 100;
%   data    = rand(N,2);
%   cf      = @(t) cf2D_Empirical(t,data);
%   tt      = linspace(-4,4,51);
%   [t1 t2] = meshgrid(tt,tt);
%   t = [t1(:) t2(:)];
%   figure; 
%   contour(t1,t2,real(reshape(cf(t),51,51)),'ShowText','on')
%   title('Contour Plot of the Bivariate Mixture Distribution')
%
% EXAMPLE 2:
% % PDF/CDF of the Bivariate Empirical Distributions smoothed by 
% % Gaussian kernel  
%   N         = 100;
%   data      = rand(N,2);
%   bandwidth = 0.1;
%   cf        = @(t) cf2D_Empirical(t,data) .* cf2D_Normal(bandwidth*t);
%   clear options;
%   options.isInterp = true;
%   result = cf2Dist2D(cf,[],options)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Jan-2023 21:43:54

%% ALGORITHM
% cf = cf2D_Empirical(t,data)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, data = []; end

weight = [];
cf = cf2D_DiracMixture(t,data,weight);

end