function cf = cfND_Empirical(t,data)
%% cfND_Empirical  
%  Characteristic function of the N-VARIATE EMPIRICAL DISTRIBUTION
%  specified by the (n x N) data matrix, data = [d_11,...,d_1N; ... ;
%  d_n1,...,d_nN].
%  
%  For given data and t = [t1,...,tN] the characteristic function is given
%  as a weighted N-variate Mixture Distribution with equal weights, given
%  as cf(t) = (1/n)*exp(1i*(d_11*t1 +...+ d_1N*tN)) +...+
%  (1/n)*exp(1i*(d_n1*t1 +...+ d_nN*tN)), where exp(1i*(d_1*t1 +...+
%  d_N*tN)) represents the characteristic function of the N-dimensional
%  DIRAC RV concentrated at the vector d = [d_1,...,d_N].
%
%  
% SYNTAX:
%  cf = cfND_Empirical(t,data)
%
% INPUTS:
%  t      - (M x N)-matrix t = [t1,...,tN] of real values, where the CF is
%           evaluated.   
%  data   - (n x N)-matrix data = [d1,...,dN] of real values.
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
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

%% ALGORITHM
% cf = cfND_Empirical(t,data)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, data = []; end

weight = [];
cf = cfND_DiracMixture(t,data,weight);

end