function cf = cfS_Laplace(t,beta,coef,niid)
%% cfS_Laplace
%  Characteristic function of a symmetric LAPLACE distribution with scale
%  parameter beta > 0. 
%
%  cfS_Laplace evaluates the characteristic function of a linear
%  combination of independent (symmetric) LAPLACE sistributed random
%  variables. 
%
%  That is, cfS_Laplace evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Laplace (0,beta_i) are
%  inedependent zero-mean RVs with the scale parameters beta_i > 0, for i =
%  1,...,N.  
%
%  The characteristic function of Y is defined by 
%   cf(t) = Prod( 1 ./ (1 + (t*coef(i)*beta(i)).^2 )
%
% SYNTAX
%  cf = cfS_Laplace(t)
%  cf = cfS_Laplace(t,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  beta  - vector of the scale parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          LAPLACE random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.  
% 
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Normal_distribution.
%
%
% EXAMPLE 1:
% % CF of the the linear combination of symmetric Laplace RV
%   beta = 1;
%   t    = linspace(-10,10,201);
%   cf   = cfS_Laplace(t,beta);
%   figure; plot(t,cf);grid on
%   title('Characteristic function of the symmetric Laplace RV')
%
% EXAMPLE 2:
% % PDF/CDF of the symmetric Laplace RVs
%   beta = 1;
%   x    = linspace(-5,5,101);
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cfS_Laplace(t,beta);
%   result = cf2DistGP(cf,x,prob);
%
% EXAMPLE 3:
% % PDF/CDF of the linear combination of symmetric Laplace RVs 
%   beta = [0.1 0.2 0.3 0.4];
%   coef = [1 2 3 4];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cfS_Laplace(t,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);

% (c) 2018 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Aug-2018 16:00:43
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
narginchk(1, 4);
if nargin < 4, niid = []; end
if nargin < 3, coef = []; end
if nargin < 2, beta = []; end

 cf = cf_Laplace(t,[],beta,coef,niid);

end