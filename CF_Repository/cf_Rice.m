function cf = cf_Rice(t,scale,delta,coef,niid)
%cf_Rice 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent non-central Rice distributed random variables, with the
%  scale parameters scale_i > 0, and the non-centrality parameters delta_i
%  >= 0, for i  = 1,...,N.
%  
%  The non-central Rice distribution is a continuous probability
%  distribution for positive-valued random variables. It is a special case
%  of the non-central chi distribution in three degrees of freedom. The
%  Rice distribution is also known as the non-central distribution Rayleigh
%  distribution.
%
%  cf_Rice evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i~ RiceNC(scale_i,delta_i) are
%  inedependent non-central Rice distributed RVs with the scale
%  parameters scale_i > 0,  and the non-centrality parameters delta_i >= 0,
%  for i  = 1,...,N.
%
%  The characteristic function of X ~ RiceNC(scale,delta) is defined by
%   cf_Rice(t) = cf_ChiNC(scale*t,df=2,delta/scale), 
%  where by cf_ChiNC(t,df,delta) we denote the characteristic function of
%  the noncentral chi distribution with df degrees of freedom and the
%  non-centrality parameter delta. Hence, the characteristic function of Y
%  is  
%   cf(t) = Prod ( cf_Rice(t,scale_i,delta_i) )
%
% SYNTAX:
%  cf = cf_Rice(t,scale,delta,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  scale - vector of the scale parameters of the Rice distributed
%          random variables.  If scale is scalar, it is assumed that all
%          scale parameters are equal. If empty, default value is scale =
%          1.
%  delta - vector of the non-centrality parameters delta >= 0. If empty,
%          default value is delta = 0. Notice that the noncentrality
%          parameter delta can be interpreted as a square root of the sum
%          of squared means, delta = sqrt(sum_{i=1}^df mu_i^2).
%  coef  - vector of the coefficients of the linear combination of the
%          Rice distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Noncentral_chi_distribution
%  https://en.wikipedia.org/wiki/Rice_distribution   
%
% EXAMPLE 1:
% % CF of the distribution of Rice RV with scale = 3
%   scale = 1;
%   delta = 5;
%   t     = linspace(-5,5,501);
%   cf    = cf_Rice(t,scale,delta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of the non-central Rice RV')
%
% EXAMPLE 2:
% % PDF/CDF of the distribution of Rice RV with scale = 3
%   scale = 1;
%   delta = 5;
%   cf    = @(t) cf_Rice(t,scale,delta);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,10,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% EXAMPLE 3: 
% % CF of a linear combination of independent Rice RVs
%   scale = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   t     = linspace(-2,2,501);
%   cf    = cf_Rice(t,scale,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of independent Rice RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a linear combination of independent Rice RVs
%   scale = [1 2 3];
%   delta = [1 1 1];
%   coef  = [1 1 1];
%   cf    = @(t) cf_Rice(t,scale,delta,coef);
%   clear options
%   options.N = 2^10;
%   options.xMin = 0;
%   x = linspace(0,20,201);
%   prob = [0.9 0.95 0.975 0.99];
%   result = cf2DistGP(cf,x,prob,options);
%
% See also: cf_ChiNC, cf_RayleighNC

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 5-Oct-2018 17:08:51

%% ALGORITHM
%  cf_Rice is an alias name for cf_RayleighNC
%  cf = cf_RayleighNC(t,scale,delta,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid  = []; end
if nargin < 4, coef  = []; end
if nargin < 3, delta = []; end
if nargin < 2, scale = []; end

cf = cf_RayleighNC(t,scale,delta,coef,niid);

end