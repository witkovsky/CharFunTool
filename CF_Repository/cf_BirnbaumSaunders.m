function cf = cf_BirnbaumSaunders(t,alpha,beta,coef,niid)
%% cf_BirnbaumSaunders 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent random variables with Birnbaum–Saunders distributions.
%
%  The Birnbaum–Saunders distribution is known as the fatigue life
%  distribution. This probability distribution is used extensively in
%  reliability applications to model failure times. The standard
%  Birnbaum–Saunders distribution is a two parameter continuous probability
%  distributions with support on [0,Inf], with alpha > 0 (the shape
%  parameter) and beta > 0 (the scale parameter).
%
%  There is a mixture relationship between the Birnbaum–Saunders
%  distribution and the Inverse-Gaussian distributions. Using this mixture
%  representation, the moment generating function MGF(t) and the
%  characteristic function cf(t) of the Birnbaum–Saunders distribution
%  BS(alpha,beta) can be obtained easily. For more details the
%  Birnbaum–Saunders distribution see Balakrishnan and Kundu (2017).
%
%  The algorithm cf_BirnbaumSaunders evaluates the characteristic function
%  cf(t) of Y = sum_{i=1}^N coef_i * X_i, where X_i ~ BS(alpha_i,beta_i)
%  are inedependent Birnbaum-Saunders distributed RVs, with alpha_i > 0 and
%  beta_i > 0, for i = 1,...,N.
%
%  The characteristic function of X ~ BS(alpha,beta) is given by
%   cf(t) = exp( 1/alpha^2 - 1/alpha * (1/alpha^2-2*i*t*beta)^(1/2) ) ...
%           *( 1 + (1 - 2*i*t*alpha^2*beta)^(-1/2) ) / 2.
%  Hence, the characteristic function of Y is
%   cf(t) = Prod ( cf_BirnbaumSaunders(coef_i*t,alpha_i,beta_i) ).
%
% SYNTAX:
%  cf = cf_BirnbaumSaunders(t,alpha,beta,coef,niid)
% 
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of shape parameters, alpha_i > 0. If empty, default
%          value is alpha = 1. 
%  beta  - vector of shape parameters, beta_i > 0. If empty, default
%          value is beta = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient n, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef_i * X_i is independently and identically distributed random
%          variable. If empty, default value is niid = 1. 
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution
%
% EXAMPLE 1:
%  % CF of Birnbaum-Saunders RVs with alpha = 1 and beta = 1
%  alpha  = 1;
%  beta   = 1;
%  t      = linspace(-20,20,201);
%  cf     = cf_BirnbaumSaunders(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('CF of Birnbaum-Saunders RVs with alpha = 1 and beta = 1')
%
% EXAMPLE 2:
%  % CDF/PDF of Birnbaum-Saunders RVs with alpha = 1 and beta = 1
%  alpha  = 1;
%  beta   = 1;
%  cf     = @(t) cf_BirnbaumSaunders(t,alpha,beta);
%  x      = linspace(0,8);
%  prob   = [0.9 0.95 0.975 0.99];
%  clear options;
%  options.N = 2^12;
%  options.xMin = 0;
%  result = cf2DistGP(cf,x,prob,options);
%  disp(result)
%
% EXAMPLE 3:
%  % CF of a linear combination of independent Birnbaum-Saunders RVs
%  alpha  = [1 2 3 4]/4;
%  beta   = [4 3 2 1];
%  t      = linspace(-2,2,201);
%  cf     = cf_BirnbaumSaunders(t,alpha,beta);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('CF of a linear combination of Birnbaum-Saunders RVs')
%
% EXAMPLE 4:
%  % CDF/PDF of a linear combination of independent Birnbaum-Saunders RVs
%  alpha  = [1 2 3 4]/4;
%  beta   = [4 3 2 1];;
%  cf   = @(t) cf_BirnbaumSaunders(t,alpha,beta);
%  x    = linspace(0,30);
%  prob = [0.9 0.95 0.975 0.99];
%  clear options;
%  options.N = 2^12;
%  options.xMin = 0;
%  result = cf2DistGP(cf,x,prob,options);
%  disp(result)
% 
% REFERENCES:
%
% [1] Balakrishnan, N. and Kundu, D., 2019. Birnbaum‐Saunders distribution:
%     A review of models, analysis, and applications. Applied Stochastic
%     Models in Business and Industry, 35(1), 4-49.    
%     https://onlinelibrary.wiley.com/doi/epdf/10.1002/asmb.2348

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 30-Aug-2021 22:30:54

%% ALGORITHM
% cf = cf_BirnbaumSaunders(t,alpha,beta,coef,niid);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5,   niid  = []; end
if nargin < 4,   coef  = []; end
if nargin < 3,   beta  = []; end
if nargin < 2,   alpha = []; end

if isempty(alpha), alpha = 1; end
if isempty(beta),  beta  = 1; end
if isempty(coef),  coef  = 1; end
if isempty(niid),  niid  = 1; end

%% Equal size of the parameters
[errorcode,coef,alpha,beta] = ...
    distchck(3,coef(:)',alpha(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function 
szt   = size(t);
t     = t(:);
o     = ones(size(t));
inva  = 1./alpha;

cf    = prod( exp( (o*inva.^2) - (o*inva) ...
        .* sqrt( (o*inva.^2) - 2*1i*t*(coef.*beta)) ) ...
        .* (1 + (1 - 2*1i*t*(coef.*alpha.^2.*beta)).^(-1/2))/2 , 2 );
cf    = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end