function cf = cf_Beta(t,alpha,beta,coef,niid)
%% cf_Beta 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent BETA random variables defined on the interval (0,1).  
%
%  That is, cf_Beta evaluates the characteristic function cf(t)  of  Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Beta(alpha_i,beta_i) are
%  independent RVs, with the shape parameters alpha_i > 0 and beta_i >0,
%  and with the mean = alpha_i / (alpha_i + beta)_i and the variance =
%  (alpha_i*beta_i) / ((alpha_i+beta_i)^2*(alpha_i+beta_i+1)), for i =
%  1,...,N.  
%
%  The characteristic function of X ~ Beta(alpha,beta) is 
%   cf(t) = cf_Beta(t,alpha,beta) = 1F1(alpha; alpha +beta; i*t),
%  where 1F1(.;.;.) is the Confluent hypergeometric function. Hence,the
%  characteristic function of Y  = coef(1)*X_1 + ... + coef(N)*X_N is  
%   cf(t) =  cf_X_1(coef(1)*t) * ... * cf_X_N(coef(N)*t), 
%  where X_i ~ Beta(alpha(i),beta(i)) with cf_X_i(t).
%
% SYNTAX
%  cf = cf_Beta(t,alpha,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  alpha - vector of the 'shape' parameters alpha > 0. If empty, default
%          value is alpha = 1.  
%  beta  - vector of the 'shape' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          Beta distributed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Beta_distribution
%
% EXAMPLE 1:
% % CF of a Beta RV
%   alpha = 1;
%   beta  = 3;
%   t = linspace(-50,50,501);
%   cf = cf_Beta(t,alpha,beta);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of a Beta RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a Beta RV
%   alpha = 1;
%   beta  = 3;
%   cf = @(t)cf_Beta(t,alpha,beta);
%   clear options
%   options.N = 2^12;
%   options.xMin = 0;
%   options.xMax = 1;
%   x = linspace(0,1,201);
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE 3:
% % CF of a linear combination of independent Beta RVs
%   alpha = 1;
%   beta  = 3;
%   coef  = 1./(((1:50) - 0.5)*pi).^2;
%   weights = coef/sum(coef);
%   t = linspace(-100,100,501);
%   cf = cf_Beta(t,alpha,beta,weights);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('CF of a weighted linear combination of independent Beta RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a weighted linear combination of independent Beta RVs
%   alpha = 1;
%   beta  = 3;
%   coef  = 1./(((1:50) - 0.5)*pi).^2;
%   weights = coef/sum(coef);
%   cf = @(t)cf_Beta(t,alpha,beta,weights);
%   clear options
%   options.xMin = 0;
%   options.xMax = 1;
%   x = linspace(0,1,201);
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,x,prob,options)

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 14-May-2017 12:08:24
% Rev.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf_Beta(t,alpha,beta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, alpha = []; end

%%
if isempty(beta) && ~isempty(alpha)
    beta = 1;
elseif isempty(beta) && ~isempty(coef)
    beta = 1;
elseif ~any(beta)
    beta = 1;
end

if isempty(alpha) && ~isempty(coef)
    alpha = 1;
elseif isempty(alpha) && ~isempty(beta)
    alpha = 1;
end

if isempty(coef) && ~isempty(beta)
    coef = 1;
elseif isempty(coef) && ~isempty(alpha)
    coef = 1;
end

if isempty(niid)
    niid = 1;
end

%% Equal size of the parameters   
if ~isempty(coef) && isscalar(alpha) && isscalar(beta) && isempty(niid)
    coef = sort(coef);
    m    = length(coef);
    [coef,idx] = unique(coef);
    alpha = alpha * diff([idx;m+1]);
end

[errorcode,coef,alpha,beta] = distchck(3,coef(:)',alpha(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szc = length(coef);
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:szc
    cf = cf .* Hypergeom1F1(alpha(i),alpha(i)+beta(i),1i*coef(i)*t);
end
cf  = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end