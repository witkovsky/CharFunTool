function cf = cf_GeneralizedLindley(t,alpha,beta,gamma,coef,niid)
%% cf_GeneralizedLindley 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent GENERALIZED LINDLEY random variables (RVs) where X ~
%  GL(alpha,beta,gamma). Here gamma > 0 is a shape parameter, alpha > 0 is
%  scale parameter and beta >= 0 is mixing parameter.  
%  
%  The original Lindley distribution was proposed in Lindley (1958). The
%  Lindley distribution is expressed as a mixture of two random variables
%  according to the exponential distribution Exp(alpha) and the gamma
%  distribution Gamma(2,alpha) with probabilities alpha/(1 + alpha) and
%  1/(1 + alpha). 
%
%  The generalized (hree-parameter) Lindley distribution was proposed by
%  Zakerzadeh and Dolati (2009), for present notation see . The generalized
%  Lindley distribution is the mixture of Gamma(gamma,alpha) and
%  Gamma(gamma+1,alpha) with probabilities alpha/(alpha+beta) and
%  beta/(alpha+beta), where where alpha > 0, gamma > 0, and beta >= 0. The
%  probability density function (PDF) of GL(alpha,beta,gamma) distribution
%  is given (for x > 0) by
%  pdf(x;alpha,beta,gamma) = (alpha^2*(alpha*x)^(gamma-1)*(gamma +
%                  beta*x)*exp(-alpha*x)) / ((alpha+beta)*Gamma(gamma+1)).
%  The characteristic function of X ~ GL(alpha,beta,gamma) is given by
%  cf_X(t) = (alpha/(alpha-1i*t))^(gamma+1)*(alpha+beta-1i*t)/(alpha+beta).
%
%  That is, cf_GeneralizedLindley evaluates the characteristic function
%  cf(t) of  Y =  coef_1*X_1 +...+ coef_N*X_N, where X_i ~
%  GL(alpha_i,beta_i,gamma_i) are inedependent RVs, with the parameters
%  alpha_i > 0, gamma_i > 0, and beta_i >= 0, for i = 1,...,N. Hence,the
%  characteristic function of Y  = coef(1)*X1 + ... + coef(N)*XN  is
%  cf_Y(t) =  cf_X1(coef(1)*t) * ... * cf_XN(coef(N)*t), where cf_Xi(t) 
%  is evaluated with the parameters alpha(i), beta(i) and gamma(i).
%
% SYNTAX
%  cf = cf_GeneralizedLindley(t,alpha,beta,gamma,coef,niid)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  alpha  - vector of the 'scale' parameters, alpha > 0. If empty, default
%           value is alpha = 1.  
%  beta   - vector of the 'mixing' parameters, beta >= 0. If empty, default
%           value is beta = 1. 
%  gamma  - vector of the 'shape' parameters, gamma > 0. If empty,
%           default value is gamma = 1.  
%  coef   - vector of the coefficients of the linear combination of the
%           random variables. If coef is scalar, it is assumed that all
%           coefficients are equal. If empty, default value is coef = 1.
%  niid   - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%           sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%           coef(i) * X_i is independently and identically distributed
%           random variable. If empty, default value is niid = 1.  
%
% EXAMPLE (CF of a linear combination of generalized Lindley RVs)
%   alpha  = [2.680, 1.640];
%   beta   = [0.208, 2.852];
%   gamma  = [1.903, 0.552];
%   t      = linspace(-10,10,201);
%   cf     = cf_GeneralizedLindley(t,alpha,beta,gamma);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('Characteristic function of a linear combination of GL RVs')
%
% EXAMPLE (See Kitani et al (2021). Case 1 / n = 2)
% % PDF/CDF of a linear combination of GL RVs from the CF by cf2DistGP
%   alpha  = [2.680, 1.640];
%   beta   = [0.208, 2.852];
%   gamma  = [1.903, 0.552];
%   cf     = @(t) cf_GeneralizedLindley(t,alpha,beta,gamma);
%   prob   = [0.6 0.7 0.8 0.9 0.95 0.975 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp([prob;result.qf])
%
% EXAMPLE (See Kitani et al (2021). Case 1 / n = 5)
% % PDF/CDF of a linear combination of GL RVs from the CF by cf2DistGP
%   alpha  = [2.984, 0.510, 2.241, 1.286, 1.641];
%   beta   = [1.705, 2.727, 1.569, 1.419, 2.292];
%   gamma  = [1.548, 0.735, 2.754, 0.975, 2.651];
%   cf     = @(t) cf_GeneralizedLindley(t,alpha,beta,gamma);
%   prob   = [0.6 0.7 0.8 0.9 0.95 0.975 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp([prob;result.qf])
%
% EXAMPLE (See Kitani et al (2021). Case 1 / n = 10)
% % PDF/CDF of a linear combination of GL RVs from the CF by cf2DistGP
%   alpha  = [0.171,1.765,0.595,0.504,0.481,1.564,1.527,1.821,1.678,1.276];
%   beta   = [2.373,2.037,0.190,1.052,2.938,0.656,1.647,1.221,2.511,2.752];
%   gamma  = [2.198,2.728,2.991,1.985,0.503,2.778,1.123,0.265,1.851,2.154];
%   cf     = @(t) cf_GeneralizedLindley(t,alpha,beta,gamma);
%   prob   = [0.6 0.7 0.8 0.9 0.95 0.975 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp([prob;result.qf])
%
% EXAMPLE (See Kitani et al (2021). Case 2 / n = 2)
% % PDF/CDF of a linear combination of GL RVs from the CF by cf2DistGP
%   alpha  = [1.322, 1.042];
%   beta   = [0.406, 0.310];
%   gamma  = [1.185, 0.994];
%   cf     = @(t) cf_GeneralizedLindley(t,alpha,beta,gamma);
%   prob   = [0.6 0.7 0.8 0.9 0.95 0.975 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp([prob;result.qf])
%
% EXAMPLE (See Kitani et al (2021). Case 2 / n = 5)
% % PDF/CDF of a linear combination of GL RVs from the CF by cf2DistGP
%   alpha  = [0.655, 1.608, 0.322, 0.798, 1.507];
%   beta   = [1.359, 0.273, 1.873, 0.870, 0.291];
%   gamma  = [0.448, 3.524, 0.605, 1.564, 0.925];
%   cf     = @(t) cf_GeneralizedLindley(t,alpha,beta,gamma);
%   prob   = [0.6 0.7 0.8 0.9 0.95 0.975 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp([prob;result.qf])
%
% EXAMPLE (See Kitani et al (2021). Case 2 / n = 10)
% % PDF/CDF of a linear combination of GL RVs from the CF by cf2DistGP
%   alpha  = [1.398,4.367,1.030,0.916,1.347,2.019,0.596,0.279,1.602,0.569];
%   beta   = [0.331,0.418,0.885,1.088,1.767,1.675,0.722,0.550,1.938,0.668];
%   gamma  = [2.758,0.215,1.252,0.505,2.212,4.133,0.370,1.022,2.270,1.243];
%   cf     = @(t) cf_GeneralizedLindley(t,alpha,beta,gamma);
%   prob   = [0.6 0.7 0.8 0.9 0.95 0.975 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp([prob;result.qf])
%
% REFERENCES:
% [1] Lindley, D.V., 1958. Fiducial distributions and Bayes' theorem.
%     Journal of the Royal Statistical Society: Series B (Methodological),
%     20(1), 102-107.  
% [2] Zakerzadeh, H. and Dolati, A., 2009. Generalized Lindley
%     distribution. Journal of Mathematical Extension 3, 13–25. 
% [3] Kitani, M., Murakami, H. and Hashiguchi, H., 2021. The distribution
%     of the sum of independent and non-identically generalized Lindley
%     random variables. Submitted to Communications in Statistics - Theory
%     and Methods. Manuscript ID: LSTA-2021-0355.
%
% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 13-Jun-2021 21:21:43

%% ALGORITHM
% cf = cf_GeneralizedLindley(t,alpha,beta,gamma,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid   = []; end
if nargin < 5, coef   = []; end
if nargin < 4, gamma  = []; end
if nargin < 3, beta   = []; end
if nargin < 2, alpha  = []; end

%%
if isempty(alpha)
    alpha = 1;
end

if isempty(beta)
    beta = 1;
end

if isempty(gamma)
    gamma = 1;
end

if isempty(coef)
    coef = 1;
end

%% Check size of the parameters
[errorcode,coef,alpha,beta,gamma] = ...
    distchck(4,coef(:)',alpha(:)',beta(:)',gamma(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);
cf  = 1;
for i = 1:length(coef)
    cf = cf .* (alpha(i)./(alpha(i) - 1i*t)).^(gamma(i)+1) .* ...
        (alpha(i) + beta(i) - 1i*t) ./ (alpha(i) + beta(i));
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