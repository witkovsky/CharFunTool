function cf = cf_LogRV_Exponential(t,lambda,coef,niid)
%%cf_LogRV_Exponential 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED EXPONENTIAL random variables (RVs) log(X),
%  where X ~ Exponential(lambda), and lambda > 0 represent the rate
%  parameter of the EXPONENTIAL distribution. 
%  
%  That is, cf_LogRV_Exponential evaluates the characteristic function
%  cf(t) of Y = coef_i*log(X_1) +...+ coef_N*log(X_N), where X_i ~
%  Exponential(lambda_i), with the parameters lambda_i > 0, for i =
%  1,...,N.
%
%  The characteristic function of Y = log(X), with X ~ Exponential(lambda) is
%  defined by cf_Y(t) = E(exp(1i*t*Y)) = E(exp(1i*t*log(X))) = E(X^(1i*t)). 
%  That is, the characteristic function can be derived from expression for
%  the r-th moment of X, E(X^r) by using (1i*t) instead of r. In
%  particular, the characteristic function of Y = log(X) is
%   cf_Y(t) = lambda^(-1i*t) .* gamma(1 + 1i*t).
%
%  Hence,the characteristic function of Y  = coef(1)*Y1 + ... + coef(N)*YN
%  is  cf_Y(t) =  cf_Y1(coef(1)*t) * ... * cf_YN(coef(N)*t), where cf_Yi(t)
%  is evaluated with the parameters alpha(i) and lambda(i).
%
% SYNTAX
%  cf = cf_LogRV_Exponential(t,lambda,coef,niid)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  lambda - vector of the 'rate' parameters, lambda > 0. If empty, default
%           value is lambda = 1.  
%  coef   - vector of the coefficients of the linear combination of the
%           log-transformed Exponential random variables. If coef is
%           scalar, it is assumed that all coefficients are equal. If
%           empty, default value is coef = 1.
%  niid   - scalar convolution coeficient niid, such that Z = Y +...+ Y is
%           sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%           coef(i) * log(X_i) is independently and identically distributed
%           random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Exponential_distribution
%
% EXAMPLE 1:
% % CF of a log-transformed Exponential RVs
%   lambda = 3/2;
%   t      = linspace(-10,10,501);
%   cf     = cf_LogRV_Exponential(t,lambda);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a log-transformed Exponential RVs')
%
% EXAMPLE 2:
% % PDF/CDF of a combination of independent log-transformed Exponential RVs
%   lambda = 3/2;
%   cf     = @(t) cf_LogRV_Exponential(t,lambda);
%   x      = linspace(-7,3);
%   prob   = [0.9 0.95 0.99];
%   clear options
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 3:
% % CF of a combination of independent log-transformed Exponential RVs
%   coef   = [1 2 3 4 5];
%   coef   = -coef/sum(coef);
%   lambda = 3/2;
%   t      = linspace(-10,10,501);
%   cf     = cf_LogRV_Exponential(t,lambda,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a combination of log-transformed Exponential RVs')
%
% EXAMPLE 4:
% % PDF/CDF of a combination of independent log-transformed Exponential RVs
%   coef   = [1 2 3 4 5];
%   coef   = -coef/sum(coef);
%   lambda = 3/2;
%   cf     = @(t) cf_LogRV_Exponential(t,lambda,coef);
%   x      = linspace(-4,2);
%   prob   = [0.9 0.95 0.99];
%   clear options
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Oct-2018 15:22:43

%% ALGORITHM
% cf = cf_LogRV_Exponential(t,lambda,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 4, niid = []; end
if nargin < 3, coef = []; end
if nargin < 2, lambda = []; end

% CF of the linear combination of independent LOG-TRANSFORMED EXPONENTIAL
% RVs expressed as a special case of the cf_LogRV_Gamma
cf = cf_LogRV_Gamma(t,1,lambda,coef,niid);

end