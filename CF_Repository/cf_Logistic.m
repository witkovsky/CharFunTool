function cf = cf_Logistic(t,mu,beta,coef,niid)
%% cf_Logistic  
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOGISTIC random variables with location parameter mu (real)
%  and scale parameter beta > 0. 
%
%  That is, cf_Logistic evaluates the characteristic function cf(t) of Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ Logistic (mu_i,beta_i) are
%  inedependent RVs, with real location parameters mu_i and the scale
%  parameters beta_i > 0, for i = 1,...,N. 
%
% PDF of X ~ Logistic(mu,beta) is given by 
%   pdf(x) = (1/beta)*exp(-(x-mu)/beta) / (1 + exp(-(x-mu)/beta))^2
% and the characteristic function of X is defined by
%   cf(t) = (pi*t*beta) * exp(1i*t*mu) * csch(pi*t*beta)
% where csch(.) is the hyperbolic cosecant function. Hence, the
% characteristic function of Y is specified by
%   cf(t) = Prod( cf(coef_i*t|mu_i,beta_i) ).
% For more details see WIKIPEDIA, and for possible applications also [1].
%
% SYNTAX:
%  cf = cf_Logistic(t,mu,beta,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of real location parameters. If empty, default value is 
%          mu = 0.   
%  beta  - vector of the scale parameters beta > 0. If empty, default
%          value is beta = 1.  
%  coef  - vector of the coefficients of the linear combination of the
%          LOGISTIC random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.   
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Logistic_distribution.
%
% EXAMPLE 1:
% % CF of the Logistic RV
%   mu   = 0;
%   beta = 1;
%   t    = linspace(-10,10,201);
%   cf   = cf_Logistic(t,mu,beta);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('Characteristic function of the Logistic RVs')
%
% EXAMPLE 2:
% % PDF/CDF of the Logistic RV
%   mu   = 0;
%   beta = 1;
%   x    = linspace(-10,10,101);
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Logistic(t,mu,beta);
%   result = cf2DistGP(cf,x,prob);
%
% EXAMPLE 3:
% % PDF/CDF of the linear combination of Logistic RVs 
%   mu   = [-4 -1 2 3];
%   beta = [0.1 0.2 0.3 0.4];
%   coef = [1 2 3 4];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Logistic(t,mu,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%
% EXAMPLE 4:
% % PDF/CDF of the linear combination of Logistic RVs 
%   mu   = [-10 10 20 30 40];
%   beta = [1 2 3 4 5];
%   coef = [1/2 1 3/4 5 1];
%   prob = [0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 0.99, 0.995, 0.999];
%   cf   = @(t) cf_Logistic(t,mu,beta,coef);
%   clear options
%   options.N = 2^12;
%   result = cf2DistGP(cf,[],prob,options);
%
% REFERENCES
% [1] Popovic, B.V., Mijanovic, A. and Genç, A.I., 2020. On linear
%     combination of generalized logistic random variables with an
%     application to financial returns. Applied Mathematics and
%     Computation, 381, p.125314.   

% (c) 2020 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Apr-2020 13:47:42

%% ALGORITHM
% cf = cf_Logistic(t,mu,beta,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 5);
if nargin < 5, niid = []; end
if nargin < 4, coef = []; end
if nargin < 3, beta = []; end
if nargin < 2, mu = []; end

%%
if isempty(mu), mu = 0; end
if isempty(beta), beta = 1; end
if isempty(coef), coef = 1; end


%% Equal size of the parameters   

[errorcode,coef,mu,beta] = distchck(3,coef(:)',mu(:)',beta(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function
szt = size(t);
t   = t(:);
cf  = prod( pi*t*(coef.*beta) .* exp(1i*t*(coef.*mu)) .* ...
      csch(pi*t*(coef.*beta)),2);    
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