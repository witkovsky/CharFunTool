function cf = cf_TsallisQGaussian(t,mu,beta,q,coef,niid)
%% cf_TsallisQGaussian 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent random variables with the Tsallis q-Gaussian distribution
%  with real mean mu, scale parameter beta > 0, and deformation (shape)
%  parameter q < 3. 
% 
%  This parametrization is consistent with Wolfram MATHEMATICA, where beta
%  > 0 represents the scale parameter. Note that this is different from the
%  original parametrization defined by Tsallis and explained, e.g., in
%  Wikipedia, which uses beta as the rate parameter. 
%  
%  The parameter beta used here is related to the original parametrization,
%  say BETA. Namely, BETA = 1/(2*beta^2) or beta = sqrt(1/(2*BETA)). This
%  parametriztion ensures that for q = 1, cf_TsallisQGaussian(mu,beta,q) ~
%  N(mu,beta), i.e. with the normal distribution with mean mu and variance
%  sigma^2 = beta^2.
%
%  cf_TsallisQGaussian evaluates the characteristic function cf(t) of  Y =
%  sum_{i=1}^N coef_i * X_i, where X_i ~ TsallisQGaussian(mu_i,beta_i,q_i)
%  are inedependent RVs, with parameters mu_i, beta_i and q_i, for i =
%  1,...,N.
%
%  For q = 1, the characteristic function of the random variable with the
%  standard Tsallis Q-Gaussian distribution, i.e. X ~
%  TsallisQGaussian(0,1,q), where q = 1, is equal to the characteristic
%  function of the standard normal distribution defined by
%   cf(t) = exp(-t^2/2).
%
%  For q < 1, the characteristic function of the random variable with
%  standard Tsallis Q-Gaussian distribution, i.e. X ~
%  TsallisQGaussian(0,1,q), where q < 1, is defined by 
%  cf(t) = Hypergeometric0F1(3/2 + 1/(1-q), -t^2/(2 * (1-q)))
%        = 2^(nu/4) * (t/sqrt(1-q))^(-nu/2) * ...
%          besselj(nu/2, sqrt(2/(1-q))*t) * gamma(3/2 + 1/(1-q)),
%  where nu = (3-q)/(1-q). Note that the support of the distribution is
%  [- beta*sqrt(2/(1-q)), beta*sqrt(2/(1-q))].
%
%  For 1 < q < 3, the characteristic function of the random variable with
%  standard Tsallis Q-Gaussian distribution, i.e. X ~
%  TsallisQGaussian(0,1,q), where 1 < q < 3, is defined by 
%   cf(t) = 2^( (7-5*q)/(4-4*q) ) * (q-1)^(-df/4) * |t|^(df/2) * ...
%           besselk(df/2, sqrt(2/(q-1))*|t|) / gamma(nu/2)
%  where df = (3-q)/(q-1). 
%  
%  This CF is related to the CF of the Student t-distribution. In
%  particular, the characteristic function of the random variable with
%  Tsallis Q-Gaussian distribution, i.e. X ~ TsallisQGaussian(mu,beta,q),
%  where 1 < q < 3, can be evaluated as
%   cf(t) = cf_Student(t,df,mu,sigma),
%  where df = (3-q)/(q-1); sigma = beta*sqrt(2/(3-q)).
%
% SYNTAX:
%  cf = cf_TsallisQGaussian(t,mu,beta,q,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the 'location' parameters mu in R. If empty, default
%          value is mu = 0.  
%  beta  - vector of the 'scale' parameters beta > 0. If empty, default
%          value is beta = 1.  
%  q     - vector of the 'shape' parameters q < 3. If empty, default
%          value is q = 1.
%  coef  - vector of the coefficients of the linear combination of the
%          Normal random variables. If coef is scalar, it is assumed
%          that all coefficients are equal. If empty, default value is
%          coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid iid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * X_i is independently and identically distributed
%          random variable. If empty, default value is niid = 1.     
%
% REMARKS:
%  The current version shows numerical problems for q ~ 1 (i.e. when the
%  distribution should be close to normal distribution) and for q -> 3 and
%  q -> -inf. This is due to overflow/underflow problems, respectively, of
%  the besselj and besselk functions used. In particular, when q -> 1 the
%  algorithm requires the computation of besselj(nu,x) and/or besselk(nu,x)
%  for nu -> inf.
% 
%  Note that X ~ TsallisQGaussian(mu,beta,q) with 1 < q < 2 is a heavy tail
%  distribution. In particular, if q = 2, we get the Cauchy distribution.
%  If 2 < q < 3, X ~ TsallisQGaussian(mu,beta,q) is in a region of
%  extremely heavy tails, corresponding to a Student t- distribution with
%  degrees of freedom df such that df -> 0 when q -> 3 and df -> 1 when q
%  -> 2. The algorithm may have some numerical problems for 2 < q < 3 that
%  require special treatment (option parameters). 
%
% WOLFRAM MATEMATICA
%   https://reference.wolfram.com/language/ref/TsallisQGaussianDistribution.html
% 
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Q-Gaussian_distribution
%
% EXAMPLE 1:
% % CF of a linear combination of Q-Gaussian random variables
%   mu   = [0 1 2];
%   beta = [1 1 0.5];
%   q    = [-1 0.5 1.25];
%   coef = [1 1 1]/3;
%   t = linspace(-10,10,201);
%   cf = cf_TsallisQGaussian(t,mu,beta,q,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of Q-Gaussian RVs')
%
% EXAMPLE 2:
% % PDF/CDF from the CF of a linear combination of Q-Gaussian RVs
%   mu   = [0 1 2];
%   beta = [1 1 0.5];
%   q    = [-1 0.5 1.5];
%   coef = [1 1 1]/3;
%   cf   = @(t) cf_TsallisQGaussian(t,mu,beta,q,coef);
%   clear options
%   options.N = 2^10;
%   prob = [0.01 0.05 0.1 0.5 0.9 0.950 .99];
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%
% EXAMPLE 3:
% % PDF/CDF from the CF of a linear combination of Q-Gaussian RVs
%   mu   = [0 0 0 0 0];
%   beta = [5 4 3 2 1];
%   q    = [-5 -1 0 1 2];
%   coef = [1 1 1 1 1]/5;
%   cf   = @(t) cf_TsallisQGaussian(t,mu,beta,q,coef);
%   clear options
%   options.N = 2^14;
%   x = linspace(-10,10,301);
%   prob = [0.01 0.05 0.1 0.5 0.9 0.950 .99];
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 4:
% % PDF/CDF from the CF of a linear combination of Q-Gaussian RVs using the
% % Tsallis parametrization
%   mu   = [0 0 0 0 0];
%   BETA = [5 4 3 2 1];
%   beta = sqrt(1./(2*BETA));
%   q    = [-5 -1 0 1 2];
%   coef = [1 1 1 1 1]/5;
%   cf   = @(t) cf_TsallisQGaussian(t,mu,beta,q,coef);
%   clear options
%   options.N = 2^14;
%   x = linspace(-10,10,301);
%   prob = [0.01 0.05 0.1 0.5 0.9 0.950 .99];
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 14-Feb-2023 15:12:06

%% ALGORITHM
% cf = cf_TsallisQGaussian(t,mu,beta,q,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, q = []; end
if nargin < 3, beta = []; end
if nargin < 2, mu = []; end

if isempty(mu), mu = 0; end
if isempty(beta), beta = 1; end
if isempty(q), q = 1; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Equal size of the parameters   
[errorcode,coef,mu,beta,q] = distchck(4,coef(:)',mu(:)',beta(:)',q(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end
szt      = size(t);
t        = t(:);
o        = ones(size(t));

idE1 = q==1;
idL1 = q<1;
idG1 = q>1 & q<3;
idG3 = q >= 3;

cf = 1;

if any(idE1)
    MU   = mu(idE1);
    SIG  = beta(idE1);
    COEF = coef(idE1);
    cf   = cf .* cf_Normal(t,MU,SIG,COEF);
end

if any(idG1)
    DF   = (3-q(idG1))./(q(idG1)-1); 
    MU   = mu(idG1);
    SIG  = beta(idG1).*sqrt(2./(3-q(idG1)));
    COEF = coef(idG1);
    cf   = cf .* cf_Student(t,DF,MU,SIG,COEF);
end

if any(idL1)
    Q    = q(idL1);
    NU   = (3-Q)./(1-Q);
    MU   = mu(idL1);
    BETA = beta(idL1);
    COEF = coef(idL1);
% %  cf(t) = 2^(nu/4) * (t/sqrt(1-q))^(-nu/2) * ...
% %          besselj(nu/2, sqrt(2/(1-q))*t) * gamma(3/2 + 1/(1-q)),
    aux = 2.^(o*NU/4);
    aux = aux .* (t*(BETA .* COEF ./sqrt(1-Q))).^(-o*NU/2);
    aux = aux .* besselj((o*NU)/2, t*(BETA .* COEF .* sqrt(2./(1-Q))));
    aux = aux .* gamma(3/2 + o *(1./(1-Q)));
    aux = prod(aux .* exp((1i*t * (BETA .* COEF .* MU))),2);
    cf  = cf .* aux;
end

if any(idG3)
     warning('VW:CharFunTool:cf_TsallisQGaussian', ...
            'Random variables with q > 3 have been ingnored');
end

cf = reshape(cf,szt);
cf(t==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end