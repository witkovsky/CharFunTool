function cf = cf_TsallisQGaussian(t,mu,sigma,q,coef,niid)
%% cf_TsallisQGaussian 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent random variables with the Tsallis q-Gaussian distribution
%  with the location parameter mu (real value), the scale parameter sigma
%  (sigma > 0), and the shape parameter q (q < 3). The random variable with
%  q-Gaussian distribution is denoted by X ~ TsallisQGaussian(mu,sigma,q).
%
%  The parametrization used here is consistent with the parametrizatin used
%  in Wolfram MATHEMATICA, where sigma > 0 represents the scale parameter.
%  Note that this is different from the parametrization defined by Tsallis
%  (see, e.g., Wikipedia), which uses beta as the rate parameter.
%  
%  The parameter sigma used here is related to the parameter beta used in
%  the original Tsallis parametrization. Namely, beta = 1/(2*sigma^2) or
%  sigma = sqrt(1/(2*beta)). 
% 
%  The present parametriztion ensures that for q = 1, X ~
%  TsallisQGaussian(mu,sigma,q) <=> X ~ N(mu,sigma^2), i.e. it represents
%  RV with the normal distribution, with mean mu and variance sigma^2.
%
%  The algorithm cf_TsallisQGaussian evaluates the characteristic function
%  cf(t) of  Y = sum_{i=1}^N coef_i * X_i, where X_i ~
%  TsallisQGaussian(mu_i,sigma_i,q_i) are inedependent RVs, with parameters
%  mu_i, sigma_i and q_i, for i = 1,...,N, see Witkovsky (2023).
%
%  For q = 1, the characteristic function of the random variable with the
%  standard Tsallis Q-Gaussian distribution, i.e. X ~
%  TsallisQGaussian0,1,q), where q = 1, is equal to the characteristic
%  function of the standard normal distribution defined by
%   cf(t) = exp(-t^2/2).
%
%  For q < 1, the characteristic function of the random variable with
%  standard Tsallis q-Gaussian distribution, i.e. X ~
%  TsallisQGaussian0,1,q), where q < 1, is defined by
%   cf(t) = Hypergeometric0F1(theta + 1/2, -(c1*t)^2/4)
%         = gamma(theta+1/2) * 2^(theta-1/2) * (c1*t)^(-(theta+1/2)) 
%           * besselj(theta-1/2,c1*t)
%         = cf_BetaSymmetric(c1*t,theta)
%  where c1 = sqrt(2/(1-q)) and theta = (2-q)/(1-q). 
%  Note that the support of the standard Tsallis q-Gaussian distribution
%  with q < 1 is limitted to the interval [-c1,c1] = [- sqrt(2/(1-q)),
%  sqrt(2/(1-q))].  
%
%  For 1 < q < 3, the characteristic function of the random variable with
%  standard Tsallis q-Gaussian distribution, i.e. X ~
%  TsallisQGaussian0,1,q), where 1 < q < 3, is defined by
%   cf(t) = besselk( nu/2, sqrt(nu)*c2*abs(t) ) 
%           * ( sqrt(nu)*c2*abs(t) )^(nu/2)
%           / ( 2^(nu/2 - 1) * gamma(nu/2) )
%         = cf_Student(c2*t,nu)
%  where c2 = sqrt(2/(3-q)) and nu = (3-q)/(q-1).
%
%  In general, the characteristic function of the random variable with the
%  Tsallis q-Gaussian distribution, X ~ TsallisQGaussian(mu,sigma,q), where
%  q < 3, is defined by
%   cf_TsallisQGaussian(t) = exp(1i * t * mu) * cf(sigma * t),
%  where cf(t) is defined as above.
%
% SYNTAX:
%  cf = cf_TsallisQGaussian(t,mu,sigma,q,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  mu    - vector of the 'location' parameters mu in R. If empty, default
%          value is mu = 0.  
%  sigma - vector of the 'scale' parameters sigma > 0. If empty, default
%          value is sigma = 1.  
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
%  As seen, the CF of the standard Tsallis q-Gaussian, X ~
%  TsallisQGaussian0,1,q), is related the CFs of other distributions:
%  * If q = 1, cf_TsallisQGaussian(t) = cf_Normal(t).
%  * If q < 1, the Tsallis q-Gaussian distribution equals to the shifted
%    and scaled symmetric Beta distribution on [-1,1],
%    cf_TsallisQGaussian(t) = cf_BetaSymmetric(c1*t,theta) with c1 =
%    sqrt(2/(1-q)) and theta = (2-q)/(1-q). 
%  * If 1 < q < 3, the Tsallis q-Gaussian distribution equals to the
%    shifted and scaled Student t-distribution, cf_TsallisQGaussian(t) =
%    cf_Student(c2*t,nu) with c2 = sqrt(2/(3-q)) and nu = (3-q)/(q-1). 
%
%  The current version shows numerical problems for q ~ 1 (i.e. when the
%  distribution should be close to normal distribution) and for q -> 3 and
%  q -> -inf. This is due to overflow/underflow problems, respectively, of
%  the besselj and besselk functions used. In particular, when q -> 1 the
%  algorithm requires the computation of besselj(nu,x) and/or besselk(nu,x)
%  for nu -> inf.
% 
%  The algorithm for numerical inversion may have some numerical problems
%  for CFs with 2 < q < 3. Such cases require special treatment, e.g.
%  setting special option parameters in the algorithm cf2DistGP or using a
%  the GP accelerated algorithms, e.g., cf2Dist_GPA, cf2CDF_GPA,
%  cf2CDF_GPA. 
%
% WOLFRAM MATEMATICA
%   https://reference.wolfram.com/language/ref/TsallisQGaussianDistribution.html
% 
% WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Q-Gaussian_distribution
%
% EXAMPLE 1:
% % CF of a linear combination of q-Gaussian random variables
%   mu     = [0 1 2];
%   sigma  = [1 1 1];
%   q      = [-1 0.5 1.5];
%   coef   = [1 1 1]/3;
%   t      = linspace(-10,10,201);
%   cf     = cf_TsallisQGaussian(t,mu,sigma,q,coef);
%   figure; plot(t,real(cf),t,imag(cf));grid on
%   title('CF of a linear combination of q-Gaussian RVs')
%
% EXAMPLE 2:
% % PDF/CDF from the CF of a linear combination of q-Gaussian RVs
%   mu     = [0 1 2];
%   sigma  = [1 1 1];
%   q      = [-1 0.5 1.5];
%   coef   = [1 1 1]/3;
%   cf     = @(t) cf_TsallisQGaussian(t,mu,sigma,q,coef);
%   clear options
%   options.N = 2^10;
%   prob   = [0.01 0.05 0.1 0.5 0.9 0.950 .99];
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%
% EXAMPLE 3:
% % PDF/CDF from the CF of a linear combination of q-Gaussian RVs
%   mu     = [0 0 0 0 0];
%   sigma  = [5 4 3 2 1];
%   q      = [-5 -1 0 1 2];
%   coef   = [1 1 1 1 1]/5;
%   cf     = @(t) cf_TsallisQGaussian(t,mu,sigma,q,coef);
%   clear options
%   options.N = 2^14;
%   x      = linspace(-10,10,301);
%   prob   = [0.01 0.05 0.1 0.5 0.9 0.950 .99];
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 4:
% % PDF/CDF from the CF of a linear combination of q-Gaussian RVs using the
% % Tsallis parametrization
%   mu     = [0 0 0 0 0];
%   beta   = [5 4 3 2 1];
%   sigma  = sqrt(1./(2*beta));
%   q      = [-5 -1 0 1 2];
%   coef   = [1 1 1 1 1]/5;
%   cf     = @(t) cf_TsallisQGaussian(t,mu,sigma,q,coef);
%   clear options
%   options.N = 2^14;
%   x      = linspace(-10,10,301);
%   prob   = [0.01 0.05 0.1 0.5 0.9 0.950 .99];
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 5:
% % CDF from the CF of a linear combination of q-Gaussian RVs with
% % extreme values of q (here max q = 2.9) computed with using cf2CDF_GPA 
%   mu     = [0 0 0];
%   sigma  = [1 0.5 0.1];
%   q      = [ 0 1 2.9];
%   coef   = [1 1 1]/3;
%   cf     = @(t) cf_TsallisQGaussian(t,mu,sigma,q,coef);
%   clear options
%   options.isAccelerated = true; 
%   % Plot the integrand functions used for accelerated computation of CDF
%   % options.isPlot = true;
%   x = [1e+10 1e+20 1e+30 1e+40 1e+50 1e+60 1e+70 1e+80 1e+90]';
%   cdf = cf2CDF_GPA(cf,x,options);
%   Table = table(x,cdf);
%   disp(Table)
%
% REFERENCES:
% [1] WITKOVSKY. V. (2023). Characteristic Function of the Tsallis
%     q-Gaussian and Its Applications in Measurement Science and Metrology.
%     In:  MEASUREMENT 2023, Proceedings of the 14th International
%     Conference on Measurement. Smolenice, Slovakia, May 28-31.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-Mar-2023 10:11:54

%% ALGORITHM
% cf = cf_TsallisQGaussian(t,mu,sigma,q,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, q = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

if isempty(mu), mu = 0; end
if isempty(sigma), sigma = 1; end
if isempty(q), q = 1; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Equal size of the parameters   
[errorcode,coef,mu,sigma,q] = distchck(4,coef(:)',mu(:)',sigma(:)',q(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end
szt      = size(t);
t        = t(:);

idE1 = q==1;
idL1 = q<1;
idG1 = q>1 & q<3;
idG3 = q >= 3;

cf = 1;

if any(idE1)
    MU    = mu(idE1);
    SIGMA = sigma(idE1);
    COEF  = coef(idE1);
    cf    = cf .* cf_Normal(t,MU,SIGMA,COEF);
end

if any(idG1)
    Q     = q(idG1);
    C1    = sqrt(2./(3-Q));
    NU    = (3-Q)./(Q-1); 
    MU    = mu(idG1);
    SIGMA = C1.*sigma(idG1);
    COEF  = coef(idG1);
    cf    = cf .* cf_Student(t,NU,MU,SIGMA,COEF);
end

if any(idL1)
    Q     = q(idL1);
    C2    = sqrt(2./(1-Q));
    THETA = (2-Q)./(1-Q);
    MU    = mu(idL1);
    SIGMA = sigma(idL1);
    COEF  = coef(idL1);
    cf    = cf .* cf_BetaSymmetric(t,THETA,SIGMA.*C2.*COEF) ...
               .* exp(1i*t*sum(MU.*COEF));
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