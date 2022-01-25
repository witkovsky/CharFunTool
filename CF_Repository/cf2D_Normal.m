function cf = cf2D_Normal(t,mu,sigma,rho,coef,niid)
%% cf2D_Normal  
%  Characteristic function of a linear combination (resp. convolution) of
%  independent BIVARIATE NORMAL (GAUSSIAN) random variables, say X =
%  (X1,X2), with location parameters mu = (mu1,mu2) (real), and scale
%  parameters specified by the covariance matrix Sigma = [sigma1^2,cov;
%  cov,sigma2^2], where (sigma1,sigma2) > 0 are standard deviations of X1
%  and X2 and the covariance cov is specified by the correlation
%  coefficient, rho = cov(X1,X2) / sigma1*sigma2. 
% 
%  Notice, that the bivariate normal distribution is a special case of the
%  multivariate normal distribution. Here, we consider the multivariate
%  normal distribution and its characteristic function as specified at
%  WIKIPEDIA https://en.wikipedia.org/wiki/Multivariate_normal_distribution
%
%  That is, cf2D_Normal evaluates the characteristic function cf(t1,t2) of
%  the random vector  Y = (Y1,Y2) =  sum_{i=1}^N (coef1_i * X1_i,coef2_i
%  *X2_i), where X_i = (X1_i,X2_i) ~ Normal((mu1_i,nu2_i),
%  (sigma1_i,sigma2_i)) are inedependent bivariate random vectors, with
%  location parameters mu_i = (mu1_i,nu2_i) and the scale parameters
%  sigma_i = (sigma1_i,sigma2_i) > 0, and the correlation coefficients
%  rho_i for i = 1,...,N.
% 
%  In fact, the characteristic function of Y, cf_Y(t1,t2), is specified as
%  cf_Y(t1,t2) = cf_X1(coef11*t1,coef21*t2)*...*cf_XN(coef1N*t1,coef2N*t2).
%  
% SYNTAX:
%  cf = cf2D_Normal(t,mu,sigma,rho,coef,niid)
%
% INPUTS:
%  t      - (n x 2)-matrix t = [t1 t2] or a cell t = {t1 t2} of real
%           values, where the CF is evaluated. If t is cell of two vectors,
%           then cf is a (n1 x n2)-matrix where n1 = length(t1), n2 =
%           length(t2). 
%  mu     - matrix mu = [mu1 mu2] of real location parameters. If empty,
%           default value is mu = [0 0].   
%  sigma  - sigma = [sigma1 sigma2], standard deviations of the marginal
%           random variables X1 and X2, sigma > 0. If empty, default value
%           is sigma = [1 1].   
%  rho    -  correlation coefficient, rho = cov(X1,X2)/ (sigma1*sigma2).
%            -1 < rho < 1. If empty, default value is rho = 0.   
%  coef   - matrix of the coefficients, coef = [coef1 coef2], of the linear
%           combination of the BIVARIATE LOGISTIC random vectors. If coef
%           is (1x2)-vector, it is assumed that all coefficients are equal.
%           If empty, default value is coef = [1 1].
%  niid   - scalar convolution coeficient niid, such that Z = Y + ... + Y
%           is sum of niid iid random variables Y, where each Y =
%           sum_{i=1}^N coef(i) * X_i is independently and identically
%           distributed random variable. If empty, default value is niid =
%           1.
%
% WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
%
% EXAMPLE 1:
% % CF of the Bivariate Normal RV
%   mu    = [0 0];
%   sigma = [1 1];
%   rho   = 0.5;
%   cf    = @(t) cf2D_Normal(t,mu,sigma,rho);
%   tt    = linspace(-4,4,51);
%   [t1 t2] = meshgrid(tt,tt);
%   t = [t1(:) t2(:)];
%   figure; 
%   contour(t1,t2,real(reshape(cf(t),51,51)),'ShowText','on')
%   title('Contour Plot of the CF of the Bivariate Normal Distribution')
%
% EXAMPLE 2:
% % PDF/CDF of Bivariate Normal Distribution by numerical inversion of CF
%   mu    = [0 0];
%   sigma = [1 1];
%   rho   = 0;
%   cf    = @(t) cf2D_Normal(t,mu,sigma,rho);
%   clear options;
%   options.isInterp = true;
%   result = cf2Dist2D(cf,[],options)
%
% EXAMPLE 3:
% % PDF/CDF of a convolution of Bivariate Normal Distributions
%   mu    = [0 0; 0 1];
%   sigma = [1 1; 2 2];
%   rho   = [0.5; 0.9];
%   cf    = @(t) cf2D_Normal(t,mu,sigma,rho);
%   clear options;
%   options.isInterp = true;
%   options.chebyPts = 51;
%   result = cf2Dist2D(cf,[],options) 
%
% EXAMPLE 4:
% % PDF/CDF of a mixture of Bivariate Normal Distributions
%   mu    = [-1 0; 0 1; 1 0 ];
%   sigma = [0.5 1; 2 3; 0.5 1];
%   rho   = [0.5; 0.9; 0];
%   w     = [0.3 0.5 0.2];
%   cf1   = @(t) cf2D_Normal(t,mu(1,:),sigma(1,:),rho(1));
%   cf2   = @(t) cf2D_Normal(t,mu(2,:),sigma(2,:),rho(2));
%   cf3   = @(t) cf2D_Normal(t,mu(3,:),sigma(3,:),rho(3));
%   cf    = @(t) w(1)*cf1(t) +  w(2)*cf2(t) +  w(3)*cf3(t);
%   clear options;
%   options.isInterp = true;
%   options.chebyPts = 101;
%   result = cf2Dist2D(cf,[],options) 

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 25-Jan-2022 12:47:42

%% ALGORITHM
% cf = cf2D_Normal(t,mu,sigma,rho,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, rho = []; end
if nargin < 3, sigma = []; end
if nargin < 2, mu = []; end

%%
if isempty(mu), mu = [0,0]; end
if isempty(sigma), sigma = [1,1]; end
if isempty(rho), rho = 0; end
if isempty(coef), coef = [1,1]; end


%% Equal size of the parameters   

[err1,m1,b1,rho,c1] = distchck(4,mu(:,1)',sigma(:,1)',rho(:)',coef(:,1)');

[err2,m2,b2,c2] = distchck(3,mu(:,2)',sigma(:,2)',coef(:,2)');
if err1 + err2 > 0
    error(message('InputSizeMismatch'));
end
mu   = [m1; m2];
sigma = [b1; b2];

coef = [c1; c2];
N = length(c1);

%% Characteristic function
if iscell(t)
    sz = [length(t{1}),length(t{2})];    
    [t2,t1] = meshgrid(t{2},t{1});
    t  = [t1(:),t2(:)];
else
    szt = size(t);
    sztMin = min(szt(1),szt(2));
    sztMax = max(szt(1),szt(2));    
    switch sztMin
        case 1 % if t is N-vector then it is assumed that t1 = t and t2 = t
               % and we create new t = [t1 t2]
            if sztMax > 2
                sz = size(t);
                t = t(:);
                t = [t t];
            elseif sztMax == 1
                sz = [1,1];
                t = [t t];
            else
                sz = [1,1];
                t = t(:)';
            end
        case 2 % t is 2xN or Nx2 matrix
            if sztMax > 2
                if szt(1) == 2 % if t is 2xN matrix transpose it to the Nx2 matrix
                    t = t';
                    sz = [1,szt(2)];
                else
                    sz = [szt(1),1];
                end
            end
        otherwise
            error(message('InputSizeMismatch'));
    end
end

cf = 0;
for i = 1:N
    cf = cf + cfLog(t(:,1),t(:,2),mu(1,i),mu(2,i),...
        sigma(1,i),sigma(2,i),rho(i),coef(1,i),coef(2,i));
end
cf = exp(cf);
cf = reshape(cf,sz);

%cf(t(:,1)==0,t(:,2)==0) = 1;

if ~isempty(niid)
    if isscalar(niid)
        cf = cf .^ niid;
    else
        error('niid should be a scalar (positive integer) value');
    end
end

end
%% Function cfLog
function cf = cfLog(t1,t2,mu1,mu2,sigma1,sigma2,rho,coef1,coef2)

z = (sum([t1*coef1*sigma1^2 + t2*coef2*rho*sigma1*sigma2, ...
    t1*coef1*rho*sigma1*sigma2 + t2*coef2*sigma2^2] ...
      .* [t1*coef1, t2*coef2],2))/2;
cf = 1i*mu1*coef1*t1 + 1i*mu2*coef2*t2 - z;
end