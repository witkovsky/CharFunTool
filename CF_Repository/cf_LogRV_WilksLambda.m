function cf = cf_LogRV_WilksLambda(t,p,n,q,coef,niid)
%% cf_LogRV_WilksLambda 
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED WILK's LAMBDA distributed random variables.
%  
%  That is, cf_LogRV_WilksLambda evaluates the characteristic function of a
%  random variable Y  = coef_1*W_1 +...+ coef_N*W_N, such that cf_Y(t) =
%  cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t), where cf_W_i(t) is CF of W_i =
%  log(Lambda_i), and each Lambda_i has central WILK's LAMBDA distribution,
%  Lambda_i ~ Lambda(p_i,n_i,q_i), for i = 1,...,N.  
%
%  In particular Lambda_i = det(E_i)/det(E_i + H_i), with  E_i and H _i  
%  being independent random matrices with central Wishart distributions,
%  E_i ~  Wp_i(n_i,Sigma_i) and H_i ~ Wp_i(q_i,Sigma_i) with n_i >= p_i and
%  Sigma_i > 0 (an unknown positive definite symmetric covariance matrix),
%  for all i = 1,...,N. 
%
%  Each particular Lambda_i statistic can be used to test specific null
%  hypothesis (by measuring its effect expressed by the matrix H_i). In
%  this case, the null hypothesis is rejected for small values of the
%  observed statistic Lambda_i, or large value of -log(Lambda_i).    
%
%  The central Wilks’ distribution of Lambda_i ~ Lambda(p_i,n_i,q_i), with
%  Lambda_i in (0,1), is defined by Lambda_i ~ Prod_{j=1}^p_i B_{i,j}, 
%  where B_{i,j} ~ Beta{(n_i+1-j)/2, q_i/2)}, i.e. B_{i,j} follow
%  independent Beta distributions for all i = 1,...,N and j = 1,...,p_i.   
%
% SYNTAX
%  cf = cf_LogRV_WilksLambda(t,p,n,q,coef,niid)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  p     - vector of the dimension parameters p = (p_1,...,p_N). If empty,
%          default value is p = (1,...,1).  
%  n     - vector of degrees of freedom of the Wishart matrices E_i, n =
%          (n_1,...,n_N). If empty, default value is n = (1,...,1).
%  q     - vector of degrees of freedom of the Wishart matrices H_i, q =
%          (q_1,...,q_N). If empty, default value is q = (1,...,1).
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  niid  - scalar convolution coeficient niid, such that Z = Y + ... + Y is
%          sum of niid random variables Y, where each Y = sum_{i=1}^N
%          coef(i) * log(X_i) is independently and identically
%          distributed random variable. If empty, default niid = 1.  
%
%  WIKIPEDIA: 
%   https://en.wikipedia.org/wiki/Wilks%27s_lambda_distribution
%
% EXAMPLE 1:
% % CF of log Wilks Lambda RV distribution Lambda(p,n,q)
%   p  = 5;
%   n  = 10; % d.f. of within SS&P
%   q  = 3;  % d.f. of between SS&P
%   t  = linspace(-10,10,201);
%   cf = cf_LogRV_WilksLambda(t,p,n,q);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of log Wilks Lambda RV with Lambda(p,n,q), p=5, n=10, q=3')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of minus log Wilks Lambda RVs 
%   p    = [5 5 5];
%   n    = [10 15 20];
%   q    = [3 2 1];
%   coef = -[10 15 20]/45;
%   t    = linspace(-20,20,201);
%   cf   = cf_LogRV_WilksLambda(t,p,n,q,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a weighted linear combination of -log Wilks Lambda RVs')
%
% EXAMPLE 3:
% % PDF/CDF of minus log Wilks Lambda RV (p=5, n=10, q=3) from its CF
%   p    = 5;
%   n    = 10;
%   q    = 3;
%   coef = -1;
%   cf   = @(t) cf_LogRV_WilksLambda(t,p,n,q,coef);
%   x    = linspace(0,5)';
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 4:
% % PDF/CDF of Wilks Lambda (p=5, n=10, q=3) from PDF of LOG-TRANSFORMED RV
%   p    = 5;
%   n    = 10;
%   q    = 3;
%   cf   = @(t) cf_LogRV_WilksLambda(t,p,n,q);
%   clear options
%   options.isPlot = 0;
%   options.xMax = 0;
%   result = cf2DistGP(cf,[],[],options);
%   [x,pdf] = LogPDF2PDF(result);
%   figure; plot(x,pdf)
%   title('PDF of the Wilks Lambda (p=5, n=10, q=3)')
%   [x,cdf] = LogCDF2CDF(result);
%   figure; plot(x,cdf)
%   title('CDF of the Wilks Lambda (p=5, n=10, q=3)')
%
% EXAMPLE 5:
% % Compare the exact distribution with the Bartlett's approximation
% % The Bartlett's approximation (see e.g. Wikipedia) is given by:
% % ((p-q+1)/2 - n)*log(Lambda(p,n,q)) ~ chi^2_{q*p}
%   p    = 15;
%   n    = 30;
%   q    = 3;
%   coef = (p-q+1)/2 - n;
%   cf   = @(t) cf_LogRV_WilksLambda(t,p,n,q,coef);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,[],prob,options);
%   disp(result)
%   x = result.x;
%   figure;plot(x,result.cdf,x,chi2cdf(x,q*p));grid
%   title('Exact CDF vs. the Bartlett approximation')
%   xlabel('transformed \Lambda(p,n,q)')
%   ylabel('CDF')
%   disp(prob)
%   disp(result.qf)
%   disp(chi2inv(prob,q*p))
%
% REFERENCES:
% [1] Witkovský, V., 2018. Exact distribution of selected multivariate test
%     criteria by numerical inversion of their characteristic functions. 
%     arXiv preprint arXiv:1801.02248.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 17-Oct-2018 18:14:11

%% ALGORITHM
% cf = cf_LogRV_WilksLambda(t,p,n,q,coef,niid)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, niid = []; end
if nargin < 5, coef = []; end
if nargin < 4, q = []; end
if nargin < 3, n = []; end
if nargin < 2, p = []; end

%%
if isempty(p), p = 1; end
if isempty(n), n = 1; end
if isempty(q), q = 1; end
if isempty(coef), coef = 1; end
if isempty(niid), niid = 1; end

%% Check size of the parameters
[errorcode,p,n,q,coef] = distchck(4,p(:)',n(:)',q(:)',coef(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination 
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:length(coef)
    alpha = (n(i)+1-(1:p(i)))/2;
    beta  = q(i)/2;
    cf = cf .* cf_LogRV_Beta(coef(i)*t,alpha,beta);
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