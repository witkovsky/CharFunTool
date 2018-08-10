function cf = cf_LogRV_WilksLambdaNC(t,p,n,q,delta,coef,MAX)
%% cf_LogRV_WilksLambdaNC
%  Characteristic function of a linear combination (resp. convolution) of
%  independent LOG-TRANSFORMED WILK's LAMBDA distributed random variables,
%  with non-central distributions specified by the parameters p, n, q, the
%  non-centrality parameters delta, and the coefficients coef.
%
%  That is, cf_LogRV_WilksLambdaNC evaluates the characteristic function of
%  a random variable Y  = coef_1*W_1 +...+ coef_N*W_N, such that cf_Y(t) =
%  cf_W_1(coef_1*t) *...* cf_W_N(coef_N*t), where cf_W_i(t) is CF of W_i =
%  log(Lambda_i), and each Lambda_i has non-central WILK's LAMBDA
%  distribution, Lambda_i ~ Lambda(p_i,m_i,n_i,delta_i), for i = 1,...,N.
%
%  In particular, Lambda_i = det(E_i)/det(E_i + H_i), with  E_i and H _i
%  being independent random matrices with Wishart distributions, E_i ~
%  Wp_i(m_i,Sigma_i) with central Wishart distribution, and H_i ~
%  Wp_i(n_i,Sigma_i,delta) with non-central Wishart distribution, with n_i
%  >= p_i and  Sigma_i > 0 (an unknown positive definite symmetric
%  covariance matrix), for all i = 1,...,N. 
%
%  The non-central distribution of MINUS LOG-TRANSFORMED WILK's LAMBDA
%  STATISTIC, say L ~ Lambda(p,n,q,delta), with L in (0,1), is specified by
%  its characteristic function   
%  cf_{-log(L)}(t) = cf_LogRV_Beta(-t,n/2,q/2)
%          * ... * cf_LogRV_Beta(-t,(n+1-p)/2,q/2)
%          * Hypergeom1F1Mat(-i*t,-i*t+(n+q)/2,-delta/2),
%  i.e. the characteristic function of the non-central distribution differs
%  from the central distribution by a factor specified by the confluent
%  hypergeometric function 1F1(a;b;X) of a matrix argument X, with the
%  complex (vector) parameters specified as a = -i*t, b = -i*t+(n+q)/2, and
%  the matrix argument X = -delta/2, where delta is the non-centrality
%  (matrix) parameter of the distribution.  
%
% SYNTAX
%  cf = cf_LogRV_WilksLambdaNC(t,p,n,q,delta,coef,MAX)
%
% INPUTS:
%  t     - vector or array of real values, where the CF is evaluated.
%  p     - vector of the dimension parameters p = (p_1,...,p_N). If empty,
%          default value is p = (1,...,1).
%  n     - vector of degrees of freedom of the Wishart matrices E_i, n =
%          (n_1,...,n_N). If empty, default value is n = (1,...,1).
%  q     - vector of degrees of freedom of the Wishart matrices H_i, q =
%          (q_1,...,q_N). If empty, default value is n = (1,...,1).
%  delta - p.s.d. matrix or vector of nonnegative eigenvalues or cell array
%          of matrices or vectors of eigenvalues of the non-centrality
%          parameters, delta = {delta_1,...,delta_N). Default value is
%          delta = [].
%  coef  - vector of the coefficients of the linear combination of the
%          log-transformed random variables. If coef is scalar, it is
%          assumed that all coefficients are equal. If empty, default value
%          is coef = 1.
%  MAX   - the maximum number of partitions used for computing the
%          hypergeometric 1F1 function with matrix argument, for more
%          details see HypergeompFqMat.m. If MAX is empty, default value is
%          set to MAX = 0. If MAX = 0, then the algorithm uses the fast
%          (approximate) method for computing the confluent hypergeometric
%          function 1F1(a;b;X), see Hypergeom1F1MatApprox, otherwise the
%          algorithm uses the algorithm Hypergeom1F1Mat based on computing
%          the truncated expansion series of 1F1(a;b;X) with MAX number of
%          partitions. 
%
% WIKIPEDIA:
%  https://en.wikipedia.org/wiki/Wilks%27s_lambda_distribution
%
% EXAMPLE 1:
% % CF of log of noncentral Wilks Lambda RV distribution Lambda(p,n,q,delta)
%   p    = 5;
%   n    = 10; % d.f. of within SS&P
%   q    = 3;  % d.f. of between SS&P
%   delta = sort(rand(1,p));
%   coef = 1;
%   t    = linspace(-10,10,201);
%   MAX  = 0;
%   cf1  = cf_LogRV_WilksLambdaNC(t,p,n,q,delta,coef,MAX);
%   figure; plot(t,real(cf1),t,imag(cf1)); grid on; hold on;
%   MAX  = 20;
%   cf2  = cf_LogRV_WilksLambdaNC(t,p,n,q,delta,coef,MAX);
%   plot(t,real(cf2),t,imag(cf2)); hold off;
%   title('CF of log of non-central Wilks Lambda RV')
%
% EXAMPLE 2:
% % CF of a weighted linear combination of minus log Wilks Lambda RVs 
%   p     = [5 5 5];
%   n     = [10 15 20];
%   q     = [3 2 1];
%   delta = {sort(rand(1,p(1))),sort(rand(1,p(2))),sort(rand(1,p(3)))};
%   coef  = -[10 15 20]/45;
%   t     = linspace(-20,20,201);
%   cf    = cf_LogRV_WilksLambdaNC(t,p,n,q,delta,coef);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of a weighted linear combination of -log Wilks Lambda RVs')
%
% EXAMPLE 3: 
% % PDF/CDF of minus log Wilks Lambda RV (p=10, n=30, q=5) from its CF
%   p      = 10;
%   n      = 30;
%   q      = 5;
%   delta  = [1 2 3 10 30];
%   coef   = -1;
%   cf0    = @(t) cf_LogRV_WilksLambdaNC(t,p,n,q,[],coef);
%   cf     = @(t) cf_LogRV_WilksLambdaNC(t,p,n,q,delta,coef);
%   prob   = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result0 = cf2DistGP(cf0,[],prob,options);
%   figure
%   result  = cf2DistGP(cf,[],prob,options);
%   disp(result)
%   figure
%   plot(result0.x,result0.cdf,result.x,result.cdf);grid on
%   xlabel('x')
%   ylabel('CDF')
%   title('CDFs of -log(\Lambda) under null and alternative hypothesis')
%   figure
%   plot(result0.x,result0.pdf,result.x,result.pdf);grid on
%   xlabel('x')
%   ylabel('PDF')
%   title('PDFs of -log(\Lambda) under null and alternative hypothesis')
%
% EXAMPLE 4: (Compare exact vs. simulated non-central Wilk's distribution)
%   p = 10;             % p - length of the sample vectors (dimensionality)
%   n = 30;             % n - sample size / degrees of freedon
%   q = 5;              % q - degrees of freedom due to hypothesis model
%   N = 10000;          % N - number of simulation samples
%   M = rand(p,q);      % M - the (true) mean (p x q)-matrix
%   delta = eig(M*M');  % delta - the eigenvalues of non-centrality matrix
%   L = zeros(N,1);
%   for i = 1:N
%     X = randn(p,n);
%     E = X*X';
%     Y = randn(p,q) + M;
%     H = Y*Y';
%     L(i) = -log(det((E+H)\E));
%   end
%   % Exact and the empirical CDF of -log(L)
%   cf = @(t) cf_LogRV_WilksLambdaNC(t,p,n,q,delta,-1);
%   clear options
%   options.xMin = 0;
%   options.SixSigmaRule = 6;
%   prob = [0.9 0.95 0.99];
%   result = cf2DistGP(cf,[],prob,options);
%   hold on; ecdf(L); hold off
%   title('Exact vs. empirical CDF of the non-central distribution')
%
% REMARK:
%  Computing CF of the LOG-TRANSFORMED NON-CENTARL WILK's LAMBDA random
%  variable depends on computing the generalized hypergeometric function of
%  matrix argument. 
%  cf_LogRV_WilksLambdaNC uses modified implementation for computing the
%  truncated hypergeometric function, see the algorithm HypergeompFqMat.m,
%  as originaly suggested in Koev and Edelman (2006). The truncation of the
%  hypergeometric series is controled by the parameter MAX (start with MAX
%  = 20). If necessary, the parameter MAX should be set to sufficiently
%  large value but this can be numerically intractable. 
%  If MAX = 0 (now it is the default value), the algorithm uses the fast
%  (approximate) method for computing the confluent hypergeometric function
%  1F1(a;b;X), see Hypergeom1F1MatApprox, otherwise the algorithm uses the
%  algorithm Hypergeom1F1Mat based on computing the truncated expansion
%  series of 1F1(a;b;X) with MAX number of partitions.
%  The approximate alternative method (if MAX = 0) is based on heuristic
%  arguments (not formally proven), and the value of the confluent
%  hypergeometric function 1F1(a,b,X) of a matrix argument is calculated
%  (approximately) as 1F1(a;b;X) ~ 1F1(a;b;x(1)) * ... * 1F1(a;b;x(p)),
%  where 1F1(a;b;x(1)) is the scalar value confluent hypergeometric
%  function 1F1(a,b,x(i)) with [x(1),...,x(p)] = eig(X). 
%
%  By default (or if MAX is set to value MAX = 0), cf_LogRV_WilksLambdaNC
%  uses this alternative method for computing 1F1(a;b;X).
%
% REFERENCES:
% [1] Koev, P. and Edelman, A., 2006. The efficient evaluation of the
%     hypergeometric function of a matrix argument. Mathematics of
%     Computation, 75(254), 833-846.
% [2] Witkovský, V., 2018. Exact distribution of selected multivariate test
%     criteria by numerical inversion of their characteristic functions. 
%     arXiv preprint arXiv:1801.02248.
%
% SEE ALSO: cf_LogRV_WilksLambda, Hypergeom1F1Mat, Hypergeom1F1MatApprox

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 10-Aug-2018 15:46:49

%% ALGORITHM
% function cf = cf_LogRV_WilksLambdaNC(t,p,n,q,delta,coef,MAX)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 7);
if nargin < 7, MAX   = []; end
if nargin < 6, coef  = []; end
if nargin < 5, delta = []; end
if nargin < 4, q = []; end
if nargin < 3, n = []; end
if nargin < 2, p = []; end

%%
if isempty(p), p = 1; end
if isempty(n), n = 1; end
if isempty(q), q = 1; end
if isempty(delta), delta = cell(1); end
if isempty(coef),  coef  = 1; end
if isempty(MAX),   MAX   = 0; end

if ~iscell(delta)
    [p1,p2] =size(delta);
    if p1 == p2
        delta{1} = eig((delta+delta')/2);
    elseif min(p1,p2) == 1
        delta = mat2cell(delta(:)',1);
    elseif min(p1,p2) ~= 1
        error(message('InputSizeMismatch'));
    end
else
    for i = 1:length(delta)
        [p1,p2] =size(delta{i});
        if p1 == p2 && p1 > 0
            delta{i} = eig((delta{i}+delta{i}')/2);
        elseif min(p1,p2) > 1
            error(message('InputSizeMismatch'));
        end
    end
end

%% Check size of the parameters
[errorcode,p,n,q,coef] = distchck(4,p(:)',n(:)',q(:)',coef(:)');
if errorcode > 0
    error(message('InputSizeMismatch'));
end

if length(delta) == 1
    delta = repmat(delta,length(coef));
end

%% Characteristic function of a linear combination
szt = size(t);
t   = t(:);

cf  = 1;
for i = 1:length(coef)
    alpha = (n(i)+1-(1:p(i)))/2;
    beta  = q(i)/2;
    cf = cf .* cf_LogRV_Beta(coef(i)*t,alpha,beta);
    if ~isempty(delta{i})
        if MAX == 0
            cf = cf .* Hypergeom1F1MatApprox(1i*coef(i)*t, ...
                1i*coef(i)*t + (n(i)+q(i))/2, -delta{i}/2);
        elseif MAX > 0
            cf = cf .* Hypergeom1F1Mat(1i*coef(i)*t, ...
                1i*coef(i)*t + (n(i)+q(i))/2, -delta{i}/2,ceil(MAX));
        else
            cf = cf .* Hypergeom1F1Mat(1i*coef(i)*t, ...
                1i*coef(i)*t + (n(i)+q(i))/2, -delta{i}/2,0);
        end
    end
end
cf  = reshape(cf,szt);
cf(t==0) = 1;

end