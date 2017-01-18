%% EXPERIMETS with APPROXIMATE CFs
% Approximate CF by weighted mixture ef Dirac CFS
% cf(t) ~ exp(1i * t * x' ) * w
% - x is vector of selected nodes from the domain of X
% - w is vector of (optimum integration) weights for selected nodes x
% Any optimum quadrature rule can be used (Gauss, Clenshaw-Curtis, ...)
% Here we illustrate using chebyshew poits and Clenshaw-Curtis weights

clear

%% Standard Normal (Gaussian) distribution
N = 2^15;

% Chebpoints - nodes and weights of the Clenshaw-Curtis quadrature rule
[s,w] = chebpts(N);

% probability p-nodes and weights
p = (s+1)/2; 
p = p(2:end-1);
p = p(:);
w = w(2:end-1);
w = w(:)/2;

% Gauss quadrature
%[p,w] = grule(N);


% x-nodes from probability p-nodes for given distribution (requires
% quantile function) 
x = norminv(p);

% Approximate characteristic function:
cfApprox = @(t) real(reshape(exp(1i*t(:)*x')*w,size(t)));

% Here we know the true characteristic function:
cfTrue = @(t) cfS_Gaussian(t);

%% Plot and compare the True and the Approximate CFs
t = linspace(-10,10,1001)';

figure
plot(t,real(cfTrue(t)),t,imag(cfTrue(t)));
hold on
plot(t,real(cfApprox(t)),t,imag(cfApprox(t)));
grid
hold off

figure
plot(t,abs(cfApprox(t)-cfTrue(t)))

%% PDF / CDF by numerical inversion of the APPROXIMATE CF

clear options
options.N = 2^10;
 options.tolDiff = 1e-4; 
result = cf2DistGP(cfApprox,[],[],options);
disp(result)

%% Distribution of a transformed random variable f(X)

fun = @(x) x.^4 - 2*x.^2;
xFun = fun(x);
cfxFun = @(t) reshape(exp(1i*t(:)*fun(x'))*w,size(t));

figure
plot(t,real(cfxFun(t)),t,imag(cfxFun(t)));


clear options
options.N = 2^10;
options.xMin = 0;
options.SixSigmaRule = 15;
result = cf2DistGP(cfxFun,[],[],options);
disp(result)
