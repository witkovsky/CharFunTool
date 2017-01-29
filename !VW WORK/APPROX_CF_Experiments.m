%% EXPERIMETS with APPROXIMATE CFs
% Approximate CF by weighted mixture ef Dirac CFS
% cf(t) ~ exp(1i * t * x' ) * w
% - x is vector of selected nodes from the domain of X
% - w is vector of (optimum integration) weights for selected nodes x
% Any optimum quadrature rule can be used (Gauss, Clenshaw-Curtis, ...)
% Here we illustrate using chebyshew poits and Clenshaw-Curtis weights

clear

%% Standard Normal (Gaussian) distribution
nPts = 2^10;

% Chebpoints - nodes and weights of the Clenshaw-Curtis quadrature rule
precission = eps(1);
[p,w] = chebpts(nPts,[precission/2,1-precission/2]);

% Gauss quadrature
%[p,w] = grule(N);

% x-nodes from probability p-nodes for given distribution (requires
% quantile function) 
x = norminv(p(:));

% Approximate characteristic function:
cfApprox = @(t) real(reshape(exp(1i*t(:)*x')*w(:),size(t)));

% Here we know the true characteristic function:
cfTrue = @(t) cfS_Gaussian(t);

%% Plot and compare the True and the Approximate CFs
t = linspace(-30,30,1001)';

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
options.N = 2^8;
options.tolDiff = 1e-4;
options.isInterp = true;
options.xMin = min(x);
options.xMax = max(x);

result = cf2DistGP(cfApprox,[],[],options);
disp(result)

qf = result.QF;

%% Distribution of a transformed random variable f(X)

fun = @(x) x.^4 - 2*x.^2;
xFun = fun(x);
cfxFun = @(t) reshape(exp(1i*t(:)*fun(x'))*w(:),size(t));

figure
plot(t,real(cfxFun(t)),t,imag(cfxFun(t)));


clear options
options.N = 2^10;
options.xMin = 0;
options.SixSigmaRule = 15;
result = cf2DistGP(cfxFun,[],[],options);
disp(result)


%% Compound Log-normal-Poisson distribution

clear

%% LogNormal distribution  with parameters (0,2)
nPts = 2^10;

% Chebpoints - nodes and weights of the Clenshaw-Curtis quadrature rule
precission = eps(1);
[p,w] = chebpts(nPts,[precission/2,1-precission/2]);

% x-nodes from probability p-nodes for given distribution (requires
% quantile function) 
x = logninv(p(:),0,2);

% Approximate characteristic function:
cfX_LN = @(t) reshape(exp(1i*t(:)*x')*w(:),size(t));

% Compound Poisson with lambda = 1000
lambda = 1000;
cf = @(t) cfN_Poisson(t,lambda,cfX_LN);

%% Invert PDF/CDF/QF
clear options
options.N = 2^10;
options.isInterp = true;
options.isCompound = true;
options.SixSigmaRule = 18;
prob = 0.999;

result = cf2DistGP(cf,[],prob,options);
disp(result)

qf = result.QF;q = qf(prob);
disp(q)

%% Sinc Expansion
N = 2^10;
m = -N:N;
h = result.dt;
a = 1e-8;

t = m*h +1i*a;

%plot(real(t),real(cf(t)),real(t),imag(cf(t)))

%fSincExp = @(z) sum(cf(m*h + 1i*a) .* sin(pi*(z-(m*h+1i*a))/h) ./ (pi*(z-(m*h+1i*a))/h));

prob = 0.95;
plot(real(t),real(exp(-1i*t*qf(prob)).*cf(t)),real(t),imag(exp(-1i*t*qf(prob)).*cf(t)))

Ipdf = real(sum(exp(-1i*t*qf(prob)).*cf(t)*h)/(2*pi))

%%

plot(real(t),real((exp(-1i*t*qf(prob)).*cf(t) .* (1-(-1).^m*exp(-pi*a*h)) ./ (pi*t/h))/2),real(t),imag((exp(-1i*t*qf(prob)).*cf(t) .* (1-(-1).^m*exp(-pi*a/h)) ./ (pi*t/h))/2))

Icdf = -imag( sum(exp(-1i*t*qf(prob)).*cf(t) .* (1-(-1).^m*exp(-pi*a/h)) ./ (pi*t/h))/2)

