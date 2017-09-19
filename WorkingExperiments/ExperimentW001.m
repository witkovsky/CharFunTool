%% EXPERIMENTW001
%  19-Sep-2017 18:15:50

% Gil-Pelaez integral calculated by the Fourier Integral (for non-negative
% distributions) 

clear 

%% LogNormal distribution

mu = 0;
sigma = 1;
cf = @(t)cfX_LogNormal(t,mu,sigma);
fun = @(t) imag(cf(t))./t;
A = 1e-12;
B = 50;
x = linspace(0,7,501);

[FI,coefs] = FourierIntegral_BV(x,fun,A,B);
CDF = 1-2*real(FI)/pi;
plot(x,CDF,'.-');grid
hold on
plot(x,logncdf(x,0,1))
hold off

%% ChiSquare distribution

nu = 3
cf = @(t)cf_ChiSquare(t,nu);
fun = @(t) imag(cf(t))./t;
A = 1e-12;
B = 50;
x = linspace(0,15,501);
nPts = 50;

[FI,coefs] = FourierIntegral_BV(x,fun,A,B,nPts);
CDF = 1-2*real(FI)/pi;
plot(x,CDF,'.-');grid
hold on
plot(x,chi2cdf(x,nu))
hold off

figure
plot(coefs,'o-');grid

%% Compound LogNormal-Poisson distribution

mu = 0;
sigma = 1;
lambda = 5;
cfX = @(t)cfX_LogNormal(t,mu,sigma);
cf = @(t) cfN_Poisson(t,lambda,cfX);
fun = @(t) imag(cf(t))./t;
A = 1e-12;
B = 50;
x = linspace(0,50,501);

[FI,coefs] = FourierIntegral_BV(x,fun,A,B);
CDF = 1-2*real(FI)/pi;
plot(x,CDF,'.-');grid
