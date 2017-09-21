%% EXPERIMENTW001
%  19-Sep-2017 18:15:50

% Gil-Pelaez integral (CDF) calculated by using the BAKHVALOV-VASILEVA
% method for computing the Fourier Integral (for non-negative
% distributions).
% The CDF is calculated by the GIL-PELAEZ formula:
%  F(x) =  2 * (1/2 - 1/pi * Int_0^Inf exp(i*t*x) * Imag(cf(t)/t) dt)
%       =  1 - 2/pi * Int_0^Inf exp(i*t*x) * Imag(cf(t)/t) dt
% 

clear 

%% LogNormal distribution

mu = 0;
b = 1;
cf = @(t)cfX_LogNormal(t,mu,b);
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

figure
plot(coefs,'o-');grid

%% ChiSquare distribution

nu = 3;
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

%% Weibull distribution

a = 1.5;
b = 1;
cf = @(t)cfX_Weibull(t,a,b);
fun = @(t) imag(cf(t))./t;
A = 1e-13;
B = 70;
nPts = 30;
x = linspace(0,5,501);

[FI,coefs] = FourierIntegral_BV(x,fun,A,B,nPts);
CDF = 1-2*real(FI)/pi;
plot(x,CDF,'.-');grid
hold on
plot(x,wblcdf(x,1/a,b))
hold off

figure
plot(coefs,'o-');grid


%% Compound LogNormal-Poisson distribution

mu = 0;
b = 1;
lambda = 5;
cfX = @(t)cfX_LogNormal(t,mu,b);
cf = @(t) cfN_Poisson(t,lambda,cfX);
fun = @(t) imag(cf(t))./t;
A = 1e-12;
B = 50;
x = linspace(0,50,501);

[FI,coefs] = FourierIntegral_BV(x,fun,A,B);
CDF = 1-2*real(FI)/pi;
plot(x,CDF,'.-');grid

figure
plot(coefs,'o-');grid

%% Compound Weibull-Poisson distribution

a = 1.5;
b = 1;
lambda = 5;
cfX = @(t)cfX_Weibull(t,a,b);
cf = @(t) cfN_Poisson(t,lambda,cfX);
fun = @(t) imag(cf(t))./t;
A = 1e-12;
B = 50;
nPts = 75;
x = linspace(0,12,501);

tic;[FI,coefs] = FourierIntegral_BV(x,fun,A,B,nPts);toc
CDF = 1-2*real(FI)/pi;
plot(x,CDF,'.-');grid

figure
plot(coefs,'o-');grid

%% LogNormal distribution

mu = 0;
b = 1;
cf = @(t)cfX_LogNormal(t,mu,b);
fun = @(t) imag(cf(t)./t);
x = linspace(0,5,1001);
nMax = 150;
tol = 1e-12;

A = 1e-15; B = 0.1;
[coefs,scale,shift,xNodes,w,P] = LegendreSeries(fun,A,B,nMax,tol);
FI1 = LegendreSeriesFourierIntegral(x,coefs,scale,shift);

A = B; B = 1;
[coefs,scale,shift] = LegendreSeries(fun,A,B,nMax,tol,xNodes,w,P);
FI2 = LegendreSeriesFourierIntegral(x,coefs,scale,shift);

A = B; B = 10;
[coefs,scale,shift] = LegendreSeries(fun,A,B,nMax,tol,xNodes,w,P);
FI3 = LegendreSeriesFourierIntegral(x,coefs,scale,shift);

A = B; B = 100;
[coefs,scale,shift] = LegendreSeries(fun,A,B,nMax,tol,xNodes,w,P);
FI4 = LegendreSeriesFourierIntegral(x,coefs,scale,shift);

A = B; B = 1000;
[coefs,scale,shift] = LegendreSeries(fun,A,B,nMax,tol,xNodes,w,P);
FI5 = LegendreSeriesFourierIntegral(x,coefs,scale,shift);

A = B; B = 10000;
[coefs,scale,shift] = LegendreSeries(fun,A,B,nMax,tol,xNodes,w,P);
FI6 = LegendreSeriesFourierIntegral(x,coefs,scale,shift);

A = B; B = 1e+300;
[coefs,scale,shift] = LegendreSeries(fun,A,B,nMax,tol,xNodes,w,P);
FI7 = LegendreSeriesFourierIntegral(x,coefs,scale,shift);

FI = FI1 + FI2 + FI3 + FI4 + FI5 + FI6 + FI7;

CDF = 1-2*real(FI)/pi;
figure
plot(x,CDF,'.-');grid

hold on
plot(x,logncdf(x,0,1))
hold off