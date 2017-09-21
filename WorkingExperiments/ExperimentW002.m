%% EXPERIMENTW002
%  20-Sep-2017 10:12:13

% Gil-Pelaez integral (PDF) calculated by using the BAKHVALOV-VASILEVA
% method for computing the Fourier Integral (for non-negative
% distributions).
% The CDF is calculated by the GIL-PELAEZ formula:
%  f(x) =  2 * (1/pi * Int_0^Inf exp(i*t*x) * Real(cf(t)) dt)
%       =  2/pi * Int_0^Inf exp(i*t*x) * Real(cf(t)) dt
% 

clear 
%% LogNormal distribution

mu = 0;
b = 1;
cf = @(t)cfX_LogNormal(t,mu,b);
fun = @(t) real(cf(t));
A = 1e-12;
B = 50;
xNodes = linspace(0,7,501);

[FI,coefs] = FourierIntegral_BV(xNodes,fun,A,B);
PDF = 2*real(FI)/pi;
plot(xNodes,PDF,'.-');grid
hold on
plot(xNodes,lognpdf(xNodes,0,1))
hold off

figure
plot(coefs,'o-');grid

%% ChiSquare distribution

nu = 2;
cf = @(t)cf_ChiSquare(t,nu);
fun = @(t) real(cf(t));
A = 1e-12;
B = 500;
xNodes = linspace(0,15,501);
nPts = 150;

[FI,coefs] = FourierIntegral_BV(xNodes,fun,A,B,nPts);
PDF = 2*real(FI)/pi;
plot(xNodes,PDF,'.-');grid
hold on
plot(xNodes,chi2pdf(xNodes,nu))
hold off

figure
plot(coefs,'o-');grid
%% Weibull distribution

a = 1.5;
b = 1;
cf = @(t)cfX_Weibull(t,a,b);
fun = @(t) real(cf(t));
A = 1e-13;
B = 70;
nPts = 70;
xNodes = linspace(0,5,501);

[FI,coefs] = FourierIntegral_BV(xNodes,fun,A,B,nPts);
PDF = 2*real(FI)/pi;
plot(xNodes,PDF,'.-');grid
hold on
plot(xNodes,wblpdf(xNodes,1/a,b))
hold off

figure
plot(coefs,'o-');grid

%% Compound LogNormal-Poisson distribution

mu = 0;
b = 1;
lambda = 5;
cfX = @(t)cfX_LogNormal(t,mu,b);
cf = @(t) cfN_Poisson(t,lambda,cfX);
fun = @(t) real(cf(t));
A = -1e-12;
B = 2.5;
xNodes = linspace(0,50,501);
nPts = 150;

[FI,coefs] = FourierIntegral_BV(xNodes,fun,A,B,nPts);
PDF = 2*real(FI)/pi;
plot(xNodes,PDF,'.-');grid

figure

plot(coefs,'o-');grid

%% Compound Weibull-Poisson distribution

a = 1.5;
b = 1;
lambda = 5;
cfX = @(t)cfX_Weibull(t,a,b);
cf = @(t) cfN_Poisson(t,lambda,cfX);
fun = @(t) real(cf(t));
A = 1e-12;
B = 5;
nPts = 75;
xNodes = linspace(0,15,501);

tic;[FI,coefs] = FourierIntegral_BV(xNodes,fun,A,B,nPts);toc
PDF = 2*real(FI)/pi;
plot(xNodes,PDF,'.-');grid

figure
plot(coefs,'o-');grid

%% LogNormal distribution

mu = 0;
b = 1;
cf = @(t)cfX_LogNormal(t,mu,b);
fun = @(t) real(cf(t));
xNodes = linspace(0,5,501);
nPts = 10;

A1 = 1e-12;
B1 = 1;
[FI1,coefs1] = FourierIntegral_BV(xNodes,fun,A1,B1,nPts);
plot(coefs1,'o-');grid

A2 = 1;
B2 = 10;
[FI2,coefs2] = FourierIntegral_BV(xNodes,fun,A2,B2,nPts);
figure
plot(coefs2,'o-');grid

A3 = 10;
B3 = 100;
[FI3,coefs3] = FourierIntegral_BV(xNodes,fun,A3,B3,nPts);
figure
plot(coefs3,'o-');grid

A4 = 100;
B4 = 1000;
[FI4,coefs4] = FourierIntegral_BV(xNodes,fun,A4,B4,nPts);
figure
plot(coefs4,'o-');grid

A5 = 1000;
B5 = 10000;
[FI5,coefs5] = FourierIntegral_BV(xNodes,fun,A5,B5,nPts);
figure
plot(coefs5,'o-');grid

FI = FI1 + FI2 + FI3 + FI4 + FI5;

PDF = 2*real(FI)/pi;
figure
plot(xNodes,PDF,'.-');grid
hold on
plot(xNodes,lognpdf(xNodes,0,1))
hold off

%% LogNormal distribution

mu = 0;
b = 1;
cf = @(t)cfX_LogNormal(t,mu,b);
fun = @(t) real(cf(t));
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

PDF = 2*real(FI)/pi;
figure
plot(x,PDF,'.-');grid

hold on
plot(x,lognpdf(x,0,1))
hold off