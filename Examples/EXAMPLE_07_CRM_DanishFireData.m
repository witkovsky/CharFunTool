%% CharFunTool / EXAMPLES
%  CHARACTERISTIC FUNCTIONS OF THE EMPIRICAL COLLECTIVE RISK DISTRIBUTIONS
%
%  (c) 2017 Viktor Witkovsky, witkovsky@gmail.com
%  Version: 10-May-2017

clear
close all

%% EMPIRICAL CRM DISTRIBUTION (The Danish Insurance Data)

% Data by Alexander McNeil: http://www.macs.hw.ac.uk/~mcneil/data.html
%
% The dataset is comprised of Danish fire losses analyzed in McNeil (1997).
% These data represent Danish fire losses in million Danish Krones and were
% collected by a Danish reinsurance company. The dataset contains
% individual losses above 1 million Danish Krones, a total of 2167
% individual losses, covering the period from January 3, 1980 to December
% 31, 1990. The data is adjusted for inflation to reflect 1985 values. The
% dataset can be found in the R packages fEcofin and fExtremes. These data
% are considered in Resnick (1997), Cooray and Ananda (2005), and
% DellAquila and Embrechts (2006), among others.
%

clear
close all

%% The Danish Insurance Data
load('DanishFireData.mat')

% Histogram
figure1 = figure;
axes1 = axes('Parent',figure1);
hist(log10(Severity),50);
xlim(axes1,[0 3]);
set(axes1,'CLim',[1 3],'XGrid','on',...
    'XTick',[0  1 2  3],...
    'XTickLabel',{'1','10','100','1000'});
title('Histogram of Danish Fire Losses 1980-1990')
xlabel('individual losses (millions DKK)')
ylabel('count')

%% *Empirical Compound Distribution*

%% Empirical Characteristic Function (ECF) of the severity distribution

cfX = @(t) cfE_Empirical(t,Severity); 

t = linspace(-25,25,501);
figure
plot(t, real(cfX(t)),t,imag(cfX(t)))
axis([-25,25,-1,1]);
title('Empirical CF of Severity Distribution')
xlabel('t')
ylabel('CF')
grid on

%% Empirical Characteristic Function (ECF) of the frequency distribution

cfN = @(t) cfE_Empirical(t,Frequency);

t = linspace(-0.5,0.5,501);
figure
plot(t, real(cfN(t)),t,imag(cfN(t)))
axis([-0.5,0.5,-1,1]);
title('Empirical CF of Frequency Distribution')
xlabel('t')
ylabel('CF')
grid on

%% Empirical Characteristic Function (ECF)of the compound distribution

cf   = @(t) cfN(-1i*log(cfX(t)));

t = linspace(-.05,.05,501);
figure
plot(t, real(cf(t)),t,imag(cf(t)))
axis([-.05,.05,-1,1]);
title('Empirical CF of Compound Distribution')
xlabel('t')
ylabel('CF')
grid on

%% PDF/CDF by numerical inversion of the compound ECF
%  Empirical Aggregate Loss Distribution / Collective Risk Distribution of
%  the Danish Insurance Data

prob = [0.9 0.95 0.99];
loss    = linspace(0,2000,201)';
clear options
options.isCompound = true;
result = cf2DistGP(cf,loss,prob,options);
disp(result)

%% PLOT CRM Data Matrix

% Histogram of insurance data: Danish fire losses 1980 -- 1990
axes1 = subplot(2,3,1);
hist(log10(Severity),50);
xlim(axes1,[0 3]);
set(axes1,'CLim',[1 3],'XGrid','on',...
    'XTick',[0 1  2 3],...
    'XTickLabel',{'1','10','100','1000'});
title('Histogram of Danish Fire Losses 1980-1990')
xlabel('individual loss (millions DKK)')
ylabel('count')

% Aggregate loss ditribution (PDF)
axes2 = subplot(2,3,2);
plot(loss,result.pdf,'Linewidth',2)
grid on
title('PDF of Aggregate Loss Distribution')
xlabel('aggregate loss (millions DKK)')
ylabel('PDF')

% Aggregate loss distribution (CDF)
axes3 = subplot(2,3,3);
plot(loss,result.cdf,'Linewidth',2);
grid on
title('CDF of Aggregate Loss Distribution')
xlabel('aggregate loss (millions DKK)')
ylabel('CDF')

% ECF of the severity distribution
axes4 = subplot(2,3,4);
t = linspace(-25,25,501);
plot(t, real(cfX(t)),t,imag(cfX(t)))
axis([-25,25,-1,1]);
title('Empirical CF of Severity Distribution')
xlabel('t')
ylabel('CF')
set(axes4,'YGrid','on');

% ECF of the frequency distribution
axes5 = subplot(2,3,5);
t = linspace(-0.5,0.5,501);
plot(t, real(cfN(t)),t,imag(cfN(t)))
axis([-0.5,0.5,-1,1]);
title('Empirical CF of Frequency Distribution')
xlabel('t')
ylabel('CF')
set(axes5,'YGrid','on');

% ECF of the compound distribution
axes6 = subplot(2,3,6);
t = linspace(-.05,.05,501);
plot(t, real(cf(t)),t,imag(cf(t)))
axis([-.05,.05,-1,1]);
title('Empirical CF of Compound Distribution')
xlabel('t')
ylabel('CF')
set(axes6,'YGrid','on');

%% AGGREGATE LOSS DISTRIBUTION with CF TOOLBOX

% Danish fire losses data
load('DanishFireData.mat')

% Empirical characteristic functions               
cfN  = @(t) cfE_Empirical(t,Frequency);
cfX  = @(t) cfE_Empirical(t,Severity);
cf   = @(t) cfN(-1i*log(cfX(t)));

% Parameters/options
prob = [0.9 0.99 0.999];
loss = linspace(0,2000,201)';
options.isCompound = true;

% Numerical inversion of CF by cf2DistGP
result = cf2DistGP(cf,loss,prob,options); 
disp(result)

%% Semi-parametric approach with fitted generalized Pareto tail distribution

% Danish fire losses data
load('DanishFireData.mat')

% Set the threshold parameter theta
p      = 0.95;
theta  = quantile(Severity,p);

% Fit the GP (Generalized Pareto) distribution
GPfit  = paretotails(Severity,0,p);
Pars   = GPfit.UpperParameters;
xi     = Pars(1);
sigma  = Pars(2);

% CF of the fitted tail GP distribution 
pdfGP  = @(x) gppdf(x,xi,sigma);
method = 'fit';
cfGP   = @(t) cfX_PDF(t,pdfGP,method) .* exp(1i*t*theta);

% CF of the mixture severity distribution
XL     = Severity(Severity <= theta);
% Smoothed empirical CF
% cfXL  = @(t) cfE_Empirical(t,XL) .* cfS_Gaussian(t*0.2);
cfXL   = @(t) cfE_Empirical(t,XL);
cfX    = @(t) p * cfXL(t) + (1-p) * cfGP(t);

% Empirical CF of the frequency distribution
cfN    = @(t) cfE_Empirical(t,Frequency);

% Compound CF of the aggregate loss distribution
cf     = @(t) cfN(-1i*log(cfX(t)));

% Parameters
prob   = [0.9 0.99 0.999];
loss   = linspace(0,2500,201)';

% Options
clear options
options.N = 2^15;
options.SixSigmaRule = 15;
options.isCompound = true;

% Numerical inversion of CF by cf2DistGP
result = cf2DistGP(cf,loss,prob,options); 

%% Plot the combined compound CF

t = linspace(-.05,.05,501);
figure
plot(t, real(cf(t)),t,imag(cf(t)))
axis([-.05,.05,-1,1]);
title('CF of Compound Distribution')
xlabel('t')
ylabel('CF')
grid on