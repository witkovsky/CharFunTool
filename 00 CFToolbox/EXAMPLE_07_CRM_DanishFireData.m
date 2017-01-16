%% EMPIRICAL CRM DISTRIBUTION (The Danish Insurance Data)
%  (c) 2016 Viktor Witkovsky, witkovsky@gmail.com
%  Version: 15-Nov-2016

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
title('Histogram of Danish insurance data: Fire losses from 1980 to 1990')
xlabel('Claim size (in millions DKK)')
ylabel('Count')

%% *Empirical Compound Distribution*

%% Empirical Characteristic Function (ECF) of the severity distribution

cfX = @(t) cfE_DiracMixture(t,Severity); 

t = linspace(-25,25,501);
figure
plot(t, real(cfX(t)),t,imag(cfX(t)))
axis([-25,25,-1,1]);
title('ECF of the severity distribution')
xlabel('t')
ylabel('CF')
grid on

%% Empirical Characteristic Function (ECF) of the frequency distribution

cfN = @(t) cfE_DiracMixture(t,Frequency);

t = linspace(-0.5,0.5,501);
figure
plot(t, real(cfN(t)),t,imag(cfN(t)))
axis([-0.5,0.5,-1,1]);
title('ECF of the frequency distribution')
xlabel('t')
ylabel('CF')
grid on

%% Empirical Characteristic Function (ECF)of the compound distribution

cf   = @(t) cfN(-1i*log(cfX(t)));

t = linspace(-.05,.05,501);
figure
plot(t, real(cf(t)),t,imag(cf(t)))
axis([-.05,.05,-1,1]);
title('ECF of the compound distribution')
xlabel('t')
ylabel('CF')
grid on

%% PDF/CDF by numerical inversion of the compound ECF
%  Empirical Aggregate Loss Distribution / Collective Risk Distribution of
%  the Danish Insurance Data

prob = [0.9 0.95 0.99];
x    = linspace(0,2000,201)';
clear options
options.isCompound = true;
result = cf2DistGP(cf,x,prob,options);
disp(result)

%% PLOT CRM Data Matrix

% Histogram of insurance data: Danish fire losses 1980 -- 1990
axes1 = subplot(2,3,1);
hist(log10(Severity),50);
xlim(axes1,[0 3]);
set(axes1,'CLim',[1 3],'XGrid','on',...
    'XTick',[0 1  2 3],...
    'XTickLabel',{'1','10','100','1000'});
title('Insurance data: Danish fire losses 1980 -- 1990')
xlabel('Claim loss (millions DKK)')
ylabel('count')

% Aggregate loss ditribution (PDF)
axes2 = subplot(2,3,2);
plot(x,result.pdf,'Linewidth',2)
grid on
title('Aggregate loss distribution')
xlabel('Aggregate loss (in millions DKK)')
ylabel('PDF')

% Aggregate loss distribution (CDF)
axes3 = subplot(2,3,3);
plot(x,result.cdf,'Linewidth',2);
grid on
title('Aggregate loss distribution')
xlabel('Aggregate loss (in millions DKK)')
ylabel('CDF')

% ECF of the severity distribution
axes4 = subplot(2,3,4);
t = linspace(-25,25,501);
plot(t, real(cfX(t)),t,imag(cfX(t)))
axis([-25,25,-1,1]);
title('ECF of the severity distribution')
xlabel('t')
ylabel('CF')
set(axes4,'YGrid','on');

% ECF of the frequency distribution
axes5 = subplot(2,3,5);
t = linspace(-0.5,0.5,501);
plot(t, real(cfN(t)),t,imag(cfN(t)))
axis([-0.5,0.5,-1,1]);
title('ECF of the frequency distribution')
xlabel('t')
ylabel('CF')
set(axes5,'YGrid','on');

% ECF of the compound distribution
axes6 = subplot(2,3,6);
t = linspace(-.05,.05,501);
plot(t, real(cf(t)),t,imag(cf(t)))
axis([-.05,.05,-1,1]);
title('ECF of the compound distribution')
xlabel('t')
ylabel('CF')
set(axes6,'YGrid','on');
