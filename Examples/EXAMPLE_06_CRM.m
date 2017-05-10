%% CharFunTool / EXAMPLES
%  CHARACTERISTIC FUNCTIONS OF THE COLLECTIVE RISK DISTRIBUTIONS
%
%  (c) 2017 Viktor Witkovsky, witkovsky@gmail.com
%  Version: 10-May-2017

clear
close all

%% Collective Risk Model (CRM) Distributions

% A typical model for insurance risk is the Collective Risk Model (CRM).
% The CRM mathematically describes the aggregate loss of an insurance
% portfolio in a certain period of time (e.g. 1 year). Insurance portfolio
% is regarded as a process that produces claims over time. The size of
% these claims (severity) are iid (independent identically distributed)
% random variables, independent also of the number of claims generated
% (frequency) in this time period. Here we illustrate method of calculation
% based on the numerical inversion of the associated CRM’s compound
% characteristic function of the frequency distribution and the severity
% distribution. The used numerical algorithm (cf2DistGP) is based on the
% Gil-Pelaez inversion formulae for computing the probability distribution
% functions (PDF and/or CDF) of the univariate continuous random variables.

% Here, the sizes of the claims are taken to be iid  random variables
% X1,...,XN, independent also of the number of claims N generated in this
% time period. The (unknown) amount Xi, i = 1,...,N, is called single-loss
% random variable or severity and N is called claim count random variable
% or frequency. This makes the total claim the sum of a random number N of
% iid individual random claims X1,...,XN: S = X1 + ... + XN.
%
% NOTE
% The examples can be calculated easily by using the CRMTool application.
% CRMTool is an calculator suggested for calculation of probability density
% functions (PDF) and cumulative distribution functions (CDF) of the
% aggregate loss distributions specified by the CRMs compound
% characteristic functions of frequency distribution and the severity
% distribution. 
%
% REFERENCE
% Witkovsky, V., WIMMER, G., DUBY, T. (2017). Computing the aggregate loss
% distribution based on numerical inversion of the compound empirical
% characteristic function of frequency and severity. Working Paper. 

clear 
close all

% Set the auxiliary values
x    = (0:0.1:5)';
x    = [x;(5.5:0.5:15)'];
x    = [x;(15.5:0.5:50)'];
prob = [0.9 0.95 0.99];

format shorte

%% *COMPOUND BINOMIAL DISTIBUTIONS*

%% CDF/PDF Compound Binomial-Pareto1 distribution

alpha = 2;  
cfX   = @(t) cfX_Pareto(t,alpha);

n     = 4;
p     = 0.75;
cf    = @(t) cfN_Binomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoBinoPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoBinoPar1.x CoBinoPar1.pdf CoBinoPar1.cdf])
disp(CoBinoPar1)

%% CDF/PDF Compound Binomial-ParetoA distribution

alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

n     = 3;
p     = 0.75;
cf    = @(t) cfN_Binomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoBinoParA = cf2DistGP(cf,x,prob,options);

%disp([CoBinoParA.x CoBinoParA.pdf CoBinoParA.cdf])
disp(CoBinoParA)

%% CDF/PDF Compound Binomial-ParetoE distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

n     = 3;
p     = 0.75;
cf    = @(t) cfN_Binomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoBinoParE = cf2DistGP(cf,x,prob,options);

%disp([CoBinoParE.x CoBinoParE.pdf CoBinoParE.cdf])
disp(CoBinoParE)

%% CDF/PDF Compound Binomial-PearsonVI distribution
alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

n     = 3;
p     = 0.75;
cf    = @(t) cfN_Binomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoBinoPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoBinoPeaVI.x CoBinoPeaVI.pdf CoBinoPeaVI.cdf])
disp(CoBinoPeaVI)

%% *COMPOUND DELAPORTE DISTIBUTIONS*

%% CDF/PDF Compound Delaporte-Pareto1 distribution
alpha = 2;
cfX   = @(t) cfX_Pareto(t,alpha);

b     = 2.2;
a     = 3.3;
c     = 4;
cf    = @(t) cfN_Delaporte(t,a,b,c,cfX);

clear options
options.isCompound = 1;
CoDelaPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoDelaPar1.x CoDelaPar1.pdf CoDelaPar1.cdf])
disp(CoDelaPar1)

%% CDF/PDF Compound Delaporte-ParetoA distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

b     = 2.2;
a     = 3.3;
c     = 4;
cf    = @(t) cfN_Delaporte(t,a,b,c,cfX);

clear options
options.isCompound = 1;
CoDelaParA = cf2DistGP(cf,x,prob,options);

%disp([CoDelaParA.x CoDelaParA.pdf CoDelaParA.cdf])
disp(CoDelaParA)

%% CDF/PDF Compound Delaporte-ParetoE distribution

alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

b     = 2.2;
a     = 3.3;
c     = 4;
cf    = @(t) cfN_Delaporte(t,a,b,c,cfX);

clear options
options.isCompound = 1;
CoDelaParE = cf2DistGP(cf,x,prob,options);

%disp([CoDelaParE.x CoDelaParE.pdf CoDelaParE.cdf])
disp(CoDelaParE)

%% CDF/PDF Compound Delaporte-PearsonVI distribution

alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

b     = 2.2;
a     = 3.3;
c     = 4;
cf    = @(t) cfN_Delaporte(t,a,b,c,cfX);

clear options
options.isCompound = 1;
CoDelaPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoDelaPeaVI.x CoDelaPeaVI.pdf CoDelaPeaVI.cdf])
disp(CoDelaPeaVI)

%% *COMPOUND GENERALIZED POISSON DISTIBUTIONS*

%% CDF/PDF Compound GeneralizedPoisson-Pareto1 distribution

alpha = 2;
cfX   = @(t) cfX_Pareto(t,alpha);

a     = 2.2;
p     = 0.75;
cf    = @(t) cfN_GeneralizedPoisson(t,a,p,cfX);

clear options
options.isCompound = 1;
CoGePoPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoGePoPar1.x CoGePoPar1.pdf CoGePoPar1.cdf])
disp(CoGePoPar1)

%% CDF/PDF Compound GeneralizedPoisson-ParetoA distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2.2;
p     = 0.75;
cf    = @(t) cfN_GeneralizedPoisson(t,a,p,cfX);

clear options
options.isCompound = 1;
CoGePoParA = cf2DistGP(cf,x,prob,options);

%disp([CoGePoParA.x CoGePoParA.pdf CoGePoParA.cdf])
disp(CoGePoParA)

%% CDF/PDF Compound GeneralizedPoisson-ParetoE distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2.2;
p     = 0.75;
cf    = @(t) cfN_GeneralizedPoisson(t,a,p,cfX);

clear options
options.isCompound = 1;
CoGePoParE = cf2DistGP(cf,x,prob,options);

%disp([CoGePoParE.x CoGePoParE.pdf CoGePoParE.cdf])
disp(CoGePoParE)

%% CDF/PDF Compound GeneralizedPoisson-PearsonVI distribution

alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

a     = 2.2;
p     = 0.75;
cf    = @(t) cfN_GeneralizedPoisson(t,a,p,cfX);

clear options
options.isCompound = 1;
CoGePoPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoGePoPeaVI.x CoGePoPeaVI.pdf CoGePoPeaVI.cdf])
disp(CoGePoPeaVI)

%% *COMPOUND NEGATIVE BINOMIAL DISTIBUTIONS*

%% CDF/PDF Compound NegativeBinomial-Pareto1 distribution

alpha = 2;
cfX   = @(t) cfX_Pareto(t,alpha);

n     = 4;
p     = 0.75;
cf    = @(t) cfN_NegativeBinomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoNeBiPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoNeBiPar1.x CoNeBiPar1.pdf CoNeBiPar1.cdf])
disp(CoNeBiPar1)

%% CDF/PDF Compound NegativeBinomial-ParetoA distribution

alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

n     = 4;
p     = 0.75;
cf    = @(t) cfN_NegativeBinomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoNeBiParA = cf2DistGP(cf,x,prob,options);

%disp([CoNeBiParA.x CoNeBiParA.pdf CoNeBiParA.cdf])
disp(CoNeBiParA)

%% CDF/PDF Compound NegativeBinomial-ParetoE distribution

alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

n     = 4;
p     = 0.75;
cf    = @(t) cfN_NegativeBinomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoNeBiParE = cf2DistGP(cf,x,prob,options);

%disp([CoNeBiParE.x CoNeBiParE.pdf CoNeBiParE.cdf])
disp(CoNeBiParE)

%% CDF/PDF Compound NegativeBinomial-PearsonVI distribution

alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

n     = 4;
p     = 0.75;
cf    = @(t) cfN_NegativeBinomial(t,n,p,cfX);

clear options
options.isCompound = 1;
CoNeBiPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoNeBiPeaVI.x CoNeBiPeaVI.pdf CoNeBiPeaVI.cdf])
disp(CoNeBiPeaVI)

%% *COMPOUND POLYA EGGENBERGER DISTIBUTIONS*

%% CDF/PDF Compound PolyaEggenberger-Pareto1 distribution

alpha = 2;
cfX   = @(t) cfX_Pareto(t,alpha);

a     = 2.2;
b     = 3.3;
m     = 4;
cf    = @(t) cfN_PolyaEggenberger(t,a,b,m,cfX);

clear options
options.isCompound = 1;
CoPoEgPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoPoEgPar1.x CoPoEgPar1.pdf CoPoEgPar1.cdf])
disp(CoPoEgPar1)

%% CDF/PDF Compound PolyaEggenberger-ParetoA distribution

alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2.2;
b     = 3.3;
m     = 4;
cf    = @(t) cfN_PolyaEggenberger(t,a,b,m,cfX);

clear options
options.isCompound = 1;
CoPoEgParA = cf2DistGP(cf,x,prob,options);

%disp([CoPoEgParA.x CoPoEgParA.pdf CoPoEgParA.cdf])
disp(CoPoEgParA)

%% CDF/PDF Compound PolyaEggenberger-ParetoE distribution

alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2.2;
b     = 3.3;
m     = 4;
cf    = @(t) cfN_PolyaEggenberger(t,a,b,m,cfX);

clear options
options.isCompound = 1;
CoPoEgParE = cf2DistGP(cf,x,prob,options);

%disp([CoPoEgParE.x CoPoEgParE.pdf CoPoEgParE.cdf])
disp(CoPoEgParE)

%% CDF/PDF Compound PolyaEggenberger-PearsonVI distribution

alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

a     = 2.2;
b     = 3.3;
m     = 4;
cf    = @(t) cfN_PolyaEggenberger(t,a,b,m,cfX);

clear options
options.isCompound = 1;
CoPoEgPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoPoEgPeaVI.x CoPoEgPeaVI.pdf CoPoEgPeaVI.cdf])
disp(CoPoEgPeaVI)

%% *COMPOUND POISSON DISTIBUTIONS*

%% CDF/PDF Compound Poisson-Pareto1 distribution
alpha = 2;
cfX   = @(t) cfX_Pareto(t,alpha);

a     = 2;
cf    = @(t) cfN_Poisson(t,a,cfX);

clear options
options.isCompound = 1;
CoPoisPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoPoisPar1.x CoPoisPar1.pdf CoPoisPar1.cdf])
disp(CoPoisPar1)

%% CDF/PDF Compound Poisson-ParetoA distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2;
cf    = @(t) cfN_Poisson(t,a,cfX);

clear options
options.isCompound = 1;
CoPoisParA = cf2DistGP(cf,x,prob,options);

%disp([CoPoisParA.x CoPoisParA.pdf CoPoisParA.cdf])
disp(CoPoisParA)

%% CDF/PDF Compound Poisson-ParetoE distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2;
cf    = @(t) cfN_Poisson(t,a,cfX);

clear options
options.isCompound = 1;
CoPoisParE = cf2DistGP(cf,x,prob,options);

%disp([CoPoisParE.x CoPoisParE.pdf CoPoisParE.cdf])
disp(CoPoisParE)

%% CDF/PDF Compound Poisson-PearsonVI distribution 
alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

a     = 2;
cf    = @(t) cfN_Poisson(t,a,cfX);

clear options
options.isCompound = 1;
CoPoisPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoPoisPeaVI.x CoPoisPeaVI.pdf CoPoisPeaVI.cdf])
disp(CoPoisPeaVI)

%% *COMPOUND QUINKERT DISTIBUTIONS*

%% CDF/PDF Compound Quinkert-Pareto1 distribution
alpha = 2;
cfX   = @(t) cfX_Pareto(t,alpha);

a     = 2.2;
b     = 3.4;
cf    = @(t) cfN_Quinkert(t,a,b,cfX);

clear options
options.xMax = 80;
options.isCompound = 1;
CoQuinPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoQuinPar1.x CoQuinPar1.pdf CoQuinPar1.cdf])
disp(CoQuinPar1)

%% CDF/PDF Compound Quinkert-ParetoA distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2.2;
b     = 3.4;
cf    = @(t) cfN_Quinkert(t,a,b,cfX);

clear options
options.xMax = 80;
options.isCompound = 1;
CoQuinParA = cf2DistGP(cf,x,prob,options);

%disp([CoQuinParA.x CoQuinParA.pdf CoQuinParA.cdf])
disp(CoQuinParA)

%% CDF/PDF Compound Quinkert-ParetoE distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a     = 2.2;
b     = 3.4;
cf    = @(t) cfN_Quinkert(t,a,b,cfX);

clear options
options.xMax = 80;
options.isCompound = 1;
CoQuinParE = cf2DistGP(cf,x,prob,options);

%disp([CoQuinParE.x CoQuinParE.pdf CoQuinParE.cdf])
disp(CoQuinParE)

%% CDF/PDF Compound Quinkert-PearsonVI distribution
alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

a     = 2.2;
b     = 3.4;
cf    = @(t) cfN_Quinkert(t,a,b,cfX);

clear options
options.xMax = 80;
options.isCompound = 1;
CoQuinPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoQuinPeaVI.x CoQuinPeaVI.pdf CoQuinPeaVI.cdf])
disp(CoQuinPeaVI)

%% *COMPOUND WARING DISTIBUTIONS*

%% CDF/PDF Compound Waring-Pareto1 distribution
alpha = 2;
cfX   = @(t) cfX_Pareto(t,alpha);

a = 2.2;
b = 3.3;
r = 4;
cf    = @(t) cfN_Waring(t,a,b,r,cfX);

clear options
options.isCompound = 1;
CoWariPar1 = cf2DistGP(cf,x,prob,options);

%disp([CoWariPar1.x CoWariPar1.pdf CoWariPar1.cdf])
disp(CoWariPar1)

%% CDF/PDF Compound Waring-ParetoA distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'amer';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a = 2.2;
b = 3.3;
r = 4;
cf    = @(t) cfN_Waring(t,a,b,r,cfX);

clear options
options.isCompound = 1;
CoWariParA = cf2DistGP(cf,x,prob,options);

%disp([CoWariParA.x CoWariParA.pdf CoWariParA.cdf])
disp(CoWariParA)

%% CDF/PDF Compound Waring-ParetoE distribution
alpha = 1.5;
xmin  = 2.2;
type  = 'eur';
cfX   = @(t) cfX_Pareto(t,alpha,xmin,type);

a = 2.2;
b = 3.3;
r = 4;
cf    = @(t) cfN_Waring(t,a,b,r,cfX);

clear options
options.isCompound = 1;
CoWariParE = cf2DistGP(cf,x,prob,options);

%disp([CoWariParE.x CoWariParE.pdf CoWariParE.cdf])
disp(CoWariParE)

%% CDF/PDF Compound Waring-PearsonVI distribution
alpha = 1.5;
beta  = 2.2;
cfX   = @(t) cfX_PearsonVI(t,alpha,beta);

a = 2.2;
b = 3.3;
r = 4;
cf    = @(t) cfN_Waring(t,a,b,r,cfX);

clear options
options.isCompound = 1;
CoWariPeaVI = cf2DistGP(cf,x,prob,options);

%disp([CoWariPeaVI.x CoWariPeaVI.pdf CoWariPeaVI.cdf])
disp(CoWariPeaVI)

