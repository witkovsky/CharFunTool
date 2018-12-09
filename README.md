# CharFunTool: The Characteristic Functions Toolbox
MATLAB repository of characteristic functions and tools for their combinations and numerical inversion.

For current status of the MATLAB toolbox see the CharFunTool development available at

- https://github.com/witkovsky/CharFunTool

For the R version (not an identical clone) of the toolbox see the CharFun package development available at

- https://github.com/gajdosandrej/CharFunToolR

About
=====

The Characteristic Functions Toolbox (CharFunTool) consists of a set of algorithms for evaluating selected characteristic functions and algorithms for numerical inversion of the combined and/or compound characteristic functions, used to evaluate the cumulative distribution function (CDF), the probability density function (PDF), and/or the quantile function (QF).
                                                                              
The toolbox comprises different inversion algorithms, including those based on the Gil-Pelaez inversion formulae in combination with the simple trapezoidal quadrature rule, or other more sofisticated quadratures and advanced acceleration methods, used for computing the  required Fourier transform integrals of oscillatory functions.
                                                                       
Installation and requirements
=============================

CharFunTool was developed with MATLAB Version: 9.2 (R2017a).

To install, you can either clone the directory with Git or download a .zip file. 

## Option 1: Download .zip file

Download a .zip of CharFunTool from

- https://github.com/witkovsky/CharFunTool/archive/master.zip

After unzipping, you will need to add CharFunTool to the MATLAB path. You can do this either (a) by typing
```
addpath(CharFunToolRoot), savepath
```
where `CharFunToolRoot` is the path to the unzipped directory, (b) by selecting the `CharFunTool` directory with the `pathtool` command, or (c) though the File > Set Path... dialog from the MATLAB menubar.

## Option 2: Clone with Git

To clone the CharFunTool repository, first navigate in a terminal to where you want the repository cloned, then type
```
git clone https://github.com/witkovsky/CharFunTool.git
```
To use CharFunTool in MATLAB, you will need to add the `CharFunTool` directory to the MATLAB path as above.


Getting started
===============

We recommend taking a look at the Examples collection and the detailed helps of the included characteristic functions and the inversion algorithms. 

* To get a taste of what computing with CharFunTool is like, try to invert the characteristic function (CF) of a standard Gaussian distribution. For that, simply type
```
 cf = @(t) exp(-t.^2/2);                  % the standard normal (Gaussian) distribution characteristic function
 result = cf2DistGP(cf)                   % Invert the CF to get the CDF and PDF   
```
* For a more advanced distribution, based on the theory of Gaussian Processes, type 
```
 df     = 1;
 cfChi2 = @(t) (1-2i*t).^(-df/2);         % CF of the chi-squared distribution

 idx    = 1:100;
 coef   = 1./((idx-0.5)*pi).^2;
 cf     = @(t) cf_Conv(t,cfChi2,coef);    % CF of the linear combination of iid RVs 

 prob   = [0.9 0.95 0.99];
 x      = linspace(0,3,500);

 clear options
 options.N = 2^12;
 options.xMin = 0;
 options.SixSigmaRule = 10;

 result = cf2DistGP(cf,x,prob,options)    % Invert the CF to get the CDF and PDF 
```
* Alternatively, by using the included characteristic function `cfTest_EqualityCovariances` and the inversion algorithm `cf2DistGPT` (based on the Gil-Pelaez inversion formula and the simple trapezoidal quadrature rule) evaluate the exact null-distribution (PDF/CDF and the selected quantiles) of the negative log-transformed Likelihood Ratio Test (LRT) statistic for testing the hypothesis on equality of covariance matrices in q normal p-dimensional populations, based on random samples of size n (n > p for each population). 

* In particular, for testing the null hypothesis H0: Sigma_1 = ... = Sigma_q we consider the log-transformed LRT statistic W = -log(LRT), where LRT = [ q^(pq) prod_{k=1}^q det(S_k) / det(S)^q ]^(n/2) with S_k being the maximum likelihood estimators (MLEs) of the matrices Sigma_k and S = S_1 + ... + S_q. Then the null distribution of W = -log(LRT) can be estimated by 
```
  q    = 3;                               % number of independent normal populations (q >= 2)
  p    = 5;                               % dimension of each normal population (p >= 1)
  n    = 10;                              % sample size taken from each normal population (n > p) 
                                          % i.e. X_{k,j} ~ N_p(mu_k,Sigma_k), k = 1,...,q, j = 1,...,n
  type = 'standard';                      % type of considered LRT statistic ('standard' or 'modified')
  
  % CF of the null-distribution of the negative log-transformed LRT statistic
  % for testing the null hypothesis H0: Sigma_1 = ... = Sigma_q
  cf   = @(t) cfTest_EqualityCovariances(t,n,p,q,type);   
  
  x    = linspace(0,50,201);              % x values where the PDF/CDF is to be evaluated
  prob = [0.9 0.95 0.99];                 % probabilities for which the quantiles are calculated
  clear options                           % clear/set the options structure
  options.xMin = 0;                       % set the known minimal value (e.g. for non-negative distribution) 
  
  result = cf2DistGPT(cf,x,prob,options)  % Invert the CF to get the CDF/PDF/QF and other results
```
License
=======

See `LICENSE.txt` for CharFunTool licensing information.
