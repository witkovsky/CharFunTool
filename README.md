# CharFunTool: The Characteristic Functions Toolbox
MATLAB repository of characteristic functions and tools for their combinations and numerical inversion.

For current status of the MATLAB toolbox see the CharFunTool development available at

- https://github.com/witkovsky/CharFunTool

For the R version (not an identical clone) of the toolbox see the CharFun package development available at

- https://github.com/Simkova/CharFun

About
=====

The Characteristic Functions Toolbox (CharFunTool) consists of a set of algorithms for evaluating selected characteristic funcions
and algorithms for numerical inversion of the (combined and/or compound) characteristic functions, used to evaluate the probability density function (PDF) and the cumulative distribution function (CDF).
                                                                              
The toolbox includes several different inversion algorithms, including those based on simple trapezoidal rule for computing the integrals defined by the Gil-Pelaez formulae, and/or by using the FFT algorithm for computing the Fourier transform integrals.
                                                                       
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

We recommend taking a look at the Examples collection. 

To get a taste of what computing with Chebfun is like, type
```
   cf = @(t) exp(-t.^2/2);  % the standard normal characteristic function (CF)
   result = cf2DistGP(cf)   % Invert the CF to get the CDF and PDF   
```


License
=======

See `LICENSE.txt` for CharFunTool licensing information.
