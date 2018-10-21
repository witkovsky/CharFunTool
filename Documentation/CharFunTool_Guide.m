%% CharFunTool User Guide
%% The Characteristic Fuctions Toolbox for MATLAB
%% 
% 
% 
% *1st Edition (21-Oct-2018)*
% 
% *CharFunTool version 1.4 (21-Oct-2018)*
% 
% *Edited by: Viktor Witkovsky *
% 
% **
% 
% *CharfunTool was created in January 2017 by: *
% 
% *Viktor Witkovsky*
% 
% *and is continuously developing with help of the (still growing) CharFunTool 
% Development Team*
% 
% **
% 
% *Current members: *
% 
% *Viktor Witkovsky, Gejza Wimmer, Tomas Duby, Andrej Gajdoš, Jozef Hanè*
% 
% Former members: ¼udmila Šimková
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% Copyright 2018 by Viktor Witkovsky and the Institute of Measurement Science 
% of the Slovak Academy of Sciences, Bratislava, Slovakia. All rights reserved.
% 
% For more information, write to witkovsky@gmail.com.
% 
% MATLAB is a registered trademark of The MathWorks, Inc. For MATLAB product 
% information, please write to info@mathworks.com.
% 
% 
%% Preface
% 
% 
% This guide is an introduction to the use of CharFunTool, an open source 
% software package that aims to provide repository of algorithms for computing 
% the characteristic functions and tools for their combinations and numerical 
% inversion.
% 
% Application of the exact statistical inference frequently leads to a non-standard 
% probability distributions of the considered estimators or test statistics. Frequently, 
% evaluation of the probability density function (PDF), cumulative distribution 
% function (CDF), and/or the quantile function (QF) is possible from the characteristic 
% function (CF). 
% 
% In many important situations, derivation of the CFs is more simple than 
% derivation of the PDFs and/or CDFs. In particular, the exact distribution of 
% many estimators and test statistics can be structurally expressed as a linear 
% combination or product of independent random variables with known distributions 
% and characteristic functions as is the case for many standard multivariate test 
% criteria. However, analytical inversion of the characteristic function (if possible) 
% frequently leads to complicated expressions of the corresponding distribution 
% functions, CDF/PDF and the required quantiles. 
% 
% As we shall illustrate here, for many applications, the method based on 
% simple implementation of the numerical inversion of the characteristic functions 
% is fully sufficient.
% 
% I gratefully acknowledge the support of the Institute of Measurement Science 
% of the Slovak Academy of Sciences in Bratislava, Slovakia, and also the support 
% by the Slovak Research and Development Agency, and by the Scientific Grant Agency 
% of the Ministry of Education of the Slovak Republic and the Slovak Academy of 
% Sciences, projects APVV-15-0295, VEGA 2/0054/18, VEGA 2/0011/16. 
% 
% Most especially I acknowledge and praise my collaborators and especially 
% the contributors to CharFunTool. For their efficient cooperation and development 
% efforts, I thank especially to  Gejza Wimmer (Mathematical Institute  of the 
% Slovak Academy of Science in Bratislava), Tomas Duby (OAA Computing, Oxfordshire, 
% UK), Andrej Gajdoš and Jozef Hanè (Pavol Jozef Šafárik University in Košice), 
% and ¼udmila Šimková (former student at the Comenius University in Bratislava). 
% Further, for their support and discussions, I would like to thank my colleagues 
% at the Institute of Measurement Science, František Rublík, Anna Krakovská, Martina 
% Chvosteková, and Jozef Jakubík, as well as the project partners and collaborators, 
% Stanislav Ïuriš, Rudolf Palenèár (Slovak University of Technology in Bratislava), 
% Zuzana Ïurišová (Slovak Institute of Metrology).
% 
% Viktor Witkovsky
% 
% Bratislava, October 2018
%% *CharFunTool:The Characteristic Functions Toolbox*
% 
% 
% The Characteristic Functions Toolbox (CharFunTool) is a MATLAB repository 
% of characteristic functions and tools for their combinations and numerical inversion.
% 
% CharFunTool consists of a set of algorithms for evaluating selected characteristic 
% functions and algorithms for numerical inversion of the combined and/or compound 
% characteristic functions, used to evaluate the cumulative distribution function 
% (CDF), the probability density function (PDF), and/or the quantile function 
% (QF).
% 
% 
%% Instalation
% 
% 
% To install, you can either clone the directory with Git or download a .zip 
% file. 
% 
% * Option 1: Download .zip file
% 
% Download a .zip of CharFunTool from
%% Inversion algorithms
% The toolbox comprises different inversion algorithms, including those based 
% on simple trapezoidal quadrature rule for computing the integrals defined by 
% the Gil-Pelaez formulae, and/or based on using the FFT algorithm for computing 
% the Fourier transform integrals.
%% Algorithms for numerical inversion of the characteristic functions
%% cf2DistBV
%% cf2DistFFT
%% cf2DistGP
%% Repository of characteristic funtions
% CharFunTool consists of a set of algorithms for evaluating selected characteristic 
% functions and algorithms for numerical inversion of the combined and/or compound 
% characteristic functions,
%% Characteristic functions of general probability distributions
%% cf_ArcsineSymmetric
% cf_ArcsineSymmetric evaluates the characteristic function of a linear combination 
% (resp. convolution) of independent zero-mean symmetric ARCSINE random variables 
% defined on the interval $(-1,1)$.
%% cf_Beta
%% cf_BetaNC
%% cf_BetaSymmetric
%% cf_Chi
%% cf_ChiNC
%% cf_ChiSquare
%% cf_Exponential
%% cf_FisherSnedecor
%% cf_FisherSnedecorNC
%% cf_Gamma
%% cf_GeneralizedExponential
%% cf_Gumbel
%% cf_InverseGamma
%% cf_Laplace
%% cf_MaxwellBoltzmann
%% cf_MaxwellBoltzmannNC
%% cf_Nakagami
%% cf_NakagamiNC
%% cf_Normal
%% cf_Rayleigh
%% cf_RayleighNC
%% cf_Rectangular
%% cf_RectangularSymmetric
%% cf_Rice
%% cf_Stable
%% cf_Student
%% cf_TrapezoidalSymmetric
%% cf_TriangularSymmetric
%% cf_TSPSymmetric
%% cf_Weibull
%% cf_WignerSemicircle
%% Tools for manuipulating with the characteristic function
% CharFunTool consists of a set of algorithms for manipulating and combining 
% the characteristic functions.
%% Utility functions