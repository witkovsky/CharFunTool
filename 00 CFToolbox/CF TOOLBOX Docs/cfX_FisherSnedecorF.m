%function cf = cfX_FisherSnedecorF(t,df1,df2,tol) 
% cfX_FisherSnedecorF(t,df1,df2,tol)  evaluates the characteristic function
% cf(t) of the Fisher-Snedecor F-distribution with the parameters df1 and
% df2 (degrees of freedom, df1>=1, df2>=1) computed for real vector
% argument t, i.e. 
%   cf(t) = cfX_FisherSnedecorF(t,df1,df2)  
%         = (gamma(df1/2+df2/2)/gamma(df2/2)) * ...
%           U(df1/2, 1-df2/2, -1i*(df2/df1)*t),
%  where U(a,b,z) denotes the confluent hypergeometric function of the
%  second kind. 
%  For more details see WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/F-distribution 
%
% SYNTAX:
%  cf = cfX_FisherSnedecorF(df1,df2,t,tol)    
%
% EXAMPLE1 (CF of the F-distribution with df1 = 3, df2 = 5)
%  df1 = 3;
%  df2 = 5;
%  t   = linspace(-30,30,2^10+1)';
%  cf  = cfX_FisherSnedecorF(t,df1,df2);
%  plot(t,real(cf),t,imag(cf));grid
%  title('Characteristic function of the F-distribution')
%
% EXAMPLE2 (PDF/CDF of the F-distribution with df1 = 3, df2 = 5)
%  df1 = 3;
%  df2 = 5;
%  x = linspace(0,25,101);
%  prob = [0.9 0.95 0.99];
%  cf = @(t) cfX_FisherSnedecorF(t,df1,df2);
%  clear options
%  options.xMin = 0;
%  options.xMax = 500;
%  options.N  = 2^15;
%  result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Binomial-Fisher-Snedecor distribution)
%  n = 25;  
%  p = 0.3;
%  df1 = 3;
%  df2 = 5;
%  cfX = @(t) cfX_FisherSnedecorF(t,df1,df2);
%  cf = @(t) cfN_Binomial(t,n,p,cfX);
%  x = linspace(0,80,101);
%  prob = [0.9 0.95 0.99];
%  clear options
%  options.isCompound = true;
%  result = cf2DistGP(cf,x,prob,options)
%
% REFERENCES:
% [1] WITKOVSKY, V.: On the exact computation of the density and of
%     the quantiles of linear combinations of t and F random
%     variables. Journal of Statistical Planning and Inference 94
%     (2001), 1–13.
% [2] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  
% [3] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Working Paper.
%     Insurance: Mathematics and Economics. 
% [4] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models.  Working Paper.
%     Journal of Statistical Software.

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfX_FisherSnedecorF(t,df1,df2,tol);