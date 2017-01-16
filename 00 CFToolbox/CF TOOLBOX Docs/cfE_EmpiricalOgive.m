%function cf = cfE_EmpiricalOgive(t,bins,frequencies,cfX)
%cfE_EmpiricalOgive(t,bins,frequencies,cfX) evaluates the characteristic
%   function cf(t) of the OGIVE empirical distribution function (i.e. the
%   piecewise linear approximation of the empirical CDF based on the
%   observed grouped data) as a weighted mixture of the Uniform iid RVs. In
%   particular,
%     cf(t) = cfE_EmpiricalOgive(t,bins,freq) =
%           = sum_{j=1}^n freq(j) * cf_Uniform(t,bins_{j-1},bins_{j}), 
%   where cf_Uniform(t,bin_{j-1},bin_{j}) is CF of the Uniform distribution
%   on the interval [bins_{j-1},bins_j].
% 
%   The bins and frequencies = counts/sum(counts) are based on the
%   discretized/grouped or histogram data, for more details see e.g. the
%   MATLAB function HISTC.m: N = HISTC(X,EDGES), counts the number of
%   values in X that fall between the elements in the EDGES vector (which
%   must contain monotonically non-decreasing values). N is a LENGTH(EDGES)
%   vector containing these counts, N(k) will count the value X(i) if
%   EDGES(k) <= X(i) < EDGES(k+1).  The last bin will count any values of X
%   that match EDGES(end).
% 
%   cfE_EmpiricalOgive(t,coefs,weights,cf_X) evaluates the compound
%   characteristic function  
%     cf(t) = cfE_EmpiricalOgive(-1i*log(cfX(t)),coefs,weights),
%   where cfX is function handle of the characteristic function cfX(t) of
%   the random variable X (as e.g. another empirical CF based on observed
%   data of X).
%
% SYNTAX
%   cf = cfE_EmpiricalOgive(t,bins,frequencies)
%   cf = cfE_EmpiricalOgive(t,bins,frequencies,cfX)
%
% EXAMPLE1 (Ogive ECF - a weighted mixture of independent Uniform RVs)
%   rng(101);
%   n = 1000;
%   Data = [normrnd(5,0.2,3*n,1); trnd(10,n,1)];
%   Edges = (-6:6)';
%   Counts = histc(Data,Edges);
%   Xbins = Edges;
%   Xfreq = Counts/sum(Counts);
%   t = linspace(-10,10,501)';
%   cf = cfE_EmpiricalOgive(t,Xbins,Xfreq);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('Ogive ECF of the grouped (histogram) data')
%
% EXAMPLE2 (PDF/CDF of the Ogive distribution)
%   rng(101);
%   n = 1000;
%   Data = [normrnd(5,0.2,3*n,1); trnd(10,n,1)];
%   Edges = (-6:6)';
%   Counts = histc(Data,Edges);
%   Xbins = Edges;
%   Xfreq = Counts/sum(Counts);
%   cf = @(t) cfE_EmpiricalOgive(t,Xbins,Xfreq);
%   x = linspace(-5,7,101);
%   prob = [0.9 0.95];
%   clear options
%   options.N = 2^12;
%   options.SixSigmaRule = 8;
%   result = cf2DistFFT(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Ogive distribution)
%   rng(101);
%   n = 1000;
%   Data = [normrnd(5,0.2,3*n,1); trnd(10,n,1)];
%   Edges = (-6:6)';
%   Counts = histc(Data,Edges);
%   Xbins = Edges;
%   Xfreq = Counts/sum(Counts);
%   cfX = @(t) cfE_EmpiricalOgive(t,Xbins,Xfreq);
%   lambda = 15;
%   cf = @(t) cfN_Poisson(t,lambda,cfX);
%   x = linspace(0,100,101);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.isCompound = true;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE4 (PDF/CDF of the compound DanishFireData distribution)
%   Data = [normrnd(5,0.2,3*n,1); trnd(10,n,1)];
%   Xbins = dlmread('GroupBins.txt');
%   Xfreq = dlmread('GroupFrequencies.txt');
%   cfX = @(t) cfE_EmpiricalOgive(t,Xbins,Xfreq);
%   cf = @(t) cfN_Poisson(t,lambda,cfX);
%   x = linspace(0,100,101);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.isCompound = true;
%   result = cf2DistGP(cf,x,prob,options)
% 
% REFERENCES:
% [1] WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
%     distribution based on numerical inversion of the compound empirical
%     characteristic function of frequency and severity. Preprint submitted
%     to Insurance: Mathematics and Economics.
% [2] DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
%     computing distributions of collective risk models. Preprint submitted
%     to Journal of Statistical Software.
% [3] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44. 
% [4] WIMMER G., ALTMANN G. (1999). Thesaurus of univariate discrete
%     probability distributions. STAMM Verlag GmbH, Essen, Germany. ISBN
%     3-87773-025-6. 

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfE_EmpiricalOgive(t,bins,frequencies,cfX);
