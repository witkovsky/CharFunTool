function [cf,result] = cfE_SampleMedian(t,data,counts,options)
%% cfE_SampleMedian
%   Empirical characteristic function of the sample median based on the
%   observed data - the realization of the random sample X = (X1,...,Xn) of
%   size n. 
%
%   The distribution of the sample median is estimated by its empirical
%   distribution, which is defined as the distribution of the kth-order
%   statistics derived from the observed empirical cumulative distribution
%   function (ECDF). The order of the median is given by k = (n+1)/2.
% 
%   Note that for the continuous distribution F, the kth-order statistic
%   has a beta distribution and its CDF values at x are given as CDF(x) =
%   betainc(F(x),k,n+1-k). Surprisingly, this is also true for discrete
%   distributions, such as the ECDF, when evaluated at the support points,
%   i.e. CDF (x) = betainc(ECDF(x),k,n+1-k). The empirical CDF of the
%   sample median is then CDF (x) = betainc(ECDF(x),(n+1)/2,(n+1)/2). 
%   
%   When the sample size n is an odd integer, CDF corresponds to the
%   bootstrap distribution of the sample median (i.e. the limiting
%   bootstrap distribution calculated from the large number of bootstrap
%   samples B -> infinity).
% 
%   When the sample size n is an even integer, the sample median is usually
%   defined as the mean of the kth and the (k+1)th order statistics, where
%   k = n/2. However, CDF (x) = betainc(ECDF(x),(n+1)/2,(n+1)/2) is also
%   well defined for n an even integer (but this CDF would be slightly
%   different from the bootstrap distribution as it could be concentrated
%   also on points that are out of the observed data). Therefore, we define
%   here the empirical distribution of the sample median for both cases,
%   even and odd n, as CDF (x) = betainc(ECDF(x),(n+1)/2,(n+1)/2).
% 
%   The characteristic function of the empirical distribution of the sample
%   median is derived as the characteristic function of a weighted discrete
%   mixture distribution whose support is concentrated on the observed
%   unique x-values and whose weights are equal to the empirical
%   probability mass function (PMF) derived from CDF (x) =
%   betainc(ECDF(x),(n+1)/2,(n+1)/2). In fact, we define the characteristic
%   function of the sample median as cf = cfE_DiracMixture(t,Xunique,PMF).
% 
%   The characteristic function of the sample median can be used further,
%   e.g. to test hypotheses about the equality of the medians (see the
%   Example section).
%
%   SYNTAX
%   cf = cfE_SampleMedian(t,data)
%   cf = @(t) cfE_SampleMedian(t,data)
%   [~,results] = cfE_SampleMedian([],data,counts,options)
%
%   INPUTS:
%   t       - vector or array of real values, where the CF is evaluated.
%   data   - vector of observed data. If empty, we set the default value
%             data = 0.  
%   counts - vector of the same size as data that represents the counts
%             of the repeated values. This is useful if the data are
%             discrete with many repeted values or if the data are to be
%             specified from a histogram. If counts is a scalar integer
%             value, each value specified by data is repeated equally
%             counts times. If empty, default value is counts = 1. 
%   options - structure with further optional specificaions:
%             - options.isOgive = false.
%
%   OUTPUTS:
%   cf      - characteristic function of the empirical distribution of the
%             sample median evalueted at the specified vector argument t.
%   result  - structure with further useful details.
% 
%   EXAMPLE: 
%   % Empirical CF and PMF/CDF of the sample median
%    dataA   = [1 1 2 2 2 3 3 9 105 105 106 106 106 107 107];
%    [~,result] = cfE_SampleMedian([],dataA);
%    x   = result.Xunique;
%    pmf = result.PMFmedian;
%    cdf = result.CDFmedian;
%    cf  = result.cf;
%    t   = linspace(-10,10,1001);
%    figure
%    stem(x,pmf)
%    title('Empirical PMF of the Sample Median')
%    figure
%    stairs(x,cdf)
%    title('Empirical CDF of the Sample Median')
%    figure
%    plot(t,real(cf(t)),t,imag(cf(t))) 
%    title('Empirical CF of the Sample Median')
% 
%   EXAMPLE: 
%   % Empirical distribution of the difference of the sample medians
%   % Test of the hypothesis that the population medians are equal      
%    dataA   = [1 1 2 2 2 3 3 9 105 105 106 106 106 107 107];
%    dataB   = [5 5 6 6 6 7 7 99 101 101 102 102 102 103 103];
%    cfA     = @(t) cfE_SampleMedian(t,dataA);
%    cfB     = @(t) cfE_SampleMedian(t,dataB);
%    cfDiff  = @(t) cfA(t) .* cfB(-t);
%    minDiff = min(dataA) - max(dataB);
%    maxDiff = max(dataA) - min(dataB);
%    delta   = 1;
%    result = cf2PMF_FFT(cfDiff,minDiff,maxDiff,delta)
%
%  REFERENCES
%  [1] Divine, G. W., Norton, H. J., Barón, A. E., & Juarez-Colunga, E.
%      (2018). The Wilcoxon–Mann–Whitney procedure fails as a test of
%      medians. The American Statistician, 72(3), 278-286.   

% (c) 2022 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-May-2022 17:20:15

%% ALGORITHM
%[cf,result] = cfE_SampleMedian(t,data,counts,options)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, data = []; end
if nargin < 3, counts = []; end
if nargin < 4, options = []; end

if isempty(data)
    data = 0;
end

if isempty(counts)
    counts = 1;
end

if ~isfield(options, 'isOgive')
    options.isOgive = false;
end

[errorcode,data,counts] = distchck(2,data(:),counts(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% EMPIRICAL PMF and CDF BASED ON THE OBSERVED DATA
n  = sum(counts);
nc = length(counts);
[data,id] = sort(data);
counts = counts(id);

[X,id] = unique(data);
nX = length(X);
id = [id;nc+1];

Xmin = min(X);
Xmax = max(X);
Xdiff = min(diff(X));

EPMF = zeros(size(X));
ECDF = zeros(size(X));
cdfsum = 0;
for i = 1:nX
    EPMF(i) = sum(counts(id(i):(id(i+1)-1)))/n;
    cdfsum = cdfsum + EPMF(i);
    ECDF(i) = cdfsum;
end
ECDF = max(0,min(1,ECDF));

%% EMPIRICAL PMF and CDF of the SAMPLE MEDIAN
medianOrder = (n+1)/2;
CDFmedian   = betainc(ECDF,medianOrder,medianOrder);
PMFmedian   = diff([0;CDFmedian]);

%% CHARACTERISTI FUNCTION OF THE SAMPLE MEDIAN
if options.isOgive
    if isempty(t)
        cf = @(t)cfE_EmpiricalOgive(t,X,PMFmedian);
    else
        cf = cfE_EmpiricalOgive(t,X,PMFmedian);
    end
else
    if isempty(t)
        cf = @(t)cfE_DiracMixture(t,X,PMFmedian);
    else
        cf = cfE_DiracMixture(t,X,PMFmedian);
    end
end

%% RESULTS
if nargout > 1
    result.description = 'Empirical characteristic function of the sample median';
    result.cf          = cf;
    result.PMFmedian   = PMFmedian;
    result.CDFmedian   = CDFmedian;
    result.EPMF        = EPMF;
    result.ECDF        = ECDF;
    result.Xunique     = X;
    result.Xmin        = Xmin;
    result.Xmax        = Xmax;
    result.Xdiff       = Xdiff;
    result.data       = data;
    result.counts     = counts;
    result.n           = n;
    result.options     = options;
end

end