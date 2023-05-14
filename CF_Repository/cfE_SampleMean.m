function [cf,result] = cfE_SampleMean(t,data,counts,options)
%% cfE_SampleMean
%   Empirical characteristic function of the sample mean based on the
%   observed random sample X = (X1,...,Xn) of size n.
%
%   The empirical distribution of the sample mean is estimated by the
%   bootstrap distribution of the sample mean (i.e. the bootstrap mean
%   distribution).  
% 
%   The characteristic function of the empirical distribution of the sample
%   mean is derived as the characteristic function of an equally weighted
%   linear combination (arithmetic mean) of n independent random variables
%   each having distribution defined by the observed empirical cumulative
%   distribution function (ECDF), defined from the observed data. This can
%   be expressed as a weighted mixture of the DIRAC distributions whose
%   support is concentrated on the observed unique x-values and the weights
%   are equal to the empirical probability mass function (EPMF) derived
%   from the ECDF.
% 
%   In fact, we define the characteristic function of the sample mean as 
%   cf = cfE_DiracMixture(t/n,X,EPMF).^n, where X is a vector of the
%   observed unique x-values and EPMF is the vector of the weights.
% 
%   The characteristic function of the sample mean can be used further,
%   e.g. to test hypotheses about the equality of the means (see the
%   Example section).
%
%   SYNTAX
%   cf = cfE_SampleMean(t,data)
%   cf = @(t) cfE_SampleMean(t,data)
%   [~,results] = cfE_SampleMean([],data,counts,options)
%
%   INPUTS:
%   t       - vector or array of real values, where the CF is evaluated.
%   data    - vector of observed data. If empty, we set the default value
%             data = 0.  
%   counts  - vector of the same size as data that represents the counts
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
%             sample mean evalueted at the specified vector argument t.
%   result  - structure with further useful details.
% 
%   EXAMPLE: 
%   % Empirical CF of the sample mean
%    data = [1 1 2 2 2 3 3 9 105 105 106 106 106 107 107];
%    [~,result] = cfE_SampleMean([],data);
%    x   = result.Xunique;
%    cf  = result.cf;
%    t   = linspace(-2,2,1001);
%    figure
%    plot(t,real(cf(t)),t,imag(cf(t))) 
%    title('Empirical CF of the Sample Mean')
% 
%   EXAMPLE: 
%   % Empirical CF and PMF/CDF of the sample mean
%    data   = [-2 -1 0 1 2 3 4];
%    counts = [ 2  2 3 1 2 3 2];
%    [~,resultCF] = cfE_SampleMean([],data,counts);
%    x   = resultCF.Xunique;
%    cf  = resultCF.cf;
%    t   = linspace(-100,100,1001);
%    figure
%    plot(t,real(cf(t)),t,imag(cf(t))) 
%    title('Empirical CF of the Sample Mean')
%    n = resultCF.n;
%    delta = 1/n;
%    resultFFT = cf2PMF_FFT(cf,-2,4,delta)
%    % Numerical inversion of the smoothed CF by cf2DistGP
%    cfSmooth = @(t) cf(t) .* cf_Normal(t*delta);
%    resultGP  = cf2DistGP(cfSmooth) 
% 
%   EXAMPLE: 
%   % Empirical distribution of the difference of the sample means
%   % Test of the hypothesis that the population means are equal      
%    dataA   = [1 1 2 2 2 3 3 9 105 105 106 106 106 107 107];
%    dataB   = [5 5 6 6 6 7 7 99 101 101 102 102 102 103 103];
%    cfA     = @(t) cfE_SampleMean(t,dataA);
%    cfB     = @(t) cfE_SampleMean(t,dataB);
%    cfDiff  = @(t) cfA(t) .* cfB(-t);
%    minDiff = min(dataA) - max(dataB);
%    maxDiff = max(dataA) - min(dataB);
%    delta   = 1/15;
%    result = cf2PMF_FFT(cfDiff,minDiff,maxDiff,delta)
% 
%   EXAMPLE: 
%   % Empirical distribution of the difference of the sample means
%   % Test of the hypothesis that the population means are equal      
%    dataA   = [1 1 2 2 2 3 3 9 105 105 106 106 106 107 107];
%    dataB   = [5 5 6 6 6 7 7 99 101 101 102 102 102 103 103];
%    clear options;
%    options.isOgive = true;
%    cfA     = @(t) cfE_SampleMean(t,dataA,[],options);
%    cfB     = @(t) cfE_SampleMean(t,dataB,[],options);
%    cfDiff  = @(t) cfA(t) .* cfB(-t);
%    result = cf2DistGP(cfDiff,[],[0.25,0.975])

% (c) 2022 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-May-2022 17:20:15

%% ALGORITHM
%[cf,result] = cfE_SampleMean(t,data,counts,options)

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

Xcounts = zeros(size(X));
EPMF = zeros(size(X));
ECDF = zeros(size(X));
cdfsum = 0;
for i = 1:nX
    Xcounts(i) = sum(counts(id(i):(id(i+1)-1)));
    EPMF(i) = Xcounts(i)/n;
    cdfsum = cdfsum + EPMF(i);
    ECDF(i) = cdfsum;
end
ECDF = max(0,min(1,ECDF));

%% CHARACTERISTI FUNCTION OF THE SAMPLE MEAN

if options.isOgive
    if isempty(t)
        cf = @(t)cfE_EmpiricalOgive(t/n,X,EPMF).^n;
    else
        cf = cfE_EmpiricalOgive(t/n,X,EPMF).^n;
    end
else
    if isempty(t)
        cf = @(t)cfE_DiracMixture(t/n,X,EPMF).^n;
    else
        cf = cfE_DiracMixture(t/n,X,EPMF).^n;
    end
end

%% RESULTS
if nargout > 1
    result.description = 'Empirical characteristic function of the sample mean';
    result.cf          = cf;
    result.EPMF        = EPMF;
    result.ECDF        = ECDF;
    result.Xunique     = X;
    result.Xmin        = Xmin;
    result.Xmax        = Xmax;
    result.Xdiff       = Xdiff;
    result.data        = data;
    result.counts      = counts;
    result.n           = n;
    result.options     = options;
end

end