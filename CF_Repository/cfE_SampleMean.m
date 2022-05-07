function [cf,result] = cfE_SampleMean(t,Xdata,Xcounts)
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
%   cf = cfE_DiracMixture(t/n,Xunique,EPMF).^n. 
% 
%   The characteristic function of the sample mean can be used further,
%   e.g. to test hypotheses about the equality of the means (see the
%   Example section).
%
%   SYNTAX
%   cf = cfE_SampleMean(t,Xdata)
%   cf = @(t) cfE_SampleMean(t,Xdata)
%   [~,results] = cfE_SampleMean([],Xdata,Xcounts)
%
%   INPUTS:
%   t       - vector or array of real values, where the CF is evaluated.
%   Xdata   - vector of observed data. If empty, we set the default value
%             Xdata = 0.  
%   Xcounts - vector of the same size as Xdata that represents the counts
%             of the repeated values. This is useful if the data are
%             discrete with many repeted values or if the data are to be
%             specified from a histogram. If Xcounts is a scalar integer
%             value, each value specified by Xdata is repeated equally
%             Xcounts times. If empty, default value is Xcounts = 1. 
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

% (c) 2022 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-May-2022 21:42:52

%% ALGORITHM
%[cf,result] = cfE_SampleMean(t,Xdata,Xcounts)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 4);
if nargin < 2, Xdata = []; end
if nargin < 3, Xcounts = []; end

if isempty(Xdata)
    Xdata = 0;
end

if isempty(Xcounts)
    Xcounts = 1;
end

[errorcode,Xdata,Xcounts] = distchck(2,Xdata(:),Xcounts(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% EMPIRICAL PMF and CDF BASED ON THE OBSERVED DATA
n  = sum(Xcounts);
nc = length(Xcounts);
[Xdata,id] = sort(Xdata);
Xcounts = Xcounts(id);

[Xunique,id] = unique(Xdata);
nX = length(Xunique);
id = [id;nc+1];

Xmin = min(Xunique);
Xmax = max(Xunique);
Xdiff = min(diff(Xunique));

EPMF = zeros(size(Xunique));
ECDF = zeros(size(Xunique));
cdfsum = 0;
for i = 1:nX
    EPMF(i) = sum(Xcounts(id(i):(id(i+1)-1)))/n;
    cdfsum = cdfsum + EPMF(i);
    ECDF(i) = cdfsum;
end
ECDF = max(0,min(1,ECDF));

%% CHARACTERISTI FUNCTION OF THE SAMPLE MEAN
if isempty(t)
    cf = @(t)cfE_DiracMixture(t/n,Xunique,EPMF).^n;
else
    cf = cfE_DiracMixture(t/n,Xunique,EPMF).^n;
end

%% RESULTS
if nargout > 1
    result.description = 'Empirical characteristic function of the sample mean';
    result.cf          = cf;
    result.EPMF        = EPMF;
    result.ECDF        = ECDF;
    result.Xunique     = Xunique;
    result.Xmin        = Xmin;
    result.Xmax        = Xmax;
    result.Xdiff       = Xdiff;
    result.Xdata       = Xdata;
    result.Xcounts     = Xcounts;
    result.n           = n;
end

end