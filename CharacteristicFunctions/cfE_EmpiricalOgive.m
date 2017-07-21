function cf = cfE_EmpiricalOgive(t,bins,freq)
%% cfE_EmpiricalOgive
%  Characteristic function of the OGIVE EMPIRICAL distribution based on
%  the observed histogram (given by bins and frequencies).
% 
%  That is, cf(t) is given as a weighted mixture of the UNIFORM iid RVs.
%   cf(t) = cfE_EmpiricalOgive(t,bins,freq) =
%         = sum_{j=1}^n freq_j * cf_Uniform(t,bins_{j-1},bins_{j}), 
%  where cf_Uniform(t,bin_{j-1},bin_{j}) is CF of the Uniform distribution
%  on the interval [bins_{j-1},bins_j].
% 
%  The bins and frequencies = counts/sum(counts) are based on the
%  discretized/grouped or histogram data, for more details see e.g. the
%  MATLAB function HISTC.m: COUNTS = HISTC(X,EDGES), counts the number of
%  values in X that fall between the elements in the EDGES vector (which
%  must contain monotonically non-decreasing values). COUNTS is a
%  LENGTH(EDGES) vector containing these counts, COUNTS(k) is count of the
%  X values if, such that EDGES(k) <= X < EDGES(k+1). The last bin counts
%  any values of X that match EDGES(end). 
%
% SYNTAX
%   cf = cfE_EmpiricalOgive(t,bins,freq)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  bins   - vector of data, i.e. constants where the DIRAC RVs are
%           concentrated. If empty, default value is data = 1. 
%  freq   - vector of data, i.e. constants where the DIRAC RVs are
%           concentrated. If empty, default value is data = 1. 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
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
% REFERENCES:
%  WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate
%  loss distribution based on numerical inversion of the compound empirical
%  characteristic function of frequency and severity. arXiv preprint
%  arXiv:1701.08299.  

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jun-2017 09:34:15

%% ALGORITHM
%cf = cfE_EmpiricalOgive(t,bins,freq);

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 2, bins = []; end
if nargin < 3, freq = []; end

%%
if isempty(bins)
    bins = 1;
end

if isempty(freq)
    freq = 1 / length(bins);
end

if isscalar(freq)
    bins = bins(:)';
else
    [errorcode,bins,freq] = distchck(2,bins(:)',freq(:)');
    if errorcode > 0
        error(message('InputSizeMismatch'));
    end
end

%% Characteristic function 
szt   = size(t);
t     = t(:);
cf    = 0;
if length(freq) == length(bins)
    binsUpper  = 2*bins(end)-bins(end-1);
    bins  = [bins,binsUpper];
end
n  = length(bins);
n1 = n-1;

aux  = exp(1i * t * bins);
aux  = (aux(:,2:n)-aux(:,1:n1)) ./ (1i*t.*(bins(:,2:n)-bins(:,1:n1)));
if isscalar(freq)
    cf   = cf + sum(freq*aux,2);
else
    cf   = cf + sum(bsxfun(@times,aux,freq),2);
end

cf(t==0) = 1;
cf = reshape(cf,szt);

end