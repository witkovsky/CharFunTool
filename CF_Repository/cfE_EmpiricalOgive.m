function [cf,edges] = cfE_EmpiricalOgive(t,data,counts)
%% cfE_EmpiricalOgive
%   Characteristic function of the OGIVE EMPIRICAL distribution based on
%   the observed histogram given by data and counts (the vector lengths of
%   data and counts are equal) or by edges and counts (the length of the
%   edges is equal to the length of counts + 1).  
% 
%   That is, cf(t) is given as a weighted mixture of the RECTANGULAR
%   (continuous uniform) RVs, 
%    cf(t) = cfE_EmpiricalOgive(t,data,counts) =
%         = sum_{j=1}^n freq_j * cf_Rectangular(t,edges_{j},edges_{j+1}), 
%   where cf_Rectangular(t,edges_{j},edges_{j+1}) is CF of the RECTANGULAR
%   distribution on the interval [edges_{j},edges_{j+1}]. Here, the
%   frequencies are derived from the counts, FREQ = COUNTS/SUM(COUNTS), and
%   the vector EDGES is derived automatically from the DATA (or is
%   specified by the user as an input vector EDGES of size n+1, instead of
%   DATA). In particular, for j = 2,...,n we set EDGES(j) = (DATA(j-1) +
%   DATA(j))/2 and further,  EDGES(1) = DATA(1) - (DATA(2)-DATA(1))/2 and
%   EDGES(n+1) = DATA(n) + (DATA(n)-DATA(n-1))/2.  
% 
%   The EDGES and COUNTS could be specified automatically from DATA by
%   using MATLAB function HISTCOUNTS (histogram counts).
% 
%   In particular, [COUNTS,EDGES] = HISTCOUNTS(DATA) partitions the values
%   in DATA into bins, and returns the count in each bin, as well as the
%   bin edges. HISTCOUNTS uses an automatic binning algorithm that returns
%   bins with a uniform width, chosen to cover the range of elements in
%   DATA and reveal the underlying shape of the distribution. COUNTS(k)
%   will count the value DATA(i) if EDGES(k) <= DATA(i) < EDGES(k+1). The
%   last bin will also include the right edge such that COUNTS(end) will
%   count DATA(i) if EDGES(end-1) <= DATA(i) <= EDGES(end).
%
% SYNTAX
%   cf = cfE_EmpiricalOgive(t,data,counts)
%
% INPUTS:
%  t      - vector or array of real values, where the CF is evaluated.
%  data   - vector of data (if the size is equal to the size of the vector
%           counts) or a vector of edges (if the size is n+1, where n is
%           the length of counts). If empty, default value is data = 1.   
%  counts - vector of counts (integers) related to the values given by
%           data. If empty, default value is counts = 1.  
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
%  https://en.wikipedia.org/wiki/Ogive_(statistics)
%
% EXAMPLE1 (Ogive ECF - a weighted mixture of independent Uniform RVs)
%   rng(101);
%   n = 1000;
%   Data   = [normrnd(5,0.2,3*n,1); trnd(10,n,1)];
%   Edges  = (-6:6)';
%   Counts = histcounts(Data,Edges);
%   t = linspace(-10,10,501)';
%   cf = cfE_EmpiricalOgive(t,Edges,Counts);
%   figure; plot(t,real(cf),t,imag(cf)),grid
%   title('Ogive ECF of the grouped (histogram) data')
%
% EXAMPLE2 (PDF/CDF of the Ogive distribution)
%   rng(101);
%   n = 1000;
%   Data   = [normrnd(5,0.2,3*n,1); trnd(10,n,1)];
%   Edges  = (-6:6)';
%   Counts = histcounts(Data,Edges);
%   cf = @(t) cfE_EmpiricalOgive(t,Edges,Counts);
%   x = linspace(-5,7,101);
%   prob = [0.9 0.95];
%   clear options
%   options.N = 2^12;
%   options.SixSigmaRule = 8;
%   result = cf2DistGP(cf,x,prob,options)
%
% EXAMPLE3 (PDF/CDF of the compound Poisson-Ogive distribution)
%   rng(101);
%   n = 1000;
%   Data = [normrnd(5,0.2,3*n,1); trnd(10,n,1)];
%   Edges  = (-6:6)';
%   Counts = histcounts(Data,Edges);
%   cfX = @(t) cfE_EmpiricalOgive(t,Edges,Counts);
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
%  characteristic function of countsuency and severity. arXiv preprint
%  arXiv:1701.08299.  

% (c) 2022 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-May-2022 15:45:27
% Ver.: 24-Jun-2017 09:34:15

%% ALGORITHM
% cf = cfE_EmpiricalOgive(t,data,counts)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 2, data = []; end
if nargin < 3, counts = []; end

%%
if isempty(data)
    data = 1;
end

if isempty(counts)
    counts = 1;
end

if isscalar(counts)
    [errorcode,data,counts] = distchck(2,data(:)',counts(:)');
    if errorcode > 0
        error(message('InputSizeMismatch'));
    end
end
data = data(:)';
counts = counts(:)';

n = length(counts);
freq = counts/sum(counts);
if length(data) == n
    if length(data) == 1
        edges = zeros(1,2);
        dspan = 1;
        edges(1) = data(1) - dspan;
        edges(2) = data(1) + dspan;
    elseif length(data) == 2
        edges = zeros(1,3);
        dspan = (data(2)-data(1))/2;
        edges(1) = data(1) - dspan;
        edges(2) = data(1) + dspan;
        edges(3) = data(2) + dspan;
    else
        edges = zeros(1,n+1);
        edges(1) = data(1) - (data(2)-data(1))/2;
        edges(2:n) = (data(1:n-1) + data(2:n))/2;
        edges(n+1) = data(n) + (data(n)-data(n-1))/2;
    end
elseif length(data) == n+1
    edges = data;
end

%% Characteristic function 
szt  = size(t);
t    = t(:);
cf   = 0;
aux  = exp(1i * t * edges);
aux  = (aux(:,2:n+1)-aux(:,1:n)) ./ (1i*t.*(edges(2:n+1)-edges(1:n)));
cf   = cf + sum(bsxfun(@times,aux,freq),2);

cf(t==0) = 1;
cf = reshape(cf,szt);

end