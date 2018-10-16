function rnd = RndFunction(dim,xOld,cdfOld)
%RndFunction generates array of size dim of random numbers from the
%   distribution specified by the CDF from the pre-calculated (known) 
%   values of xOld, and cdfOld. The evaluation is based on computed
%   quantiles (from the CDF by using the barycentric interpolation) at
%   probabilities generated from the uniform distribution.
%
% SYNTAX:
%   rnd = RndFunction(n,xOld,cdfOld)
%
% INPUT:
%  dim     - size of the array of generated random numbers from the
%            distribution specified by the known CDF, with values cdfOld
%            that have been pre-computed at xOld. If dim is scalar, size =
%            [dim,1]. If empty size = [1,1].
%  xOld    - vector of x-values where the values of CDF are known (have
%            been pre-computed). Alternatively, if the this input
%            parameter is missing, it is assumed that the second parameter
%            is the result structure from the cf2DistGP.
%  cdfOld  - vector of the known CDF values that have been pre-computed 
%            at xOld. If this parameter is missing, it is assumed that the
%            second parameter is the result structure from the cf2DistGP.
%
% EXAMPLE 1
%  cf      = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result  = cf2DistGP(cf,[],[],options);
%  xOld    = result.x;
%  cdfOld  = result.cdf;
%  rnd     = @(dim) RndFunction(dim,xOld,cdfOld);
%  dim = [10000,1];
%  hist(rnd(dim),51);
%
% EXAMPLE 2
%  cf      = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result  = cf2DistGP(cf,[],[],options);
%  rnd     = @(dim) RndFunction(dim,result);
%  dim = [10000,1];
%  hist(rnd(dim),51);

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Oct-2018 13:46:35

%% ALGORITHM
%  rnd = RndFunction(dim,xOld,cdfOld)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);

if nargin < 3
    if isstruct(xOld)
        result = xOld;
        xOld   = result.x;
        cdfOld = result.cdf;
    else
        error('Missing Inputs')
    end
end

if isempty(dim)
    dim = [1,1];   
end

if isscalar(dim)
    dim = [dim,1];
end

prob = rand(dim(1)*dim(2),1);
rnd  = QfFunction(prob,xOld,cdfOld);
rnd  = reshape(rnd,dim);

end