function rnd = InterpRND(dim,xGiven,cdfGiven)
%InterpRND generates array of size dim of random numbers from the
%   distribution specified by the CDF from the pre-calculated (known) 
%   values of xGiven, and cdfGiven. The evaluation is based on computed
%   quantiles (from the CDF by using the barycentric interpolation) at
%   probabilities generated from the uniform distribution.
%
% SYNTAX:
%   rnd = InterpRND(dim,xGiven,cdfGiven)
%   rnd = InterpRND(dim,result)
%
% INPUT:
%  dim      - size of the array of generated random numbers from the
%             distribution specified by the known CDF, with values cdfGiven
%             that have been pre-computed at xGiven. If dim is scalar, size
%             = [dim,1]. If empty size = [1,1].
%  xGiven   - vector of x-values where the values of CDF are known (have
%             been pre-computed). Alternatively, if this input parameter is
%             missing, it is assumed that the second parameter is the given
%             result (structure) from the inversion algorithm (cf2DistGP).
%  cdfGiven - vector of the known CDF values that have been pre-computed 
%             at xGiven. If this parameter is missing,it is assumed that
%             the second parameter is the given result (structure) from the
%             inversion algorithm (cf2DistGP).
%
% EXAMPLE 1
%  rng(101)
%  cf = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result   = cf2DistGP(cf,[],[],options);
%  xGiven   = result.x;
%  cdfGiven = result.cdf;
%  rnd      = @(dim) InterpRND(dim,xGiven,cdfGiven);
%  dim      = [10000,1];
%  hist(rnd(dim),100);
%
% EXAMPLE 2
%  rng(101)
%  cf = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  rnd    = @(dim) InterpRND(dim,result);
%  dim    = [10000,1];
%  hist(rnd(dim),100);

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Oct-2019 12:31:36

%% ALGORITHM
%  rnd = InterpRND(dim,xGiven,cdfGiven)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);

if nargin < 3
    if isstruct(xGiven)
        result = xGiven;
        xGiven   = result.x;
        cdfGiven = result.cdf;
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
rnd  = InterpQF(prob,xGiven,cdfGiven);
rnd  = reshape(rnd,dim);

end