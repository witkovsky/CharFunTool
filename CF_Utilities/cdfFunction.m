function cdf = CdfFunction(xNew,xOld,cdfOld)
%CdfFunction evaluates the interpolant for the CDF at new specified values
%  of x, calculated from the pre-calculated (known) values of xOld, and
%  cdfOld. The evaluation is based on the barycentric interpolation. 
%
% SYNTAX:
%  cdf = CdfFunction(xNew,xOld,cdfOld)
%
% INPUT:
%  xNew    - vector of x-values where the CDF is evaluated, 
%  xOld    - vector of x-values where the values of CDF are known (have
%            been pre-computed). Alternatively, if the thisr input
%            parameter is missing, it is assumed that the second parameter
%            is the result structure from the cf2DistGP. 
%  cdfOld  - vector of the known CDF values that have been pre-computed 
%            at xOld. If this parameter is missing, it is assumed that the
%            second parameter is the result structure from the cf2DistGP.
%
% EXAMPLE 1
%  cf     = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  xOld   = result.x;
%  cdfOld = result.cdf;
%  cdf = @(x) CdfFunction(x,xOld,cdfOld);
%  x = linspace(-10,10,1001)';
%  plot(x,cdf(x))
%
% EXAMPLE 2
%  cf     = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  cdf = @(x) CdfFunction(x,result);
%  x = linspace(-10,10,1001)';
%  plot(x,cdf(x))

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Oct-2018 13:46:35

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

szx  = size(xNew);
xNew = xNew(:);
id0  = xNew < min(xOld);
id1  = xNew > max(xOld);
id   = xNew >= min(xOld) & xNew <= max(xOld);
cdf = zeros(size(xNew));

if any(id0)
    cdf(id0) = 0;
end

if any(id1)
    cdf(id1) = 1;
end

if any(id)
    cdf(id) = InterpBarycentric(xOld,cdfOld,xNew(id));
end

cdf     = max(0,min(1,cdf));
cdf     = reshape(cdf,szx);

end