function cdf = InterpCDF(x,xGiven,cdfGiven)
%InterpCDF evaluates the interpolant for the CDF at specified values of x,
%  calculated from the pre-calculated (known) values of xGiven, and
%  cdfGiven. The evaluation is based on the barycentric interpolation.
%
% SYNTAX:
%  cdf = InterpCDF(x,xGiven,cdfGiven)
%  cdf = InterpCDF(x,result)
%
% INPUT:
%  x        - vector of x-values where the CDF is evaluated.
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
%  cf = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result   = cf2DistGP(cf,[],[],options);
%  xGiven   = result.x;
%  cdfGiven = result.cdf;
%  cdf = @(x) InterpCDF(x,xGiven,cdfGiven);
%  x = linspace(-10,10,1001)';
%  plot(x,cdf(x))
%
% EXAMPLE 2
%  cf = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  cdf    = @(x) InterpCDF(x,result);
%  x      = linspace(-10,10,1001)';
%  plot(x,cdf(x))

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Oct-2019 12:31:36

%% ALGORITHM
%  cdf = InterpCDF(x,xGiven,cdfGiven)

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

szx  = size(x);
x = x(:);
id0  = x < min(xGiven);
id1  = x > max(xGiven);
id   = x >= min(xGiven) & x <= max(xGiven);
cdf = zeros(size(x));

if any(id0)
    cdf(id0) = 0;
end

if any(id1)
    cdf(id1) = 1;
end

if any(id)
    cdf(id) = InterpBarycentric(xGiven,cdfGiven,x(id));
end

cdf = max(0,min(1,cdf));
cdf = reshape(cdf,szx);

end