function qf = InterpQF(prob,xGiven,cdfGiven)
%InterpQF evaluates the interpolant for the quantile function (QF) at the
%  specified values of prob, calculated from the pre-calculated (known)
%  values of xGiven, and cdfGiven. The evaluation is based on the
%  barycentric interpolation.
%
% SYNTAX:
%  qf = InterpQF(prob,xGiven,cdfGiven)
%  qf = InterpQF(prob,result)
%
% INPUT:
%  prob     - vector of probabilities from (0,1) where the QF is evaluated. 
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
%  qf       = @(prob) InterpQF(prob,xGiven,cdfGiven);
%  prob     = linspace(0,1,1001)';
%  plot(prob,qf(prob))
%
% EXAMPLE 2
%  cf = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  qf     = @(prob) InterpQF(prob,result);
%  prob   = linspace(eps,1-eps,1001)';
%  plot(prob,qf(prob))

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 04-Oct-2019 12:31:36

%% ALGORITHM
%  qf = InterpQF(prob,xGiven,cdfGiven)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);

if nargin < 3
    if isstruct(xGiven)
        result = xGiven;
        xGiven   = sort(result.x);
        cdfGiven = sort(result.cdf);
    else
        error('Missing Inputs')
    end
end

szx     = size(prob);
prob = prob(:);
id0  = prob < cdfGiven(2);
id1  = prob > cdfGiven(end-1);
id   = prob >= cdfGiven(2) & prob <= cdfGiven(end-1);
qf   = zeros(size(prob));

if any(id0)
    qf(id0) = min(xGiven);
end

if any(id1)
    qf(id1) = max(xGiven);
end

if any(id)
    qf(id) = InterpBarycentric(cdfGiven,xGiven,prob(id));   
end

qf = reshape(qf,szx);

end