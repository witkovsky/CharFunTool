function qf = QfFunction(prob,xOld,cdfOld)
%QfFunction evaluates the interpolant for the quantile function (QF) at new
%  specified values of prob, calculated from the pre-calculated (known)
%  values of xOld, and cdfOld. The evaluation is based on the barycentric
%  interpolation.  
%
% SYNTAX:
%  qf = QfFunction(prob,xOld,cdfOld)
%
% INPUT:
%  prob    - vector of probabilities from (0,1) where the QF is evaluated, 
%  xOld    - vector of x-values where the values of CDF are known (have
%            been pre-computed). Alternatively, if the this input
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
%  cdfOld  = result.cdf;
%  qf     = @(prob) QfFunction(prob,xOld,cdfOld);
%  prob   = linspace(0,1,1001)';
%  plot(prob,qf(prob))
%
% EXAMPLE 2
%  cf     = @(t) exp(-t.^2/2);
%  clear options
%  options.isPlot = false;
%  options.isInterp = true;
%  result = cf2DistGP(cf,[],[],options);
%  qf     = @(prob) QfFunction(prob,result);
%  prob   = linspace(0,1,1001)';
%  plot(prob,qf(prob))

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 16-Oct-2018 13:46:35

%% ALGORITHM
%  qf = QfFunction(prob,xOld,cdfOld)

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

szx     = size(prob);
prob = prob(:);
id0  = prob < min(cdfOld);
id1  = prob > max(cdfOld);
id   = prob >= min(cdfOld) & prob <= max(cdfOld);
qf   = zeros(size(prob));

if any(id0)
    qf(id0) = min(xOld);
end

if any(id1)
    qf(id1) = max(xOld);
end

if any(id)
    qf(id) = InterpBarycentric(cdfOld,xOld,prob(id));   
end

qf     = reshape(qf,szx);

end