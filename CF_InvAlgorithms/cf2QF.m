function [qf,pr,cdf,x,xMin,xMax,xMean,xStd,cdfMin,cdfMax] =  ...
    cf2QF(cf,prob,N,chebyPts,SixSigmaRule,tolDiff,interpMethod)
% cf2QF Fast GPR algorithm (Gil-Pelaez inversion formulae with Riemann sum
%  used for the integration) to evaluate the quantiles of the continuous 
%  univariate distribution specified by its characteristic function CF
%  based on using barycentric aor linear interpolation of the fitted CDF.
%  
% SYNTAX:
%   [qf,prob,cdf,x,xMin,xMax,xMean,xStd,cdfL,cdfU] = ...
%       cf2QF(cf,prob,N,chebyPts,SixSig,tolDiff,interpMethod)
% INPUT:
%  cf       - function handle of the characteristic function (CF).
%  prob     - vector of values from [0,1] for which the quantiles
%             function is evaluated.
%  N        - number of points used for trapezoidal integration. If empty,
%             by default N = 2^10. 
%  chebyPts - number of Chebyshev points to create the vector of x values,
%             x = [xMin, xMax], where the CDF of the distribution specified
%             by its CF is evaluated (xMin and xMax are specified
%             automatically from the CF). If empty, by default chebyPts =
%             2^7.  
%  SixSigmaRule - multiple of standard deviation which specifies the rule
%             for estimating the principal domain of the distribution. If
%             empty, by default SixSigmaRule = 6. 
%  tolDiff  - parameter used for numerical differentiation. If empty, by
%             default tolDiff = 1e-4.  
%  interpMethod - specifies the interpolation method used for evaluating
%             the quantiles from the computed CDF values at prespecified x
%             values. If empty, by default interpMethod = 'barycentric'
%             (also interpMethod = 1). Alternatively, for faster evaluation
%             use interpMethod = 'linear' (also interpMethod = 0).
% OUTPUT:
%  qf       - vector of quantiles evaluated for the probabilities
%             prespecified in prob.  
%  pr       - vector of values from [0,1] for which the quantiles
%             were evaluated. This is corrected version of prob, such that
%             pr(prob<cdfMin) = cdfMin and pr(prob>cdfMax) = cdfMax. 
%  cdf      - CDF evaluated for the automatically created vector x, used
%             further for evaluating the quantiles by interpolating its
%             values for the prespecified probabilities in prob. 
%  x        - automatically created vector of values, Chebyshev points from
%             [xMin,xMax] of size chebyPts for which the quantiles were
%             evaluated. 
%  xMin     - automatically calculated minimal value of the principal
%             domain of the distribution specified by its CF (with respect
%             to the used SixSigmaRule).
%  xMax     - automatically calculated maximal value of the principal
%             domain of the distribution specified by its CF (with respect
%             to the used SixSigmaRule).
%  xMean    - automatically calculated estimate of the mean value of the
%             the distribution specified by its CF.
%  xStd     - automatically calculated estimate of the standard deviation
%             of the the distribution specified by its CF.
%  cdfMin   - calculated minimal value of the CDF on the used principal
%             domain of the distribution, specified by [xMin, xMax],
%             derived from its CF (with respect to the used SixSigmaRule).
%  cdfMax   - calculated maximal value of the CDF on the used principal
%             domain of the distribution, specified by [xMin, xMax],
%             derived from its CF (with respect to the used SixSigmaRule).
%
% REMARK:
%  cf2QF was suggested as a fast version of the other inversion algorithm,
%  in particular cf2DistGPR, in order to calculate the approximately the
%  required quantiles. This can be used for random sampling from the
%  distribution specified its characteristic function CF.
%
% EXAMPLE (Calculate specified quantiles from distribution specified by CF)
%  cf = @(t) exp(-t.^2/2);
%  prob = rand(100,2); 
%  [qf,prob] = cf2QF(cf,prob)
%
% EXAMPLE (Generate Random Sample from distribution specified by CF)
%  cf = @(t) cf_Logistic(t);
%  n = 10000;
%  p = rand(n,1);
%  r = cf2QF(cf,p);
%  histogram(r,100,'normalization','pdf')
%
% EXAMPLE (ECDF from distribution specified by CF)
%  cf = @(t) cf_Logistic(t);
%  n = 10000;
%  p = rand(n,1);
%  [r,pr,cdf,x] = cf2QF(cf,p);
%  ecdf(r), hold on
%  plot(x,cdf), hold off

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 05-Jan-2022 17:34:23
%% CHECK THE INPUT PARAMETERS
narginchk(1, 7);
if nargin < 7, interpMethod = []; end
if nargin < 6, tolDiff = []; end
if nargin < 5, SixSigmaRule = []; end
if nargin < 4, chebyPts = []; end
if nargin < 3, N = []; end
if nargin < 2, prob = []; end

if isempty(interpMethod), interpMethod  = 'barycentric'; end
if isempty(tolDiff), tolDiff  = 1e-4; end
if isempty(SixSigmaRule), SixSigmaRule  = 6; end
if isempty(chebyPts), chebyPts = 2^7; end
if isempty(N), N = 2^10; end
if isempty(prob), prob = linspace(0,1,chebyPts)'; end

%% ALGORITHM
szp      = size(prob);
prob     = prob(:);
cft      = cf(tolDiff*(1:4));
cftRe    = real(cft);
cftIm    = imag(cft);
xMean    = (8*cftIm(1)/5 - 2*cftIm(2)/5 + 8*cftIm(3)/105 - 2*cftIm(4)/280) / tolDiff;
xM2      = (205/72 - 16*cftRe(1)/5 + 2*cftRe(2)/5 - 16*cftRe(3)/315 + 2*cftRe(4)/560) / tolDiff^2;
xStd     = sqrt(xM2 - xMean^2);
xMin     = xMean - SixSigmaRule * xStd;
xMax     = xMean + SixSigmaRule * xStd;
range    = xMax - xMin;
dt       = 2*pi / range;
t        = (0.5+(0:N-1))' * dt;
cft      = cf(t);
cft(N)   = cft(N)/2;
x        = range * (-cos(pi*(0:(chebyPts-1)) / (chebyPts-1)) + 1)/2 + xMin;
x        = x(:);
E        = exp(-1i*x*t');
cdf      = imag(E * (cft ./ t));
cdf      = 0.5 - (cdf * dt) / pi;
cdfMin   = min(cdf);
cdfMax   = max(cdf);
prob(prob<cdfMin) = cdfMin;
prob(prob>cdfMax) = cdfMax;
[ps,pID] = sort(prob);

switch lower(interpMethod)
    case {'linear',0}
        qf(pID) = interp1(cdf,x,ps);
    case {'barycentric',1}
        qf(pID) = InterpBarycentric(cdf,x,ps);
    otherwise
        qf(pID) = interp1(cdf,x,ps);
end
qf = reshape(qf,szp);
pr = reshape(prob,szp); 

end