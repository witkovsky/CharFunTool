function [pdf,result] = cf2PDF(cf,x,options)
%cf2PDF  Evaluates the PDF from the characteristic function CF
% by the Gil-Pelaz inversion formula. The algorithm uses the adaptive
% Gauss-Kronod quadrature rule for numerical integration of the oscillatory
% integrand function devided into sub-intervals (found by a fast
% root-finding algorithm) and subsequent application of the convergence
% acceleration techniques for computing the limit of the resulted
% alternating series.
%
% SYNTAX:
%  [pdf,result] = cf2PDF(cf,x,options)
%
% INPUT:
%  cf       - function handle of the characteristic function
%  x        - vector of x values where the PDF is computed
%  options  - structure with the following parameters:
%  nPeriods   - the upper integration limit: UPPER = nPeriods * pi / x.
%               The the basic integration interval [0,UPPER] is devided
%               into two subintervals [0 A] and [A UPPER].
%  isPlot     - logical indicator for plotting the integrand function and
%               calculation of their zeros, default value is options.isPlot
%               = false.
% OUTPUT:
%  pdf      - vector of PDF values evaluated at x.
%  result   - structure with PDF  further details:
%
% EXAMPLE:
% % PDF of the standard Normal distribution
%  cf  = @(t) exp(-t.^2/2);
%  x   = [1.281551565544601 1.644853626951472 2.326347874040841];
%  pdf = cf2PDF(cf,x)
%
% EXAMPLE:
% % PDF of the Chi-squared distribution with DF = 1
%  df = 1;
%  cf = @(t) (1 - 2i*t).^(-df/2);
%  x  = [2.705543454095416 3.841458820694126 6.634896601021214];
%  [pdf,result] = cf2PDF(cf,x)
%
% EXAMPLE:
% % PDF of the linear combination (convolution) of central Chi-squared RVs
%  df   = [1 2 3];
%  cf = @(t) cf_ChiSquare(t,df) ;
%  %x = [10.644640675668420  12.591587243743977  16.811893829770927];
%  x = 25;
%  clear options
%  options.isPlot = true;
%  [pdf,result] = cf2PDF(cf,x,options)
%
% REFERENCES
%  Imhof, J.: Computing the distribution of quadratic forms in normal
%  variables. Biometrika 48, 419–426 (1961).
%
% REMARKS:
%  This version of cf2PDF was optimized for computing PDF of the
%  linear combination of independent chi-squared random variables.
%  cf2PDF is a general purpose algorithm, suitable for numerical
%  inversion of other well defined characteristic functions (with carefully
%  tuned option parameters).

% (c) 2017 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 08-Apr-2017 18:10:40

%% CHECK/SET THE INPUT PARAMETERS
tic;
narginchk(2, 3);

if nargin < 3, options = []; end

if ~isfield(options, 'nPeriods')
    options.nPeriods = 25;
end

if ~isfield(options, 'isPlot')
    options.isPlot = false;
end

%% UNDOCUMENTED Options:

if ~isfield(options, 'division1')
    options.division1 = [3 5 8 11 13];
end

if ~isfield(options, 'division2')
    options.division2 = [5 8 11];
end

if ~isfield(options, 'isAccelerated')
    options.isAccelerated = true;
end

if ~isfield(options, 'idR1')
    options.idR1 = 1;
end

if ~isfield(options, 'shift')
    options.shift = 0.5;
end

if ~isfield(options, 'shiftCriterion')
    options.shiftCriterion = 0.1;
end

if ~isfield(options, 'shiftNow')
    options.shiftNow = false;
end

if ~isfield(options, 'nPoly')
    options.nPoly = 2^5;
end

if ~isfield(options, 'iterate')
    options.iterate = false;
end   

if ~isfield(options, 'tol')
    options.tol = 1e-10;
end   

if ~isfield(options, 'verbose')
    options.verbose = false;
end

%% ALGORITHM
szx = size(x);
x   = x(:);
nx  = length(x);

% Set the proper integrand function for computiong the PDF
funPDF = @(t,x) real(cf(t).*exp(-1i*t*x));
funPDFshifted = @(t,x,shift) real(exp(1i*t*shift).*cf(t).*exp(-1i*t*x));

isPlot    = options.isPlot;
division1 = options.division1;
division2 = options.division2;
pdf = zeros(nx,1);
error = zeros(nx,1);
for id = 1:nx
    if abs(x(id)) < options.shiftCriterion || options.shiftNow
        isShifted = true;
        shift = options.shift;
        xshifted = x(id)+shift;
    else
        isShifted = false;
        shift = 0;
        xshifted = x(id);
    end
    if isShifted
        fun = @(t) funPDFshifted(t,xshifted,shift);
    else
        fun = @(t) funPDF(t,xshifted);
    end
    A = (pi/2)/abs(xshifted);
    B = pi*options.nPeriods/abs(xshifted);
    try
        [roots,warnings] = FindChebyRoots(fun,A,B,options.nPoly,isPlot);
    catch
        %warning('Problem using FindChebyRoots. Set roots = [].');
        roots = [];
        warnings = 0;
    end
    nRoots = length(roots);
    if isempty(roots)
        R1 = A; Rn = A;
        INT2 = 0; ERR2 = 0; Q  = 0;
    elseif numel(roots) == 1
        R1 = roots(1);
        Rn = roots(end);
        INT2 = 0; ERR2 = 0; Q  = 0;
    else
        if options.idR1 < numel(roots)
            R1 = roots(options.idR1);
            Intervals  =  [roots(options.idR1:end-1)';roots(options.idR1+1:end)'];
        else
            R1 = roots(1);
            Intervals  =  [roots(1:end-1)';roots(2:end)'];
       end
        Rn = roots(end);        
        [INT2,ERR2,~,Q] = IntegralGK(@(t) fun(t),Intervals, ...
            division2,isPlot);
    end
    Intervals = [0;R1];
    if options.iterate
        ERR1 = 1;
        while ERR1 > options.tol
            [~,ERR1,Intervals] = ...
                IntegralGK(@(t) fun(t),Intervals,division1,0);
        end
    end
    [INT1,ERR1] = IntegralGK(@(t) fun(t),Intervals,division1,isPlot);
    [INT3,ERR3] = IntegralGK(@(t) fun(t),[Rn;B],division1,isPlot);
    if options.isAccelerated  && warnings == 0 && nRoots > 10
        INT  = INT1 + altsum(sign(Q(1))*abs(Q));
        isAcceleration = true;
    else
        INT  = INT1 + INT2 + INT3;
        isAcceleration = false;
    end
    pdf(id)    = max(eps(0),INT/pi);
    error(id) = ERR1 + ERR2;
end
pdf    = reshape(pdf,szx);
x      = reshape(x,szx); 
error  = reshape(error,szx);
% Stop the clock
tictoc = toc;

%% RESULTS
if nargout > 1
    result.PDF = pdf;
    result.x = x;
    if options.verbose
        result.ErrorBound = error;
        result.ERR1_xLast = ERR1;
        result.ERR2_xLast = ERR2;
        result.ERR3_xLast = ERR3;
        result.warnings_xLast = warnings;
        result.isAcceleration_xLast = isAcceleration;
        result.fun_xLast = fun;
        if isShifted
            result.isShifted = isShifted;
            result.shift = shift;
            result.xshifted = xshifted;
            result.funPDFshifted = funPDFshifted;
        end
    end
    result.cf = cf;
    result.funPDF = funPDF;
    result.options = options;
    result.tictoc = tictoc;
end

end
%% Function IntegralGK
function [INT,ERR,SubIntervals,Q,Q1] = IntegralGK(fun,Intervals,division,isPlot)
%INTEGRALGK (G7,K15)-Gauss-Kronod quadrature over all Sub-Intervals of
%   Intervals defined by the index vector divisionIdx (a vector with
%   indices from 1 to 15, e.g. [3 5 8 11 13]
%
% SYNTAX
%   [INT,ERR] = IntegralGK(fun,Intervals,divisionIdx)
%
% EXAMPLE
%   [INT,ERR,SubInts] = IntegralGK(@(x) sin(x).*cos(5*x),[-1 0;0 1],[5 8 11])

% (c) 2015, Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 17-Apr-2015 09:53:32

%% CHECK/SET THE INPUT PARAMETERS
narginchk(2, 4);
if nargin < 4, isPlot   = []; end
if nargin < 3, division = []; end

if isempty(isPlot)
    isPlot = false;
end
%% ALGORITHM
[nodeGK,weightGK,WG,G] = GKnodes;
[SubIntervals,t,mids] = GetSubs(Intervals,nodeGK,division);
F   = fun(t);
Q1  = mids.*sum(bsxfun(@times,weightGK,F));
Q2  = mids.*sum(bsxfun(@times,WG,F(G,:)));
Q   = sum(reshape(Q1,length(division)+1,[]));
ERR = sum(abs(Q1-Q2),2);
INT = sum(Q);

if isPlot
    %LW = 'linewidth';
    figure
    plot(t,F,'.-')
    xlabel('x')
    ylabel('function')
    grid on
end

end
%% Function GetSubs
function [SubIntervals,nodes,mids] = GetSubs(Intervals,nodeGK,nodeIdx)
% GETSUBS Sub-division of the integration intervals for adaptive
%         Gauss-Kronod quadrature
%
% EXAMPLE
%  nodeGK =  GKnodes;
%  % Set the indices from 1:15 of nodes used for division of subintervals
%  nodeIdx = [5 8 11];
%  Intervals = [0;1]; % Given as 2 x n  matrix of initial sub-intervals
%  [SubIntervals,nodes,mids] = GetSubs(Intervals,nodeGK,nodeIdx)

% (c) 2015, Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 17-Apr-2015 09:53:32

mids = 0.5*(Intervals(2,:)-Intervals(1,:));
C = 0.5*(Intervals(2,:)+Intervals(1,:));
nodes = nodeGK*mids + ones(size(nodeGK))*C;

L = [Intervals(1,:); nodes(nodeIdx,:)];
U = [nodes(nodeIdx,:); Intervals(2,:)];
SubIntervals = [reshape(L, 1, []); reshape(U, 1, [])];

mids = 0.5*(SubIntervals(2,:)-SubIntervals(1,:));
C = 0.5*(SubIntervals(2,:)+SubIntervals(1,:));
nodes = nodeGK*mids + ones(size(nodeGK))*C;
end
%% Function GKnodes
function [XK,WK,WG,G] =  GKnodes
%GKNODES The (7-15) Gauss-Kronod nodes and weights

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 20-May-2013 01:28:00

nodes = [ ...
    0.2077849550078984676006894; 0.4058451513773971669066064; ...
    0.5860872354676911302941448; 0.7415311855993944398638648;...
    0.8648644233597690727897128; 0.9491079123427585245261897; ...
    0.9914553711208126392068547];
wt = [ ...
    0.2044329400752988924141620; 0.1903505780647854099132564; ...
    0.1690047266392679028265834; 0.1406532597155259187451896;
    0.1047900103222501838398763; 0.0630920926299785532907007;...
    0.0229353220105292249637320];
wt7 = [0.3818300505051189449503698; ...
    0.2797053914892766679014678; 0.1294849661688696932706114];

XK = [-nodes(end:-1:1); 0; nodes];
WK = [wt(end:-1:1); 0.2094821410847278280129992; wt];
WG = [wt7(end:-1:1); 0.4179591836734693877551020; wt7];
G = (2:2:15)';

end
%% Function FindChebyRoots
function [roots,warning,err] = FindChebyRoots(fun,A,B,n,isplot)
%FINDCHEBYROOTS estimates the real roots (zeros) of a real (oscilatory)
% function FUN on the interval [A,B], by using adaptive nth-order (n=2^k)
% Chebyshev polynomial approximation of the function FUN.
%
% This code was adapted from the code published in Day & Romero (2005):
% Roots Of Polynomials Expressed In Terms Of Orthogonal Polynomials. SIAM
% Journal on Numerical Analysis,  43, 1969 - 1987.
%
% SYNTAX:
% roots = FindChebyRoots(fun,A,B)
% roots = FindChebyRoots(fun,A,B,n,isplot)
%
% EXAMPLE 1
% fun = @(t) sin(t.^3 + t.^2 + t)
% A   = 0; 
% B   = 5;
% roots = FindChebyRoots(fun,A,B)
%
% EXAMPLE 2
% fun = @(t)  exp(-.3*t) .* sin(10*t) .* cos(2*t.^2)
% A   = 5; 
% B   = 15;
% roots = FindChebyRoots(fun,A,B)
%
% EXAMPLE 3
% x   = 10;
% nu  = 1;
% cf_chi2 = @(t) (1 - 2i * t) .^(-nu/2);
% fun = @(t) min(4,imag(exp(-1i*t*x).*cf_chi2(t))./t)
% A   = 0.2; 
% B   = 10;
% roots = FindChebyRoots(fun,A,B)
%
% EXAMPLE 4
% nu  = 3;
% fun = @(t)  sin(0.5+5*t) .* (besselj(nu,t) - besselk(nu,t))
% A   = 150; 
% B   = 200;
% roots = FindChebyRoots(fun,A,B)
%
% References
%   Day, D., Romero, L. (2005): Roots Of Polynomials Expressed In Terms Of
%   Orthogonal Polynomials. SIAM Journal on Numerical Analysis, 43,
%   1969 - 1987.

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 12-Aug-2013 16:18:42

%% SET THE DEFAULT VALUES (input parameters)
narginchk(1, 5);
if nargin < 5, isplot = []; end
if nargin < 4, n = []; end
if nargin < 3, B = []; end
if nargin < 2, A = []; end

%% SET THE DEFAULT VALUES (input parameters)
if isempty(isplot)
    isplot = true;
end

if isempty(n)
    n = 2^5;
end

if isempty(B)
    B = 1;
end

if isempty(A)
    A = -1;
end

interval = [A;B];

if abs(B-A) > 3*pi
    divisionrule = [-0.5;0;0.5];
    interval =  GetSubs4RF(interval,divisionrule);
else
    divisionrule = 0;
end

%% Estimate the roots of the real function fun over (A,B)
[cgl, M] = setupChebyshev(n);

roots = [];
warning = 0;
while true
    [r,interval,err,isWarning] = rootFinder(cgl,M,fun,interval,n);
    roots = cat(1,roots,r);
    if isempty(interval)
        break
    else
        interval =  GetSubs4RF(interval,divisionrule);
    end
    warning = warning + isWarning;
end
roots = unique(roots);

%% Plot the function together with the roots
if isplot
    N = 1000;
    t = linspace(A,B,N);
    plot(t,fun(t)),grid
    hold on
    plot(roots,fun(roots),'or')
    hold off
end
end % End of function FindChebyRoots
%% Function setupChebyshev
function [cgl,M] = setupChebyshev(n)
% setupChebyshev evaluates the Chebyshev-Gauss-Legendre points and the
% Chebyshev Transformation Matrix
%
% Adapted from the MATLAB code published in Day & Romero (2005): Roots of
% polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
% Numerical Analysis,  43, 1969 - 1987.

% Chebyshev-Gauss-Legendre points
ind = (0:n);
cgl = cos(ind * pi/n);

% Chebyshev Transformation Matrix
M        = cos(ind'* ind * pi/n);
M(1,:)   = M(1,:)/2;
M(:,1)   = M(:,1)/2;
M(n+1,:) = M(n+1,:)/2;
M(:,n+1) = M(:,n+1)/2;
M        = M * (2/n);
end % End of function setupChebyshev
%% Function evalCheb
function V = evalCheb(n,z)
%EVALCHEB evaluates the required Vandermonde matrix
%
% Input: vector of points, z, and the polynomial degree n.
% Output: Vandermonde matrix, m by degree_max + 1, V(j+1,k)=T_j(z_k)
%
% Adapted from the MATLAB code published in Day & Romero (2005): Roots of
% polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
% Numerical Analysis,  43, 1969 - 1987.

m = size(z,1);
if m*n >= 0
    V(:,1) = ones(m,1);
    if n >= 1
        V(:,2) = z;
        if n >= 2
            index = find( log(abs(z)) >= 100/n );
            si = size(index,1);
            if si > 0
                z(index) = ones(si,1)*exp(100/n);
            end
            for i=2:n
                V(:,i+1) = V(:,i).*(2*z) - V(:,i-1);
            end
        end
    end
else
    V = [];
end
end % End of function evalCheb
%% Function rootFinder
function [roots,intervals,err,isWarning] = rootFinder(cgl,M,fun,intervals,n)
%ROOTFINDER estimates the roots of fun over intervals by using nth order
% Chebyshev polynomial approximation of fun
%
% Adapted from the MATLAB code published in Day & Romero (2005): Roots of
% polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
% Numerical Analysis,  43, 1969 - 1987.

% Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 12-Aug-2013 16:18:42

roots = [];

nint = size(intervals,2);
err  = false(nint,1);
isWarning = 0;
for i = 1:nint
    range = (intervals(2,i) - intervals(1,i))/2;
    shift = (intervals(1,i) + intervals(2,i))/2;
    x = cgl*range + shift;
    f = fun(x);
    ExpansionCoeff = f * M;
    
    if abs(ExpansionCoeff(n+1)) <  1000*eps(0)
        isWarning = 1;
        % warning('The leading expansion coefficient vanishes');        
    else
        ExpansionCoeff = ExpansionCoeff/(-2*ExpansionCoeff(n+1));
        H = diag(ones(n-1, 1)/2, 1) + diag(ones(n-1, 1)/2, -1);
        H(1, 2) = 1;
        C = H;
        C(n, :) = C(n, :) + ExpansionCoeff(1:n);
        
        Eigenvalues = eig(C);
        Vandermonde = evalCheb(n,Eigenvalues);
        Vcolsums = sum( abs(Vandermonde') );
        tube_index = find((abs(imag(Eigenvalues))<.2) ...
            & (abs(real(Eigenvalues))< 2));
        Solutions = Eigenvalues( tube_index );
        Vcolsums = Vcolsums( tube_index );
        cond_max = min( 2^(n/2), 10^6 );
        condEigs_index =  Vcolsums < cond_max ;
        Solutions = sort(Solutions( condEigs_index ));
        r = range * Solutions + shift;
        if isreal(r)
            r = r( r>=intervals(1,i) & r<=intervals(2,i) );
            roots = cat(1,roots,r);
        else
            err(i) = true;
        end
    end
end
intervals = intervals(:,err);
end % End of function rootFinder
%% Function GETSUBS (Sub-division of the integration intervals)
function SUBS = GetSubs4RF(SUBS,XK)
%GetSubs4RF Sub-division of the intervals for adaptive root finding
%
% Example
% interval = [0;1];
% divisionrule = [-.25; -.5; 0; .5; .75] % selected points from [-1;1]
% intervals = GetSubs4RF(interval,divisionrule)
% intervals = GetSubs4RF(intervals,divisionrule)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 20-May-2013 01:28:00

NK = size(XK,1);
SUBIND = 1:NK;
M = 0.5*(SUBS(2,:)-SUBS(1,:));
C = 0.5*(SUBS(2,:)+SUBS(1,:));
Z = XK*M + ones(NK,1)*C;

L = [SUBS(1,:); Z(SUBIND,:)];
U = [Z(SUBIND,:); SUBS(2,:)];
SUBS = [reshape(L, 1, []); reshape(U, 1, [])];
end % End of function GetSubs
%% Function altsum
function s = altsum(a,method)
%ALTSUM Convergence acceleration of alternating series.
%	ALTSUM(A) where A is a vector, estimates the sum of the infinite 
%   alternating series A(1) - A(2) + A(3) - A(4) + ...
%
%   ALTSUM(A,METHOD) where A is a vector of length N and METHOD is a 
%   string selects one of the following linear acceleration methods:
%   'euler', Euler's method with error ~ 2^-N
%   'binom', a binomial weighting with error ~ 3^-N
%   'cheby', a weighting based on Chebyshev polynomials, error ~ 5.828^-N
%   Default method is 'cheby'.
%
%   EXAMPLES
%
%	% log(2) = 1 - 1/2 + 1/3 - 1/4 + ...
%	altsum(1./(1:20))
%
%	% pi = 4 - 4/3 + 4/5 - 4/7 + ...
%	altsum(4./(1:2:39))
%
%   % Euler-Mascheroni constant
%   k = 2:20;
%   altsum(log(k)./k)/log(2) + log(2)/2
%
%   % Divergent series
%   altsum(1:9,'euler')
%   altsum((1:9).^2,'euler')

%   Author: Jonas Lundgren <splinefit@gmail.com> 2009

%   2009-09-07  First published.
%   2009-09-16  Two acceleration methods added.
%
% Source: 
%   http://www.mathworks.co.uk/matlabcentral/fileexchange/25200-altsum
%
% ---------------------------------------------------------------------

if nargin < 1, help altsum, return, end
if nargin < 2, method = 'cheby'; end

% Index vector
n = numel(a);
k = (1:n)';

% Select method
switch method
    case 'euler'
        % Euler's method ~ 2^-n
        p = 1/2;
        Q = (n-k+1)./k;
    case 'binom'
        % Binomial weighting ~ 3^-n
        p = 2/3;
        Q = (n-k+1)./k/2;
    otherwise
        % Chebyshev weighting ~ 5.828^-n
        p = 4/(3+sqrt(8));
        Q = (n-k+1)./k.*(n-k+1/2)./(2*n-k);
end

% Distribution
d = exp(cumsum([n*log(p); log(Q)]));        % d = cumprod([p^n; Q]);
w = cumsum(d);

% Alternate signs
a = a(:);
a(2:2:n) = -a(2:2:n);

% Reverse order
a = a(n:-1:1);

% Weighted sum
s = sum(w(1:n).*a)/w(n+1);
end