function cf = cf_FConv(t,coefs,DF1,DF2,n_iid,options)
%%CF_FCONV Characteristic function (CF) of the distribution of a linear
%  combination of independent central F random variables, i.e. of the
%  distribution of Y = sum_{i=1}^N coef_i * F{DF1_i,DF2_i}, where DF1_i and
%  DF2_i represent the degrees of freedom of the Fisher-Snedecor
%  F-distribution.
%
% SYNTAX
%  cf = cf_FConv(t,coefs,DF1,DF2,n_iid)
%
% EXAMPLE1 (CF of a linear combination of independent F RVs)
%  N = 501;
%  t = linspace(-10,10,N);
%  coefs = [1 2];
%  DF1   = [5 10];
%  DF2   = [7 5];
%  cf = cf_FConv(t,coefs,DF1,DF2);
%  figure; plot(t,real(cf),t,imag(cf))
%  title('Characteristic function of the linear combination of F RVs')
%
% EXAMPLE2 (Compute PDF/CDF from the CF by CF2DIST [required])
%  coefs = 1;
%  DF1   = 5;
%  DF2   = 7;
%  cf = @(t) cf_FConv(t,coefs,DF1,DF2);
%  %options.sixSigmaRule = 10;
%  options.isForcedSymmetric = true;
%  result = cf2DIST(cf,options);
%
% EXAMPLE3 (Compute PDF/CDF from the CF by CF2DIST [required])
%  coefs = [1 2 3 4 5 6];
%  DF1   = [5 10 3.5 6 5 5];
%  DF2   = [7 5 4 5 6 7];
%  cf = @(t) cf_FConv(t,coefs,DF1,DF2);
%  %options.n = 2^8;
%  %options.xmax = 1e4;
%  %options.sixSigmaRule = 10;
%  options.isForcedSymmetric = true;
%  result = cf2DIST(cf,options);
%
% REFERENCES
% https://en.wikipedia.org/wiki/F-distribution

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Mar-2016 13:25:54

%% CHECK THE INPUT PARAMETERS
narginchk(1, 6);
if nargin < 6, options = []; end
if nargin < 5, n_iid = []; end
if nargin < 4, DF2 = []; end
if nargin < 3, DF1 = []; end
if nargin < 2, coefs = []; end

if isempty(DF2) && ~isempty(DF1)
    DF2 = 1;
elseif isempty(DF2) && ~isempty(coefs)
    DF2 = 1;
elseif ~any(DF2)
    DF2 = 1;
end

if isempty(DF1) && ~isempty(coefs)
    DF1 = 1;
elseif isempty(DF1) && ~isempty(DF2)
    DF1 = 1;
end

if isempty(coefs) && ~isempty(DF2)
    coefs = 1;
elseif isempty(coefs) && ~isempty(DF1)
    coefs = 1;
end

%% SET THE COMMON SIZE of the parameters
[errorcode,coefs,DF1,DF2] = distchck(3,coefs(:),DF1(:),DF2(:));
if errorcode > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function of a linear combination of independent F RVs
% Set the integration method for computing the CF:
options.method    = 'standard';
options.precision = 'standard';
options.digits    = 16;

szt = size(t);
t   = t(:);
cf  = cfFun(DF1(1),DF2(1),coefs(1)*t,options);
for i = 2:length(coefs)
    cf = cf .* cfFun(DF1(i),DF2(i),coefs(i)*t,options);
end

cf = reshape(cf,szt);
cf(t==0) = 1;

if isscalar(coefs) && isscalar(DF1) && isscalar(DF2) && isscalar(n_iid)
    cf = cf .^ n_iid;
end
end

%% FUNCTION CFFUN
function [int,err,xMax,xLow,xUpp,fun] = cfFun(DF1,DF2,t,options)
% CFFUN integrates the integrand function fun to compute the characteristic
%  function of F distribution. CFFUN evaluates the integrand function,
%  evaluates the optimum integration limits [xLow,xUpp], and computes the
%  integral by the (non-adaptive resp. adaptive) Gauss-Kronod quadrature.
%
% SYNTAX
%  [int,err,xMax,xLow,xUpp,fun] = cfFun(DF1,DF2,t,options)
%
% EXAMPLE
%  [Int,Err,xMax,L,U,fun] = cfFun(10,5,1)
%  x = linspace(L,U,501);
%  plot(x,real(fun(x)),'.-',x,imag(fun(x)))

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Mar-2016 13:25:54

%% CHECK THE INPUT PARAMETERS
narginchk(3, 4);
if nargin < 4, options = []; end

% Set the integration method: 'standard', 'quadgk', 'integral'
if ~isfield(options,'method')
    options.method = 'standard';
end

% Set the integration precision for standard: 'standard', 'medium', 'high'
if ~isfield(options,'precision')
    options.precision = 'standard';
end

% Set the integration precision for quadGK: 'RelTol', 1e-10
if ~isfield(options,'RelTol')
    options.RelTol = 1e-8;
end

% Set the required significant digits: digits = 16 or log(1/eps)
if ~isfield(options,'digits')
    options.digits =  16;
end

digits = options.digits;
precision = options.precision;
method = options.method;

%% EVALUATE THE INTEGRAND FUNCTION / fun resp. log(fun)
%szt  = size(t);
a   = DF1/2;
b   = 1-DF2/2;
z   = -DF2/DF1*t(:);
abz = abs(z);
szz = size(z);
B   = a-1;
C   = a-b+1;
D   = B-C;
E   = 1-1i*sqrt(3);
Fun   = 1+1i*sqrt(3);
Gauss   = 2^(2/3);
H   = 2^(1/3);
K   = 2*B+C;
L   = D^2-3*z.^2;
M   = sqrt(27*(4*B*D^3*z.^2 +(8*B^2+20*B*C-C^2)*z.^4 +4*z.^6));
N   = (-2*D^3 -9*K*z.^2 +M).^(1/3);
x1  = (2*D -2*H*L./N -Gauss*N)/6;
x2  = (4*D +2*Fun*H*L./N +E*Gauss*N)/12;
x3  = (4*D +4*E*L./(Gauss*N) +2*Fun*N/H)/12;

xMax  = max([zeros(szz) real(x1) real(x2) real(x3)],[],2);
xLow  = zeros(szz);
xUpp  = zeros(szz);
crit  = ones(szz);
crit0 = ones(szz);
int   = zeros(szz);
err   = zeros(szz);
x0    = 1e-10;

%% abs(z) >= 1
idx = find(abz>=1);
if any(idx)
    zz  = abz(idx);
    xm  = max(xMax(idx),x0);
    crit(idx) = digits + xm + C*log(1+(xm./zz).^2)/2 - ...
        log(1./zz) - B*log(xm./zz);
    % xupp
    x1 = xm;
    f1 = digits;
    x2 = xMax(idx) + 1e4;
    for i = 1:10 % Bisection method
        x = (x1+x2)/2;
        fx = crit(idx) - x - C*log(1+(x./zz).^2)/2 + ...
            log(1./zz) + B*log(x./zz);
        if norm(x2-x1) < 1e-1
            return
        end
        id = find(sign(fx)==sign(f1));
        x1(id) = x(id);
        id = find(sign(fx)~=sign(f1));
        x2(id)= x(id);
        f1 = crit(idx) - x1 - C*log(1+(x1./zz).^2)/2 + ...
            log(1./zz) + B*log(x1./zz);
    end
    for i = 1:3 % Halley's method
        fx   = crit(idx) - x - C*log(1+(x./zz).^2)/2 + ...
            log(1./zz) + B*log(x./zz);
        fD1 = -1 +B./x -C*x./(zz.^2+x.^2);
        fD2 = -B./x.^2 + 2*C*x.^2./(zz.^2+x.^2).^2 - C./(zz.^2+x.^2);
        x   = max(xMax(idx),x -2*fx.*fD1./(2*fD1.^2-fx.*fD2));
    end
    if x == xMax(idx)
        xUpp(idx) = xMax(idx) + 1e4;
    else
        xUpp(idx) = x;
    end
    % xlow
    crit0(idx) = crit(idx) - x0 - C*log(1+(x0./zz).^2)/2 + ...
        log(1./zz) + B*log(x0./zz);
    idx = find(abz>=1 & xMax ~= 0 & xMax ~= inf & crit0 <=0);
    zz  = abz(idx);
    x   = xMax(idx)/2;
    for i = 1:3 % Halley's method
        fx   = crit(idx) - x - C*log(1+(x./zz).^2)/2 + ...
            log(1./zz) + B*log(x./zz);
        fD1 = -1 +B./x - C*x./(zz.^2+x.^2);
        fD2 = -B./x.^2 + 2*C*x.^2./(zz.^2+x.^2).^2 - C./(zz.^2+x.^2);
        x   = max(0,x -2*fx.*fD1./(2*fD1.^2-fx.*fD2));
    end
    xLow(idx) = x;
end

%% abs(z) < 1
idx = find(abz<1);
if any(idx)
    zz  = abz(idx);
    xMax(idx) = xMax(idx)./zz;
    xm  = max(xMax(idx),1e-10);
    crit(idx) = digits + xm + C*log(1+(xm./zz).^2)/2 - ...
        log(1./zz) - B*log(xm./zz);
    zz  = abz(idx);
    % xupp
    x1 = xm;
    f1 = digits;
    x2 = xMax(idx) + 1e4;
    for i = 1:10 % Bisection method
        x = (x1+x2)/2;
        fx = crit(idx) -x.*zz -C*log(1+x.^2)/2 +B*log(x);
        if norm(x2-x1) < 1e-1
            return
        end
        id = find(sign(fx)==sign(f1));
        x1(id) = x(id);
        id = find(sign(fx)~=sign(f1));
        x2(id)= x(id);
        f1 = crit(idx) -x1.*zz -C*log(1+x1.^2)/2 +B*log(x1);
    end
    for i = 1:3 % Halley's method
        fx = crit(idx) -x.*zz -C*log(1+x.^2)/2 +B*log(x);
        fD1 = -zz +B./x -C*x./(1+x.^2);
        fD2 = -B./x.^2 +2*C*x.^2./(1+x.^2).^2 -C./(1+x.^2);
        x  = max(xMax(idx),x -2*fx.*fD1./(2*fD1.^2-fx.*fD2));
    end
    if x == xMax(idx)
        xUpp(idx) = xMax(idx) + 1e4;
    else
        xUpp(idx) = x;
    end
    % xlow
    crit0(idx) = crit(idx) -x0 -C*log(1+(x0./zz).^2)/2 +log(1./zz) +B*log(x0./zz);
    idx = find(abz<1 & xMax ~= 0 & xMax ~= inf & crit0 <=0);
    zz  = abz(idx);
    x   = xMax(idx)/2;
    for i = 1:3
        fx = crit(idx) -x.*zz -C*log(1+x.^2)/2 +B*log(x);
        fD1 = -zz +B./x -C*x./(1+x.^2);
        fD2 = -B./x.^2 +2*C*x.^2./(1+x.^2).^2 -C./(1+x.^2);
        x  = max(0,x- 2*fx.*fD1./(2*fD1.^2-fx.*fD2));
    end
    xLow(idx) = x;
end

%% Integrate
const  = (gammaln(a-b+1) -gammaln(1-b) -gammaln(a));
switch method
    case lower({'standard','GK','nonadaptive'})
        [nodeGK,weightGK,WGauss,Gauss] = GKnodes;
        [Intervals,nod,mid] = GetSubs([-1;1],nodeGK,1:15);
        switch precision
            case lower('high')
                [~,nod,mid] = GetSubs(Intervals,nodeGK,1:15);
            case lower('medium')
                [~,nod,mid] = GetSubs(Intervals,nodeGK,[3 8 13]);
        end
        nnods = length(nod);
        for i = 1:szz(1)
            nods = [(xLow(i)+xMax(i))/2 + nod*(xMax(i)-xLow(i))/2, ...
                (xMax(i)+2*xMax(i))/2 + nod*(xMax(i)-2*xMax(i))/2, ...
                (2*xMax(i)+4*xMax(i))/2 + nod*(2*xMax(i)-4*xMax(i))/2, ...
                (4*xMax(i)+xUpp(i))/2 + nod*(xUpp(i)-4*xMax(i))/2];
            mids = [mid*(xMax(i)-xLow(i))/2, mid*(2*xMax(i)-xMax(i))/2,...
                mid*(4*xMax(i)-2*xMax(i))/2, mid*(xUpp(i)-4*xMax(i))/2];
            if abs(z(i)) >= 1
                d = -1i/z(i);
            else
                if z(i) >= 0
                    d = -1i;
                else
                    d = 1i;
                end
            end
            Fun = exp(const+log(d)+(a-1)*log(d*nods)- ...
                (a-b+1)*log(1+d*nods)-1i*d*z(i)*nods);
            Q1  = mids.*sum(bsxfun(@times,weightGK,Fun));
            Q2  = mids.*sum(bsxfun(@times,WGauss,Fun(Gauss,:)));
            Q   = sum(reshape(Q1,nnods,[]));
            ERR = sum(abs(Q1-Q2),2);
            INT = sum(Q);
            int(i) = INT;
            err(i) = ERR;
        end
    case lower({'adaptive','quadGK' 'integral'})
        for i = 1:szz(1)
            if abs(z(i)) >= 1
                d = -1i/z(i);
            else
                if z(i) >= 0
                    d = -1i;
                else
                    d = 1i;
                end
            end
            if strcmp(method,'quadGK')
                [INT,ERR] = quadgk(@(x)exp(const+log(d)+(a-1)*log(d*x)- ...
                    (a-b+1)*log(1+d*x)-1i*d*z(i)*x), ...
                    xLow(i),xUpp(i),'RelTol',1e-10);
            else
                [INT,ERR] = itegral(@(x)exp(const+log(d)+(a-1)*log(d*x)- ...
                    (a-b+1)*log(1+d*x)-1i*d*z(i)*x), ...
                    xLow(i),xUpp(i),'RelTol',1e-10);
            end
            int(i) = INT;
            err(i) = ERR;
        end
end

%% abs(z)==0
idx = find(abz==0);
if any(idx)
    int(idx) = 1;
    err(idx) = 0;
end

%% fun / Function handle of the last function
const  = (gammaln(a-b+1) -gammaln(1-b) -gammaln(a));
if abs(z) >= 1
    d = -1i/z;
else
    if z >= 0
        d = -1i;
    else
        d = 1i;
    end
end
fun = @(x) exp(const +log(d) +(a-1)*log(d*x) -(a-b+1)*log(1+d*x) -1i*d*z*x);
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