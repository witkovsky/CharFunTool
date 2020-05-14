function result = PolyCal(x,y,xfit,yLow,yUpp,options)
% PolylCal Comparative polynomial calibration. PolyCal estimates the
% patrameters of the polynomial calibration function together with 
% their uncertainties. 
%
% SYNTAX
% result = PolyCal(x,y,fitXpts,yLow,yUpp,options)
% 
% EXAMPLE
% load('CalData.mat')
% clear options
% options.Aijk  = Aijk;
% options.order = 3;
% options.uX    = uX;
% options.uX0   = uX0;
% options.uY    = uY;
% options.uY0   = uY0;
% options.cfX   = cfX;
% options.cfX0  = cfX0;
% options.cfY   = cfY;
% options.cfY0  = cfY0;
% xfit          = linspace(min(x),max(x),11);
% id            = 4;
% yLow          = y(id) - 2*sqrt((10*uY(id))^2+uY0^2);
% yUpp          = y(id) + 2*sqrt((10*uY(id))^2+uY0^2);
% result = PolyCal(x,y,xfit,yLow,yUpp,options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 14-May-2020 12:58:43

%% Check the inputs
narginchk(2, 6);

if nargin < 6, options = []; end
if nargin < 5, yUpp    = []; end
if nargin < 4, yLow    = []; end
if nargin < 3, xfit    = []; end

if ~isfield(options, 'uX')
    options.uX = [];
end

if ~isfield(options, 'uX0')
    options.uX0 = [];
end

if ~isfield(options, 'uY')
    options.uY = [];
end

if ~isfield(options, 'uY0')
    options.uY0 = [];
end

if ~isfield(options, 'cfX')
    options.cfX = [];
end

if ~isfield(options, 'cfX0')
    options.cfX0 = [];
end

if ~isfield(options, 'cfY')
    options.cfY = [];
end

if ~isfield(options, 'cfY0')
    options.cfY0 = [];
end

if ~isfield(options, 'alpha')
    options.alpha = 0.05;
end

if ~isfield(options, 'probLow')
    options.probLow = [];
end

if ~isfield(options, 'probUpp')
    options.probUpp = [];
end


if ~isfield(options, 'isplot')
    options.isplot = true;
end

if ~isfield(options, 'maxiter')
    options.maxiter = 100;
end

if ~isfield(options, 'tolerance')
    options.tolerance = 1e-12;
end

if ~isfield(options, 'order')
    options.order = [];
end

if ~isfield(options, 'isXY')
    options.isXY = true; 
end

if ~isfield(options, 'Aijk')
    options.Aijk = [];
end

if ~isfield(options, 'q')
    options.q = [];
end

if ~isfield(options, 'mu0')
    options.mu0 = [];
end

if ~isfield(options, 'nu0')
    options.nu0 = [];
end

if ~isfield(options, 'par0')
    options.par0 = [];
end

if ~isfield(options, 'cfIntervals')
    options.cfIntervals = true;
end




%% Set the Inputs

isXY = options.isXY;
if isXY
    uX      = options.uX;
    uX0     = options.uX0;
    uY      = options.uY;
    uY0     = options.uY0;
    cfX     = options.cfX;
    cfX0    = options.cfX0;
    cfY     = options.cfY;
    cfY0    = options.cfY0;
    x       = x(:);
    xMin    = min(x);
    xMax    = max(x);
    y       = y(:);
    yMin    = min(y);
    yMax    = max(y);
else
    uX      = options.uY;
    uX0     = options.uY0;
    uY      = options.uX;
    uY0     = options.uX0;
    cfX     = options.cfY;
    cfX0    = options.cfY0;
    cfY     = options.cfX;
    cfY0    = options.cfX0;
    x       = y(:);
    xMin    = min(y);
    xMax    = max(y);
    y       = x(:);
    yMin    = min(x);
    yMax    = max(x);
end

order   = options.order;
Aijk    = options.Aijk;
q       = options.q;
alpha   = options.alpha;
probLow = options.probLow;
probUpp = options.probUpp;
isplot  = options.isplot;
mu0     = options.mu0;
nu0     = options.nu0;
pars0   = options.par0;
% method  = options.method;
isCFInt = options.cfIntervals;

n    = length(x);
o    = ones(n,1);

if isempty(uX), uX  = o; end
if isempty(uY), uY  = o; end
if isempty(uX0), uX0 = 0; end
if isempty(uY0), uY0 = 0; end

if isempty(probUpp)
    probUpp = 1-alpha/2;
end

if isempty(probLow)
    probLow = alpha/2;
end

if isempty(cfX)
    cfX  = cell(n,1); 
    for i = 1:n
        cfX{i} = @(t) exp(-(uX(i)*t).^2/2);
    end
end

if isempty(cfY)
    cfY  = cell(n,1); 
    for i = 1:n
        cfY{i} = @(t) exp(-(uY(i)*t).^2/2);
    end
end

if isempty(cfX0)
    cfX0 = @(t) exp(-(uX0*t).^2/2);
end

if isempty(cfY0)
    cfY0 = @(t) exp(-(uY0*t).^2/2);
end

SigX    = diag(uX.^2) + uX0^2 * ones(n);
SigY    = diag(uY.^2) + uY0^2 * ones(n);

ucX     = sqrt(diag(SigX));
ucY     = sqrt(diag(SigY));

if isempty(order)
    order = 1;
end

if isempty(xfit)
    xfit = linspace(xMin,xMax)';
end

if isempty(mu0)
    mu0 = x;
end

if isempty(nu0)
    nu0 = y;
end

if isempty(q)
    if ~isempty(Aijk)
        [~,~,q] = size(Aijk);
    end
end

if isempty(pars0)
    if isempty(Aijk)
        pars0 = flip(polyfit(mu0,nu0,order))';
    else
        B0 = update(mu0,nu0,pars0,n,order,q,Aijk);
        pars0 = B0\nu0;
    end
end

%% Algorithm

% Initialization
maxiter      = options.maxiter;
tolerance    = options.tolerance;
crit         = 1;
loops        = 1;
muEstimate   = mu0;
nuEstimate   = nu0;
parsEstimate = pars0;

while (crit > tolerance) && (loops <= maxiter)
    mu0         = muEstimate;
    nu0         = nuEstimate;
    pars0       = parsEstimate;
    muOld       = muEstimate;
    nuOld       = nuEstimate;
    parsOld     = parsEstimate;
    [B0,D0]     = update(mu0,nu0,pars0,n,order,q,Aijk);
    W0          = D0*SigX*D0 + SigY;
    W0inv       = W0\eye(n);
    Q0          = [W0 B0; B0' zeros(order+1)]\eye(n+order+1);
    Q011        = Q0(1:n,1:n);
    Y           = D0*(x-mu0) - y;
    muEstimate  = x - SigX*D0*Q011*Y;
    nuEstimate  = y + SigY*Q011*Y;
    parsEstimate = -(B0'*W0inv*B0)\(B0'*W0inv*Y);
    crit = norm([muEstimate;nuEstimate;parsEstimate] - ... 
        [muOld;nuOld;parsOld])^2/(2*n+order+1); 
    loops = loops+1;
end
[B,D,A,c] = update(muEstimate,nuEstimate,parsEstimate,n,order,q,Aijk);
W         = D*SigX*D + SigY;
Q         = [W B; B' zeros(order+1)]\eye(n+order+1);
Winv      = W\eye(n);
KY        = (B'*Winv*B)\(B'*Winv);
LX        = -KY * D;
id        = (n+1):(n+order+1);
SigPars   = -Q(id,id);
uPars     = sqrt(diag(SigPars));

% calculate the CFs of the polynomial parameters

cfPars = cell(order+1,1);
for i = 1:(order+1)
    wX = LX(i,:);
    wY = KY(i,:);
    cfPars{i} = WeightedCF(wX,wY,cfX,cfX0,cfY,cfY0,parsEstimate(i));
end

% Calculate the fitted values
xfit = xfit(:);
m = length(xfit);
w = ones(m,order+1);

for i = 2:(order+1)
    w(:,i) = w(:,i-1).*xfit;
end

yFit  = zeros(m,1);
uFit  = zeros(m,1);
cfFit = cell(m,1);
for j = 1:m
    yFit(j)  = w(j,:)*parsEstimate;
    uFit(j)  = sqrt(w(j,:)*SigPars*w(j,:)');
    wX       = w(j,:)*LX;
    wY       = w(j,:)*KY;
    cfFit{j} = WeightedCF(wX,wY,cfX,cfX0,cfY,cfY0,yFit(j));
end

if isCFInt
    cfInt = zeros(m,2);
    QF    = cell(m,1);
    prob  = [alpha/2,1-alpha/2];
    options.isInterp = true;
    options.isPlot = false;
    for j = 1:m
        res = cf2DistGP(cfFit{j},[],prob,options);
        cfInt(j,:) = res.qf;
        QF{j}      = res.QF;
    end
else
    cfInt = [];
    QF    = [];
end

xLow = [];
xUpp = [];
if ~isempty(yLow) && ~isempty(yUpp)
    nyLow = length(yLow);
    nyUpp = length(yUpp);
    xLow  = zeros(nyLow,1);
    xUpp  = zeros(nyUpp,1);
    if nyLow == nyUpp
        for i = 1:nyLow
            [xL,xU] = CalInt(yLow(i),yUpp(i),xMin,xMax,probLow,probUpp,...
                order,LX,KY,cfX,cfX0,cfY,cfY0,parsEstimate,options);
            xLow(i) = xL;
            xUpp(i) = xU;
        end
    else
        warning('The size of yLow does not fit with the size of yUpp')
    end
end


%% Results
result.pars    = parsEstimate';
result.uPars   = uPars';
result.SigPars = SigPars;
result.cfPars  = cfPars;
result.xFit    = xfit;
result.yFit    = yFit;
result.uFit    = uFit;
result.cfInt   = cfInt;
result.cfFit   = cfFit;
result.qfFit   = QF;
result.mu      = muEstimate;
result.nu      = nuEstimate;
result.yCI     = [yLow yUpp];
result.xCI     = [xLow xUpp];
result.LX      = LX;
result.KY      = KY;
result.Q       = Q;
result.x       = x;
result.xMin    = xMin;
result.xMax    = xMax;
result.y       = y;
result.yMin    = yMin;
result.yMax    = yMax;
result.SigX    = SigX;
result.SigY    = SigY;
result.ucX     = ucX;
result.ucY     = ucY;
result.n       = n;
result.order   = order;
result.q       = q;
result.Aijk    = Aijk;
result.A       = A;
result.B       = B;
result.c       = c;
result.D       = D;
result.W       = W;
result.crit    = crit;
result.tol     = tolerance;
result.loops   = loops;
result.uX      = uX;
result.uY      = uY;
result.uX0     = uX0;
result.uY0     = uY0;
result.cfX     = cfX;
result.cfY     = cfY;
result.cfX0    = cfX0;
result.cfY0    = cfY0;
result.isXY    = isXY;
result.options = options;

%% Plot

if isplot
    figure
    plot(x,y,'ko')
    hold on
    for i = 1:n
        plot([x(i)-2*ucX(i),x(i)+2*ucX(i)],[y(i),y(i)],'m-')
        plot([x(i),x(i)],[y(i)-2*ucY(i),y(i)+2*ucY(i)],'m-')
    end
    plot(muEstimate,nuEstimate,'x-')
    plot(xfit,yFit,'b+')
    if isempty(cfInt)
        plot(xfit,yFit+2*uFit,'r-.')
        plot(xfit,yFit-2*uFit,'r-.')
    else
        plot(xfit,cfInt(:,1),'r-.')
        plot(xfit,cfInt(:,2),'r-.')
        %
        plot(xfit,yFit+2*uFit,'y-.')
        plot(xfit,yFit-2*uFit,'y-.')
    end
    hold off
    grid on
end

end
%% Function update
function [B,D,A,C] = update(mu0,nu0,a0,n,order,q,Aijk)
% update calculates the update of the matrices A, B and the vector c,
% for a polynomial calibration model of given order, locally at the
% specified mu0, nu0, and a0. If the coefficient (n,order+1,q)-dimensional
% array Aijk is given, then the updated matrices A, B and the vector c are
% calculated for the generalized polynomial calibration model.
%
% SYNTAX
%  [B,D,A,C] = update(mu0,nu0,a0,n,order,q,Aijk)
%
% EXAMPLE
%  mu0   = [-19.6843   -9.8142    0.0989   10.0149   19.9634   29.8131]';
%  nu0   = [92.1489   96.0965  100.0499  103.9924  107.9354  111.8316]';
%  a0    = [100.0108    0.3982   -0.0001]';
%  n     = 6;
%  order = 2;
%  q     = [];
%  Aijk  = [];
% [B,D,A,C] = update(mu0,nu0,a0,n,order,q,Aijk)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Apr-2018 13:19:57

%%
if nargin < 7
    Aijk = []; 
end

if nargin < 6
    q = []; 
end

polyPars = flip(a0);
if isempty(Aijk)
    polyParsD1 = polyder(polyPars);
    D = diag(polyval(polyParsD1,mu0));
    A = [D -eye(n)];
    B = ones(n,order+1);
    for j = 1:order
        B(:,j+1) = mu0.*B(:,j);
    end
    C = polyval(polyPars,mu0) - nu0;
else
    B = zeros(n,order+1);
    D = zeros(n,order+1);
    for i = 1:n
        for j = 1:(order+1)
            for k = 1:q
                B(i,j) = B(i,j) + Aijk(i,j,k)*mu0(i)^(k-1);
                if k > 1
                    D(i,j) = D(i,j) + Aijk(i,j,k)*(k-1)*mu0(i)^(k-2);
                end
            end
        end
    end
    if isempty(a0)
        D = [];
        A = [];
        C = [];
    else
        D = diag(D*a0);
        A = [D -eye(n)];
        C = B*a0 - nu0;
    end
end
end
%% Function WeightedCF
function cf = WeightedCF(wX,wY,cfX,cfX0,cfY,cfY0,shift)
% WeightedCF creates the combined 'calibration' characteristic function of
% the random variable W = sum_{i=1}^n (wX_i*(X_i+X_0) + wY_i*(Y_i+Y_0)), 
% defined by 
% cfW = @(t)(prod_{i=1}^n cfX{i}(wX(i)*t)).*cfX0((sum_{i=1}^I wX(i))*t) ...
%           (prod_{i=1}^n cfY{i}(wY(i)*t)).*cfY0((sum_{i=1}^I wY(i))*t), 
% where wX and wY are n-dimensional vectors of coefficients, and cfX =
% {@(t)cfX{1}(t), ..., @(t)cfX{n}(t)} and cfY = {@(t)cfY{1}(t), ...,
% @(t)cfY{n}(t)} are n-dimensional cell vectors of the characteristic
% functions of (X_1,...,X_n) and (Y_1,...,Y_n), respectively, and cfX0 and
% cfY0 are the characteristic functions of the 'common' random variables
% X_0 and Y_0, respectively.
%
% SYNTAX
% cf = WeightedCF(wX,wY,cfX,cfX0,cfY,cfY0)
%
% EXAMPLE
%  wX = [1 2 3 4 5]/15;
%  wY = [1 1 1 1 1]/5;
%  cfX  = {@(t)cfS_Gaussian(t), @(t)cfS_Rectangular(0.5*t), ...
%          @(t)cfS_Rectangular(1.2*t), @(t)cfS_Rectangular(2*t), ...
%          @(t)cfS_Arcsine(8*t)}';
%  cfX0 =  @(t)cfS_Rectangular(1.5*t);
%  cfY  =  {@(t)cfS_StudentT(0.5*t,10), @(t)cfS_Triangular(0.5*t), ...
%          @(t)cfS_Triangular(1.2*t), @(t)cfS_Triangular(2*t), ...
%          @(t)cfS_Triangular(2*t)}';
%  cfY0 =  @(t)cfS_Arcsine(3*t);
%  cf = WeightedCF(wX,wY,cfX,cfX0,cfY,cfY0)
%  figure
%  t = linspace(-5,5,201);
%  plot(t,real(cf(t)))
%  figure 
%  x      = linspace(-10,10,201);
%  prob   = [0.9 0.95 0.99];
%  result = cf2DistGP(cf,x,prob);
%  disp(result)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Apr-2018 13:19:57

%%

if nargin < 7
    shift = [];
end

n  = length(wX);
cf = @(t) cfX0(sum(wX)*t) .* cfY0(sum(wY)*t);

for i = 1:n
    cf = @(t) cf(t) .* cfX{i}(wX(i)*t) .* cfY{i}(wY(i)*t);    
end

if ~isempty(shift)
    cf = @(t) cf(t) .* exp(1i*t*shift);
end
end
%% Function
function [xLow,xUpp] = CalInt(yLow,yUpp,xMin,xMax,probLow,probUpp,...
    order,LX,KY,cfX,cfX0,cfY,cfY0,pars,options)
% CalibrationInterval calculates the calibration interval [xmin,xmax] by
% inverting the calibration function. In fact, the values xmin, xmax are
% solutions to the following optimization problems:
%   xLow = argmin (yLow - quantile(x,probUpp))^2,
%   xUpp = argmin (yUpp - quantile(x,probLow))^2,
% where the function quantileY(x,prob) computes the quantile q of the
% conditional Y-distribution, given X = x, at the probability level
% specified by prob.  
%
% SYNTAX
%  [xLow,xUpp] = CalInt(yLow,yUpp,xMin,xMax,probLow,probUpp, ...
%                       order,LX,KY,cfX,cfX0,cfY,cfY0,pars,options)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Apr-2018 19:15:53


%% Algorithm 

% Check ymin and ymax. Valid values are from the support [yMin,yMax]
yMin = quantileY(xMin,probLow,order,LX,KY,cfX,cfX0,cfY,cfY0,pars,options);
yMax = quantileY(xMax,probUpp,order,LX,KY,cfX,cfX0,cfY,cfY0,pars,options);
if yLow < yMin || yUpp >yMax
    error('ymin or ymax out of range')
end

% set the optimization function 
xfun = @(x,y,pr) (y - ...
    quantileY(x,pr,order,LX,KY,cfX,cfX0,cfY,cfY0,pars,options))^2;

% Calculate the corresponding interval limits [xmin, xmax]
xLow = fminsearch(@(x)xfun(x,yLow,probUpp),(xMax-xMin)/2);
xUpp = fminsearch(@(x)xfun(x,yUpp,probLow),(xMax-xMin)/2);

end
%% Funtion quantileY
function q = quantileY(x,prob,order,LX,KY,cfX,cfX0,cfY,cfY0,pars,options)
% quantileY computes the quantile q of the conditional Y-distribution, for
% X = x, at the probability level specified by prob.  
%
% SYNTAX
% q = quantileY(x,prob,order,LX,KY,cfX,cfX0,cfY,cfY0,pars,options)
% 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Apr-2018 19:15:53

w = ones(1,order+1);
for i = 2:(order+1)
    w(i) = w(i-1).*x;
end
cf = CombinedCF(w*LX,w*KY,cfX,cfX0,cfY,cfY0,w*pars);
[~,~,~,q] = cf2DistGP(cf,[],prob,options);
end
%% Funtion CombinedCF
function cf = CombinedCF(wX,wY,cfX,cfX0,cfY,cfY0,shift)
% CombinedCF creates the combined 'calibration' characteristic function of
% the random variable W = sum_{i=1}^n (wX_i*(X_i+X_0) + wY_i*(Y_i+Y_0)), 
% defined by 
% cfW = @(t)(prod_{i=1}^n cfX{i}(wX(i)*t)).*cfX0((sum_{i=1}^I wX(i))*t) ...
%           (prod_{i=1}^n cfY{i}(wY(i)*t)).*cfY0((sum_{i=1}^I wY(i))*t), 
% where wX and wY are n-dimensional vectors of coefficients, and cfX =
% {@(t)cfX{1}(t), ..., @(t)cfX{n}(t)} and cfY = {@(t)cfY{1}(t), ...,
% @(t)cfY{n}(t)} are n-dimensional cell vectors of the characteristic
% functions of (X_1,...,X_n) and (Y_1,...,Y_n), respectively, and cfX0 and
% cfY0 are the characteristic functions of the 'common' random variables
% X_0 and Y_0, respectively.
%
% SYNTAX
% cf = CombinedCF(wX,wY,cfX,cfX0,cfY,cfY0)
%
% EXAMPLE
%  wX = [1 2 3 4 5]/15;
%  wY = [1 1 1 1 1]/5;
%  cfX  = {@(t)cfS_Gaussian(t), @(t)cfS_Rectangular(0.5*t), ...
%          @(t)cfS_Rectangular(1.2*t), @(t)cfS_Rectangular(2*t), ...
%          @(t)cfS_Arcsine(8*t)}';
%  cfX0 =  @(t)cfS_Rectangular(1.5*t);
%  cfY  =  {@(t)cfS_StudentT(0.5*t,10), @(t)cfS_Triangular(0.5*t), ...
%          @(t)cfS_Triangular(1.2*t), @(t)cfS_Triangular(2*t), ...
%          @(t)cfS_Triangular(2*t)}';
%  cfY0 =  @(t)cfS_Arcsine(3*t);
%  cf = CombinedCF(wX,wY,cfX,cfX0,cfY,cfY0)
%  figure
%  t = linspace(-5,5,201);
%  plot(t,real(cf(t)))
%  figure 
%  x      = linspace(-10,10,201);
%  prob   = [0.9 0.95 0.99];
%  result = cf2DistGP(cf,x,prob);
%  disp(result)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Apr-2018 13:19:57

%%

if nargin < 7
    shift = [];
end

n  = length(wX);
cf = @(t) cfX0(sum(wX)*t) .* cfY0(sum(wY)*t);

for i = 1:n
    cf = @(t) cf(t) .* cfX{i}(wX(i)*t) .* cfY{i}(wY(i)*t);    
end

if ~isempty(shift)
    cf = @(t) cf(t) .* exp(1i*t*shift);
end
end