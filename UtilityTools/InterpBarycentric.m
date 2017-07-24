function [funNew,xNew] = InterpBarycentric(x,fun,xNew,options)
%   InterpBarycentric Evaluates (interpolates) function values f_new at
%   given points x_new by barycentric interpolation from function values f
%   given at chebpoints x. 
%   For more details see: 
%   https://en.wikipedia.org/wiki/Barycentric_coordinate_system
%  
% SYNTAX:
%   funNew = InterpBarycentric(x,fun,xNew)
%   [funNew,xNew] = InterpBarycentric(x,fun,xNew,options)
%
% EXAMPLE1 (Barycentric interpolant of the Sine function on (-pi,pi))
%   x    = ChebPoints(32,[-pi,pi])';
%   f    = sin(x);
%   xNew = linspace(-pi,pi,201)';
%   fNew = InterpBarycentric(x,f,xNew);
%   plot(x,f,xNew,fNew,'.')
%   disp([xNew fNew sin(xNew)])
%
% EXAMPLE2 (Barycentric interpolant of the Normal CDF)
%   x    = ChebPoints(21,[-8,8])';
%   f    = normcdf(x);
%   xNew = linspace(-5,5,201)';
%   fNew = InterpBarycentric(x,f,xNew);
%   plot(x,f,xNew,fNew,'.')
%   disp([xNew fNew normcdf(xNew)])
%
% EXAMPLE3 (Barycentric interpolant of the Normal quantile function)
%   x    = ChebPoints(2^7,[-8,8])';
%   cdf  = normcdf(x);
%   prob = linspace(1e-4,1-1e-4,101)';
%   qf   = InterpBarycentric(cdf,x,prob);
%   plot(cdf,x,prob,qf,'.')
%   disp([prob qf norminv(prob)])
%
% EXAMPLE4 (Barycentric interpolant of the Compound distribution)
%   cfX  = @(t) cfX_Exponential(t,5);
%   cf   = @(t) cfN_Poisson(t,10,cfX);
%   x    = ChebPoints(2^9,[0,10]);
%   prob = [0.9 0.95 0.99];
%   clear options
%   options.isCompound = 1;
%   result = cf2DistGP(cf,x,prob,options)
%   CDF  = @(x) InterpBarycentric(result.x,result.cdf,x);
%   PDF  = @(x) InterpBarycentric(result.x,result.pdf,x);
%   QF   = @(p) InterpBarycentric(result.cdf,result.x,p);
%   prob = linspace(1e-4,1-1e-4,201);
%   figure;plot(result.cdf,result.x,prob,QF(prob),'.');

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 10:06:48

%% FUNCTION
%  [funNew,xNew] = InterpBarycentric(x,fun,xNew,options);

%% CHECK THE INPUT PARAMETERS

narginchk(2, 4);
if nargin < 3, xNew = []; end
if nargin < 4, options = []; end

if isempty(xNew)
    xNew = linspace(min(x),max(x))';
end

if ~isfield(options, 'isChebPts')
    options.isChebPts = true;
end

%% ALGORITHM
x      = x(:);
fun    = fun(:);
[n,m]  = size(xNew);
xNew   = xNew(:);
nx     = length(x);
nxNew  = length(xNew);
funNew = zeros(nxNew,1);

w     = (-1).^(0:nx-1).';
w(1)  = w(1)/2;
w(nx) = w(nx)/2;
for i = 1:nxNew
    A = 0;
    B = 0;
    for j = 1:nx
        if xNew(i) == x(j)
            exactVal  = true;
            funNew(i) = fun(j);
        else
            exactVal = false;
            weight   = w(j) ./ (xNew(i)-x(j));
            A = A + weight * fun(j);
            B = B + weight;
        end
        if exactVal == true
            break;
        else
            funNew(i) = A / B;
        end
    end
end

xNew   = reshape(xNew,n,m);
funNew = reshape(funNew,n,m);
end