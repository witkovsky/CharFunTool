function [funNew,xNew,yNew,xyNew] = ...
    InterpBarycentric2D(xOld,yOld,funOld,xNew,yNew,xyNew,wx,wy,options)
%   InterpBarycentric2D is a 2-dimensional barycentric interpolation.
%   Evaluates (interpolates) 2D-function, based on the n x m dimensional
%   matrix of funOld values, evaluated at all combinations (meshgrid) of
%   xOld and yOld, at given new points xNew and yNew.
%
%   The algorithm is using the barycentric interpolation of the columns of
%   U and V  matrices from the SVD decomposition of funOld, i.e.  [U,S,V] =
%   svd(funOld). The result is N x M dimensional matrix of interpolated
%   values funNew, evalueted at all combinations xNew and yNew, or a vector
%   of funNew values of length as the number of pairs in the specified
%   combinations xyNew. 
%  
% SYNTAX:
%   funNew = InterpBarycentric2D(xOld,yOld,funOld,xNew,yNew)
%   [funNew,xNew,yNew,xyNew] = ...
%    InterpBarycentric2D(xOld,yOld,funOld,xNew,yNew,xyNew,wx,wy,options)
%
% EXAMPLE1 (Barycentric interpolation of bivariate distribution)
%  cf = @(t) exp(-(0.9*t(:,1).^2 + 0.3*t(:,2).^2 +2*0.4*t(:,1).*t(:,2))/2);
%  options.isInterp = true;
%  result = cf2Dist2D(cf,[],options)
%  pdfOld = result.Zpdf;
%  xOld   = result.x1;
%  yOld   = result.x2;
%  xNew   = linspace(min(xOld),max(xOld),201)';
%  yNew   = linspace(min(yOld),max(yOld),101)';
%  pdfNew = InterpBarycentric2D(xOld,yOld,pdfOld,xNew,yNew);
%  [X,Y]  = meshgrid(xNew,yNew);
%  mesh(X,Y,pdfNew)
%
% EXAMPLE2 (Barycentric interpolation of bivariate distribution)
%  mu1    = [0 2];
%  beta1  = [1 2];
%  cf1    = @(t) cf2D_Logistic(t,mu1,beta1);
%  mu2    = [2 1];
%  beta2  = [2 1];
%  cf2    = @(t) cf2D_Logistic(t,mu2,beta2);
%  cf     = @(t) 0.25*cf1(t) + 0.75*cf2(t);
%  clear options;
%  options.xN = 51;
%  options.isChebyshev = false;
%  options.SVDtol = 1e-10;
%  result = cf2Dist2D(cf,[],options)
%  pdfOld = result.Zpdf;
%  xOld   = result.x1;
%  yOld   = result.x2;
%  xNew   = linspace(min(xOld),max(xOld),201)';
%  yNew   = linspace(min(yOld),max(yOld),151)';
%  pdfNew = InterpBarycentric2D(xOld,yOld,pdfOld,xNew,yNew);
%  [X,Y]  = meshgrid(xNew,yNew);
%  mesh(X,Y,pdfNew)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 02-May-2021 21:27:52

%% FUNCTION
%  [funNew,xNew,yNew,xyNew] = ...
%   InterpBarycentric2D(xOld,yOld,funOld,xNew,yNew,xyNew,wx,wy,options)

%% CHECK THE INPUT PARAMETERS

narginchk(3, 9);
if nargin < 9, options = []; end
if nargin < 8, wy      = []; end
if nargin < 7, wx      = []; end
if nargin < 6, xyNew   = []; end
if nargin < 5, yNew    = []; end
if nargin < 4, xNew    = []; end

if isempty(xNew) && isempty(xyNew)
    xNew = xOld;
end

if isempty(yNew) && isempty(xyNew)
    yNew = yOld;
end

if ~isempty(xyNew)
    xNew = [];
    yNew = [];
    Nx   = [];
    Ny   = [];
    Nxy  = size(xyNew,1);
else
    Nxy  = [];
    Nx   = length(xNew);
    Ny   = length(yNew);
end

if ~isfield(options, 'isChebPts')
    options.isChebPts = true;
end

if ~isfield(options, 'SVDtol')
    options.SVDtol = 1e-14;
end

%% ALGORITHM
tol    = options.SVDtol;
nx     = length(xOld(:));
ny     = length(yOld(:));

if nx < 2 || ny < 2
        error('small number of nodes');
end

% SVD decomposition of funOld
[U,S,V] = svd(funOld);
id = diag(S) > tol;
n  = sum(id);
U = U(:,id);
S = S(id,id);
V = V(:,id);

% funNew using SVD decomposition
if ~isempty(Nxy)
    funNew = zeros(Nxy,1);
    u = zeros(1,n);
    v = zeros(1,n);
    for i = 1:Nxy
        for j = 1:n
            u(j) = InterpBarycentric(xOld,U(:,j),xyNew(i,1),wx,options);
            v(j) = InterpBarycentric(yOld,V(:,j),xyNew(i,2),wy,options);
        end
        funNew(i) = u * S * v';
    end
else
    u = zeros(Ny,n);
    v = zeros(Nx,n);
    for j = 1:n
        u(:,j) = InterpBarycentric(yOld,U(:,j),yNew,wy,options);
        v(:,j) = InterpBarycentric(xOld,V(:,j),xNew,wx,options);
    end
    funNew = u * S * v';
end

end