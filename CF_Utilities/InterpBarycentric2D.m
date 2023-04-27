function [funNew,xyNew] = ...
    InterpBarycentric2D(xOld,yOld,funOld,xyNew,options)
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
%   funNew = InterpBarycentric2D(xOld,yOld,funOld,xyNew)
%   [funNew,xyNew] = InterpBarycentric2D(xOld,yOld,funOld,xyNew,options)
%
% EXAMPLE1 (PDF of bivariate distribution evaluated at random poins)
%  cf = @(t) exp(-(0.9*t(:,1).^2 + 0.3*t(:,2).^2 +2*0.4*t(:,1).*t(:,2))/2);
%  clear options;
%  options.isPlot = false;
%  result = cf2Dist2D(cf,[],options);
%  xyNew = randn(100,2);
%  pdfNew = InterpBarycentric2D(result.x1,result.x2,result.Zpdf,xyNew);
%  disp([xyNew pdfNew])
%
% EXAMPLE2 (Barycentric interpolation of bivariate PDF specified by its CF)
%  cf = @(t) exp(-(0.9*t(:,1).^2 + 0.3*t(:,2).^2 +2*0.4*t(:,1).*t(:,2))/2);
%  clear options;
%  options.isInterp = true;
%  options.isPlot = false;
%  result = cf2Dist2D(cf,[],options);
%  pdfOld = result.Zpdf;
%  xOld   = result.x1;
%  yOld   = result.x2;
%  xNew   = linspace(min(xOld),max(xOld),501)';
%  yNew   = linspace(min(yOld),max(yOld),501)';
%  xyNew  = {xNew,yNew};
%  pdfNew = InterpBarycentric2D(xOld,yOld,pdfOld,xyNew);
%  [X,Y]  = meshgrid(xNew,yNew);
%  mesh(X,Y,pdfNew)
%
% EXAMPLE3 (Barycentric interpolation of bivariate PDF/CDF specified by CF)
%  cf     = @(t) cf2D_Logistic(t);
%  clear options;
%  options.isInterp = true;
%  options.isPlot = false;
%  result = cf2Dist2D(cf,[],options);
%  cdfOld = result.Zcdf;
%  pdfOld = result.Zpdf;
%  xOld   = result.x1;
%  yOld   = result.x2;
%  xNew   = linspace(min(xOld),max(xOld),501)';
%  yNew   = linspace(min(yOld),max(yOld),501)';
%  xyNew  = {xNew,yNew};
%  cdfNew = InterpBarycentric2D(xOld,yOld,cdfOld,xyNew);
%  pdfNew = InterpBarycentric2D(xOld,yOld,pdfOld,xyNew);
%  [X,Y]  = meshgrid(xNew,yNew);
%  figure;mesh(X,Y,cdfNew)
%  figure;mesh(X,Y,pdfNew)
%  figure
%  [~,c1] =contour(X,Y,cdfNew,'ShowText','on');hold on;
%  [~,c2] = contour(X,Y,pdfNew,'ShowText','on');grid;
%  c1.LevelList = [0.001 0.01 0.05 0.1000 0.2000 0.3000 0.4000 0.5000 ...
%    0.6000 0.7000 0.8000 0.9000 0.95 0.99];
%  c2.LevelList = [0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07];
%  hold off
%
% EXAMPLE4 (Barycentric interpolation of bivariate PDF/CDF specified by CF)
%  mu1    = [0 2];
%  beta1  = [1 2];
%  cf1    = @(t) cf2D_Logistic(t,mu1,beta1);
%  mu2    = [2 1];
%  beta2  = [2 1];
%  cf2    = @(t) cf2D_Logistic(t,mu2,beta2);
%  cf     = @(t) 0.25*cf1(t) + 0.75*cf2(t);
%  clear options;
%  options.isPlot   = false;
%  options.isInterp = true;
%  options.chebyPts = 51;
%  result = cf2Dist2D(cf,[],options);
%  cdfOld = result.Zcdf;
%  pdfOld = result.Zpdf;
%  xOld   = result.x1;
%  yOld   = result.x2;
%  xNew   = linspace(min(xOld),max(xOld),501)';
%  yNew   = linspace(min(yOld),max(yOld),501)';
%  xyNew  = {xNew,yNew};
%  cdfNew = InterpBarycentric2D(xOld,yOld,cdfOld,xyNew);
%  pdfNew = InterpBarycentric2D(xOld,yOld,pdfOld,xyNew);
%  [X,Y]  = meshgrid(xNew,yNew);
%  figure;mesh(X,Y,cdfNew)
%  figure;mesh(X,Y,pdfNew)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 26-Apr-2023 14:22:19

%% FUNCTION
%  [funNew,xyNew] = InterpBarycentric2D(xOld,yOld,funOld,xyNew,options)
%% CHECK THE INPUT PARAMETERS

narginchk(3,5);
if nargin < 5, options = []; end
if nargin < 4, xyNew   = []; end

if iscell(xyNew)
    xNew = xyNew{1};
    yNew = xyNew{2};
    Nxy  = [];
    Nx   = length(xNew);
    Ny   = length(yNew);
else
    xNew = [];
    yNew = [];
    Nx   = [];
    Ny   = [];
    Nxy  = size(xyNew,1);
end

if isempty(xyNew)
    xNew = xOld;
    yNew = yOld;    
    Nxy  = [];
    Nx   = length(xNew);
    Ny   = length(yNew);
end

if ~isfield(options, 'isChebPts')
    options.isChebPts = true;
end

if ~isfield(options, 'wx')
    options.wx = [];
end

if ~isfield(options, 'wy')
    options.wy = [];
end

%% ALGORITHM
%tol    = options.SVDtol;
wx     = options.wx;
wy     = options.wy;
nx     = length(xOld(:));
ny     = length(yOld(:));


if nx < 2 || ny < 2
    error('small number of nodes');
end

% SVD decomposition of funOld
[U,S,V] = svd(funOld,'econ');
n  = size(S,1);

% funNew using SVD decomposition
if ~isempty(Nxy)
    xNew = xyNew(:,1);
    yNew = xyNew(:,2);
    u = zeros(Nxy,n);
    v = zeros(Nxy,n);
    for j = 1:n
        u(:,j) = InterpBarycentric(yOld,U(:,j),yNew,wy,options) * S(j,j);
        v(:,j) = InterpBarycentric(xOld,V(:,j),xNew,wx,options);
    end
    funNew = sum(reshape(u(:).*v(:),Nxy,n),2);
else
    u = zeros(Ny,n);
    v = zeros(Nx,n);
    for j = 1:n
        u(:,j) = InterpBarycentric(yOld,U(:,j),yNew,wy,options) * S(j,j);
        v(:,j) = InterpBarycentric(xOld,V(:,j),xNew,wx,options);
    end
    funNew = u * v';
end

end