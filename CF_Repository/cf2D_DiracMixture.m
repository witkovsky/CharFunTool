function cf = cf2D_DiracMixture(t,data,weight)
%% cf2D_DiracMixture  
%  Characteristic function of the weighted Bivariate Mixture Distribution.
%  which is created by a mixture of independent bivariate DIRAC random
%  vectors D_1 = (D_11,D_12),...,D_N = (D_N1,D_N2), concentrated at N fixed
%  2D-vectors given by the N x 2 matrix, data = [d_11,d_12;...;d_N1,d_N2]. 
%  
%  For given data and t = [t1,t2] the characteristic function is given as
%  cf(t) = weight_1*exp(1i*(d_11*t1+d_12*t2)) +...+
%  weight_N*exp(1i*(d_N1*t1+d_N2*t2)), where exp(1i*(d_1*t1+d_2*t2))
%  represents the characteristic function of the bivariate DIRAC RV
%  concentrated at the vector d = [d_1,d_2].
%
%  
% SYNTAX:
%  cf = cf2D_DiracMixture(t,data,weight)
%
% INPUTS:
%  t      - (n x 2)-matrix t = [t1 t2] or a cell t = {t1 t2} of real
%           values, where the CF is evaluated. If t is cell of two vectors,
%           then cf is a (n1 x n2)-matrix where n1 = length(t1), n2 =
%           length(t2). 
%  data   - (n x 2)-matrix data = [d1 d2] or a cell data = {d1 d2} of real
%           values. If data is cell of two vectors d1 and d2, then the new
%           data matrix is created from the meshgrid, i.e. [d2,d1] =
%           meshgrid(data{2},data{1}) and hence new data = [d1(:),d2(:)]. 
%  weight - vector of weights of the distribution mixture. If empty,
%           default value is weight = 1/length(d1). 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
%
% EXAMPLE 1:
% % CF of the weighted Bivariate Mixture Distribution
%   N       = 100;
%   data    = rand(N,2);
%   weight  = 1/N;
%   cf      = @(t) cf2D_DiracMixture(t,data,weight);
%   tt      = linspace(-4,4,51);
%   [t1 t2] = meshgrid(tt,tt);
%   t = [t1(:) t2(:)];
%   figure; 
%   contour(t1,t2,real(reshape(cf(t),51,51)),'ShowText','on')
%   title('Contour Plot of the Bivariate Mixture Distribution')
%
% EXAMPLE 2:
% % PDF/CDF of weighted Bivariate Mixture Distribution 
% % smoothed by Gaussian kernel 
%   N       = 100;
%   data    = rand(N,2);
%   weight  = 1/N;
%   bandwidth = 0.1;
%   cf = @(t) cf2D_DiracMixture(t,data,weight)  .* cf2D_Normal(bandwidth*t);
%   clear options;
%   options.isInterp = true;
%   result = cf2Dist2D(cf,[],options)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 07-Jan-2023 21:43:54

%% ALGORITHM
% cf = cf2D_DiracMixture(t,data,weight)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 3, weight = []; end
if nargin < 2, data = []; end

if isempty(data), data = [0,0]; end

%% CHECK the arguments t

if iscell(t)
    sz = [length(t{1}),length(t{2})];    
    [t2,t1] = meshgrid(t{2},t{1});
    t  = [t1(:),t2(:)];
else
    szt = size(t);
    sztMin = min(szt(1),szt(2));
    sztMax = max(szt(1),szt(2));    
    switch sztMin
        case 1 % if t is N-vector then it is assumed that t1 = t and t2 = t
               % and we create new t = [t1 t2]
            if sztMax > 2
                sz = size(t);
                t = t(:);
                t = [t t];
            elseif sztMax == 1
                sz = [1,1];
                t = [t t];
            else
                sz = [1,1];
                t = t(:)';
            end
        case 2 % t is 2xN or Nx2 matrix
            if sztMax > 2
                if szt(1) == 2 % if t is 2xN matrix transpose it to the Nx2 matrix
                    t = t';
                    sz = [1,szt(2)];
                else
                    sz = [szt(1),1];
                end
            end
        otherwise
            error(message('InputSizeMismatch'));
    end
end

t1 = t(:,1);
t2 = t(:,2);

%% CHECK the data

if iscell(data) 
    [d2,d1] = meshgrid(data{2},data{1});
    data  = [d1(:),d2(:)];
else
    szd = size(data);
    szdMin = min(szd(1),szd(2));
    szdMax = max(szd(1),szd(2));    
    switch szdMin
        case 1 % if data is N-vector then it is assumed that d1 = data and d2 = data
               % and we create new data = [d1 d2]
            if szdMax > 2
                data = data(:);
                data = [data data];
            elseif szdMax == 1
                data = [data data];
            else
                data = data(:)';
            end
        case 2 % data is 2xN or Nx2 matrix
            if szdMax > 2
                if szd(1) == 2 % if data is 2xN matrix transpose it to the Nx2 matrix
                    data = data';
                end
            end
        otherwise
            error(message('InputSizeMismatch'));
    end
end

d1 = data(:,1);
d2 = data(:,2);

N  = length(d1);
if isempty(weight), weight = 1/N; end

%% Equal size of the parameters   

[err,d1,weight] = distchck(2,d1,weight);
if err > 0
    error(message('InputSizeMismatch'));
end

%% Characteristic function

cf = exp(1i * ((t1 * d1') + (t2 * d2')))* weight;
cf = reshape(cf,sz);

end