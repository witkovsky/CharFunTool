function cf = cfND_DiracMixture(t,data,weight)
%% cfND_DiracMixture  
%  Characteristic function of the weighted N-VARIATE MIXTURE DISTRIBUTION
%  which is created by a mixture of n independent N-dimensional DIRAC
%  random vectors D_1 = (D_11,...,D_1N),...,D_n = (D_n1,...,D_nN),
%  concentrated at n fixed N-dimensional vectors given by the (n x N)
%  matrix, data = [d_11,...,d_1N;...;d_n1,...,d_nN].
%  
%  For given data and t = [t1,...,tN] the characteristic function is given
%  as cf(t) = weight_1*exp(1i*(d_11*t1 +...+ d_1N*tN)) +...+
%  weight_n*exp(1i*(d_n1*t1 +...+ d_nN*tN)), where exp(1i*(d_1*t1 +...+
%  d_N*tN)) represents the characteristic function of the N-dimensional
%  DIRAC RV concentrated at the vector d = [d_1,...,d_N].
%
%  
% SYNTAX:
%  cf = cfND_DiracMixture(t,data,weight)
%
% INPUTS:
%  t      - (M x N)-matrix t = [t1,...,tN] of real values, where the CF is
%           evaluated.  
%  data   - (n x N)-matrix data = [d1,...,dN] of real values.  
%  weight - vector of weights of the distribution mixture. If empty,
%           default value is weight = 1/n. 
%
% WIKIPEDIA: 
%  https://en.wikipedia.org/wiki/Empirical_distribution_function.
%
% EXAMPLE 1:
% % CF of the weighted n-dimensional Mixture Distribution
%   n       = 100;
%   N       = 3;
%   data    = rand(n,N);
%   weight  = 1/n;
%   cf      = @(t) cfND_DiracMixture(t,data,weight);
%   tt      = linspace(-4,4,21);
%   [t1 t2 t3] = meshgrid(tt);
%   t       = [t1(:) t2(:) t3(:)];
%   disp([t cf(t)])
%
% EXAMPLE 2:
% % PDF/CDF of weighted Bivariate Mixture Distribution 
% % smoothed by Gaussian kernel 
%   n       = 100;
%   N       = 2;
%   data    = rand(n,N);
%   weight  = 1/n;
%   bandwidth = 0.1;
%   cf = @(t) cfND_DiracMixture(t,data,weight)  .* cf2D_Normal(bandwidth*t);
%   clear options;
%   options.isInterp = true;
%   result = cf2Dist2D(cf,[],options)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 01-May-2024 16:29:07'

%% ALGORITHM
% cf = cfND_DiracMixture(t,data,weight)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 3);
if nargin < 3, weight = []; end

%% CHECK the data

[n,N] = size(data);

if isempty(weight), weight = 1/n; end
d1 = data(:,1);

%% CHECK the arguments t

szt = size(t);
sztMin = min(szt(1),szt(2));
sztMax = max(szt(1),szt(2));
switch sztMin
    case 1 % if t is M-vector then it is assumed that t1 = t, ... tN = t
        % and we create new t = [t1,...,tN]
        if sztMax > N
            sz = size(t);
            t = t(:);
            t = t * ones(1,N);
        end
        if sztMax == N
            t = t(:)';
            sz = [1 1];
        end
    case N % if t is Nxn matrix transpose it to the nxN matrix
        if sztMax > N &&  szt(1) == N
            t = t';
            sz = [1,szt(2)];
        else
            sz = [szt(1),1];
        end
    otherwise
        error('InputSizeMismatch');
end

t1 = t(:,1);

%% Equal size of the parameters   


[err,d1,weight] = distchck(2,d1,weight);
if err > 0
    error('InputSizeMismatch');
end

%% Characteristic function

arg = t1*d1';
for i = 2:N
    arg = arg + t(:,i)*data(:,i)';
end
cf = exp(1i * arg ) * weight;
cf = reshape(cf,sz);

end