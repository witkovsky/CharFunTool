function [cdf,result] = cf2CDF_BTAV(cf,x,options)
% cf2CDF_BTAV Calculates CDF of a NON-NEGATIVE distribution specified by
% its CF at values x by using numerical inversion  based on the
% Bromwich-Talbot-Abate-Valko method. The possible quadrature methods
% include the simple trapezoidal rule (the default method) or the MATLAB
% integral algorithm (adaptive Gauss-Kronrod quadrature). 
% 
% SYNTAX:
% [cdf,result] = cf2CDF_BTAV(cf,x,options)
%
% INPUT:
%  cf       - function handle of the characteristic function
%  x        - vector of values where the CDF/PDF is computed
%  options  - structure with the following parameters:
%             options.quadrature       % quadrature method. Default method
%                  = 'trapezoidal'     % is quadrature = 'trapezoidal'.
%                                      % Alternatively use quadrature =
%                                      % 'matlab' (for the MATLAB built in
%                                      % adaptive Gauss-Kronrod quadrature) 
%                                      % Gauss-Kronrod integral. 
%             options.tol = 1e-12      % absolute tolerance for the MATLAB
%                                      % Gauss-Kronrod integral. 
%             options.nTerms = 100     % number of terms used in the
%                                      % trapezoidal quadrature
%                                      % estimated from the specified CF
%             options.Mpar_BTAV = 10   % parameter M for the deformed
%                                      % Bromwich curve.
%
% OUTPUT:
%  cdf      - vector of the CDF values evaluated at x.
%  result   - structure with further details.
%
% EXAMPLE 1
% % CDF of the Chi-squared distribution with DF = 1
%  df  = 1;
%  cf  = @(t) (1 - 2i*t).^(-df/2);
%  x   = linspace(0,10)';
%  cdf = cf2CDF_BTAV(cf,x);
%  figure;plot(x,cdf,'.-')
%
% EXAMPLE 2
% % CDF of the exact null distribution for testing Compound Symmetry
%  n  = 10; 
%  p  = 5;
%  type = 'modified';
%  cf = @(t) cfTest_CompoundSymmetry(t,n,p,type);
%  x  = linspace(0,5)';
%  clear options
%  options.quadrature = 'trapezoidal';
%  options.isPlot = true;
%  [cdf,result] = cf2CDF_BTAV(cf,x,options);
%  disp(result)
% 
% EXAMPLE 3
% % CDF of the exact null distribution for testing Equality Covariances
%  n  = 10; 
%  p  = 5;
%  q  = 3;
%  type = 'modified';
%  cf = @(t) cfTest_EqualityCovariances(t,n,p,q,type);
%  x  = linspace(0,10)';
%  clear options
%  options.quadrature = 'matlab';
%  options.isPlot = true;
%  [cdf,result] = cf2CDF_BTAV(cf,x,options);
%  disp(result)
%
% REFERENCES:
% [1] Talbot, A., 1979. The accurate numerical inversion of Laplace
%     transforms. IMA Journal of Applied Mathematics, 23(1), pp.97-120.
% [2] Abate, J. and Valkó, P.P., 2004. Multi-precision Laplace transform
%     inversion. International Journal for Numerical Methods in
%     Engineering, 60(5), pp.979-993.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 11-Jan-2019 15:55:06

%% CHECK/SET THE INPUT PARAMETERS
StartTime = cputime;
narginchk(2, 3);
if nargin < 3, options = []; end

if ~isfield(options, 'quadrature')
    options.quadrature = 'trapezoidal';
end

if ~isfield(options, 'tol')
    options.tol =1e-12;
end

if ~isfield(options, 'nTerms')
    options.nTerms = 100;
end

if ~isfield(options, 'Mpar_BTAV')
    options.Mpar_BTAV = 10;
end

if ~isfield(options, 'isPlot')
    options.isPlot = false;
end

%% ALGORITHM
szx = size(x);
x   = x(:);

quadrature = options.quadrature;
switch lower(quadrature)
    case 'matlab'
        tol = options.tol;
        M   = options.Mpar_BTAV;
        cdfFun = @(t) real(IntegrandFun_BTAV(t,x,cf,'cdf',M));
        cdf = integral(cdfFun,0,pi, ...
            'ArrayValued',true,'RelTol',0,'AbsTol',tol)/pi;
    case 'trapezoidal'
        n   = options.nTerms;
        M   = options.Mpar_BTAV;
        t   = linspace(0,pi,n+1)';
        cdfFun = IntegrandFun_BTAV(t,x,cf,'cdf',M);
        cdf = (cdfFun(1,:)/2 + nansum(real(cdfFun(2:n,:))))/n;
end

cdf(isnan(cdf)) = 0;
cdf    = reshape(cdf,szx);
tictoc = cputime - StartTime;

%% RESULTS
if nargout > 1
    result.Description     = 'CDF from the characteristic function CF';
    result.inversionMethod = 'Bromwich-Talbot-Abate-Valkó';
    result.quadratureMethod = options.quadrature;
    result.cdf = cdf;
    result.x   = x;
    result.cf  = cf;
    result.cdfFun  = cdfFun;
    result.options = options;
    result.tictoc  = tictoc; 
end

%% PLOT the PDF / CDF 
if length(x)==1
    options.isPlot = false;
end
if options.isPlot
    % CDF
    figure
    plot(x,cdf,'.-')
    grid
    title('CDF Specified by the CF')
    xlabel('x')
    ylabel('cdf')        
end
end