function cf = cfTest_Bartlett(t,df)
%% cfTest_Bartlett 
%  Characteristic function of the exact null distribution of the BARTLETT's
%  test statistic for testing hypothesis about homogeneity of variances
%  of k normal populations, specified by the vector of the degrees of
%  freedom df = [df_1,...,df_k].   
%  
%  The characteristic function of the the exact null distribution of the
%  BARTLETT's test statistic is given by 
%  cf(t) = exp(1i*t*C/B).*
%          k.^(1i*t*D/B).*(gamma(D/2)/gamma(D/2-1i*t*D/B)).* 
%          prod(gamma(df/2-1i*t*df/B)./gamma(df/2))
%  where B = 1+1/(3*(k-1))*(sum(1./df)-1/D) and C = D*log(k*prod(W.^W));
%  with W = df/D, and D = sum(df).  
% 
% SYNTAX
%   cf = cfTest_Bartlett(t,df)
%
% INPUTS
%   t    - vector or array of real values, where the CF is evaluated.
%   df   - k-dimensional vector of degrees of freedom df = [df_1,...,df_k]. 
%
% WIKIPEDIA
%   https://en.wikipedia.org/wiki/Bartlett%27s_test
%
% EXAMPLE 1:
% % CF of Bartlett distribution with the specified degrees of freedom
%   df    = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3];
%   t     = linspace(-1,1,201);
%   cf    = cfTest_Bartlett(t,df);
%   figure; plot(t,real(cf),t,imag(cf)); grid on;
%   title('CF of the Bartlett distribution with specified dfs')
%
% EXAMPLE 2:
% % PDF/CDF of Bartlett distribution with the specified degrees of freedom
%   df    = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3];
%   cf    = @(t) cfTest_Bartlett(t,df);
%   x     = linspace(0,40,200)';
%   prob  = [0.9 0.95 0.99];
%   clear options
%   options.xMin = 0;
%   result = cf2DistGP(cf,x,prob,options);
%   disp(result)
%
% EXAMPLE 3:
% % Quantiles of Bartlett distribution computed by the algorithm cf2DistGPA
%   k     = 5;
%   df    = [1 2 3 4 5 6 7 8 9 19 29 49 99];
%   prob  = [0.9 0.95 0.99];
%   Table = zeros(13,3);
%   clear options
%   options.isPlot = false;
%   for i = 1:13
%       disp(['df = ',num2str(df(i))])
%       nu = df(i)*ones(k,1);
%       cf = @(t) cfTest_Bartlett(t,nu);
%       [~,~,~,qf] = cf2DistGPA(cf,[],prob,options);
%       disp(qf(:)')
%       Table(i,:) = qf;
%   end
%   disp(Table)
%
% EXAMPLE 4:
% % Quantiles of Bartlett distribution computed by the algorithm cf2QF
%   k     = 5;
%   df    = [1 2 3 4 5 6 7 8 9 19 29 49 99];
%   prob  = [0.9 0.95 0.99];
%   Table = zeros(13,3);
%   clear options
%   options.crit = 1e-10;
%   options.maxiter = 10;
%   for i = 1:13
%       disp(['df = ',num2str(df(i))])
%       nu = df(i)*ones(k,1);
%       cf = @(t) cfTest_Bartlett(t,nu);
%       qf = cf2QF(cf,prob,options);
%       disp(qf(:)')
%       Table(i,:) = qf;
%   end
%   disp(Table)
%
% REFERENCES
%  Glaser, R. E. (1976a). The ratio of the geometric mean to the arithmetic
%  mean for a random sample from a gamma distribution. Journal of the
%  American Statistical Association, 71(354), 480-487.   
%
%  Glaser, R. E. (1976b). Exact critical values for Bartlett's test for
%  homogeneity of variances. Journal of the American Statistical
%  Association, 71(354), 488-490.  
%
%  Chao, M. T., & Glaser, R. E. (1978). The exact distribution of
%  Bartlett's test statistic for homogeneity of variances with unequal
%  sample sizes. Journal of the American Statistical Association, 73(362),
%  422-426.   

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 09-Dec-2018 09:42:32

%% ALGORITHM
% cf = cfTest_Bartlett(t,df)

%% CHECK THE INPUT PARAMETERS
narginchk(2, 2);

%% Characteristic function of the Bartlett's null distribution
szt = size(t);
t   = t(:);
df  = df(:);
k   = length(df);
D   = sum(df);
W   = df/D;
C   = D*log(k*prod(W.^W));
B   = 1 + (1/(3*(k-1))) * (sum(1./df)-1/D);

cf  = GammaLog(D/2);
cf  = cf + (1i*t*C/B) - log(k)*(1i*t*D/B) - GammaLog(D/2-1i*t*D/B);
for j = 1:k
    cf = cf + GammaLog(df(j)/2-1i*t*df(j)/B) - GammaLog(df(j)/2);
end
cf  = reshape(exp(cf),szt);
cf(t==0) = 1;

end