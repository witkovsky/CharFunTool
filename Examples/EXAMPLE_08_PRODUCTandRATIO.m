%% CharFunTool / EXAMPLES
%  CHARACTERISTIC FUNCTIONS OF THE PRODUCT & RATIO OF TWO INDEPENDENT RVs
%
%  (c) Viktor Witkovsky, witkovsky@gmail.com
%  Version: 13-Feb-2018 12:37:57

clear
close all

%% *PRODUCT OF TWO INDEPENDENT RVs*
% |cf_Prod| 
%  Characteristic function of a product of two independent (continuous)
%  random variables, i.e. Y = X1 * X2. The distribution of the first random
%  variable X1 is given by its characteristic function (CF), say cf_X1(t),
%  and the distribution of the second random variable is given by its
%  probability density function (PDF), say pdf_X2(x).
% 
%  Then, the characteristic function of the product of two independent
%  random variables, say cf_Y(t), is given by
%   cf_Y(t) = integral(@(x) cf_X1(x*t) .* pdf_X2(x),a,b),
%  where a and b are limits of integration (interval representing the
%  support of the distribution of X2). 
%
% SYNTAX:
%  cf = cf_Prod(t,cfX1,pdfX2,a,b)


%% I.1.a EXAMPLE: CF of a product of two independent RVs  
%  % with central symmetric triangular distribution on [-1,1] and
%  % rectangular distributions on [-1,1]

t = linspace(-100,100,2^10+1);
cfX1 = @(t) cfS_Triangular(t);
cf   = cf_Prod(t,cfX1);
figure; plot(t,real(cf),t,imag(cf));grid on
title('CF of a product of triangular and rectangular RVs ')

%% I.1.b EXAMPLE: CF of a product of two independent RVs 
%  % with chi-squared (df = 5) and symmetric rectangular distributions on
%  % [-1,1].

t = linspace(-5,5,2^10+1);
cfX1 = @(t) cfX_ChiSquare(t,5);
cf   = cf_Prod(t,cfX1);
figure; plot(t,real(cf),t,imag(cf));grid on
title('CF of a product of chi-squared and rectangular RVs ')

%% I.1.c EXAMPLE: CF of a product of two independent RVs 
%  % with chi-squared (df = 5) distribution and the standard normal
%  % (Gaussian) distribution. 

t = linspace(-5,5,2^10+1);
cfX1  = @(t) cfX_ChiSquare(t,5);
pdfX2 = @(x) normpdf(x);
a     = -8;
b     =  8;
cf    = cf_Prod(t,cfX1,pdfX2,a,b);
figure; plot(t,real(cf),t,imag(cf));grid on
title('CF of a product of chi-squared and normal RVs ')

%% I.1.d EXAMPLE: PDF/CDF of a product of two independent RVs 
%  % with central symmetric triangular distribution on [-1,1] and
%  % rectangular distributions on [-1,1]

cfX1  = @(t) cfS_Triangular(t);
a     = -1;
b     = 1;
pdfX2 = @(x) 1/(b-a);
tol   = 1e-4;
cf    = @(t) cf_Prod(t,cfX1,pdfX2,a,b,tol);
x     = linspace(-1,1,501);
prob  = [0.9 0.95 0.99];

clear options;
options.N    = 2^9;
options.xMin = -1;
options.xMax = 1;
figure;
result = cf2DistGP(cf,x,prob,options);

disp(result)

%% I.1.e EXAMPLE: PDF/CDF of a product of two independent RVs 
%  % with chi-squared (df = 5) distribution and the standard normal
%  % (Gaussian) distribution. 

cfX1  = @(t) cfX_ChiSquare(t,5);
pdfX2 = @(x) normpdf(x);
a     = -8;
b     = 8;
tol   = 1e-4;
cf    = @(t) cf_Prod(t,cfX1,pdfX2,a,b,tol);
prob  = [0.9 0.95 0.99];
x = linspace(-30,30,201);

clear options;
options.N = 2^12;
figure;
result = cf2DistGP(cf,x,prob,options);

disp(result)

%% *RATIO OF TWO INDEPENDENT RVs*
% |cf_Ratio| 
%  Characteristic function of a ratio of two independent (continuous)
%  random variables, i.e. Y = X1 / X2. The distribution of the first random
%  variable X1 is given by its characteristic function (CF), say cf_X1(t),
%  and the distribution of the second random variable is given by its
%  probability density function (PDF), say pdf_X2(x). 
% 
%  Then, the characteristic function of the ratio of two independent
%  random variables, say cf_Y(t), is given by
%   cf_Y(t) = integral(@(x) cf_X1(t/x) .* pdf_X2(x),a,b),
%  where a and b are limits of integration (interval representing the
%  support of the distribution of X2). 
%
% SYNTAX:
%  cf = cf_Ratio(t,cfX1,pdfX2,a,b)


%% I.1.a EXAMPLE: CF of a ratio of two independent RVs  
%  % with central symmetric triangular distribution on [-1,1] and
%  % rectangular distributions on [-1,1]

t    = linspace(-20,20,2^10+1);
cfX1 = @(t) cfS_Triangular(t);
cf   = cf_Ratio(t,cfX1);
figure; plot(t,real(cf),t,imag(cf));grid on
title('CF of a ratio of triangular and rectangular RVs')

%% I.1.b EXAMPLE: CF of a ratio of two independent RVs 
%  % with chi-squared (df = 5) and symmetric rectangular distributions on
%  % [-1,1].

t = linspace(-5,5,2^10+1);
cfX1 = @(t) cfX_ChiSquare(t,5);
cf   = cf_Ratio(t,cfX1);
figure; plot(t,real(cf),t,imag(cf));grid on
title('CF of a ratio of chi-squared and rectangular RVs ')

%% I.1.c EXAMPLE: CF of a ratio of two independent RVs 
%  % with chi-squared (df = 5) distribution and the standard normal
%  % (Gaussian) distribution. 

t = linspace(-5,5,2^10+1);
cfX1  = @(t) cfX_ChiSquare(t,5);
pdfX2 = @(x) normpdf(x);
a     = -8;
b     =  8;
cf    = cf_Ratio(t,cfX1,pdfX2,a,b);
figure; plot(t,real(cf),t,imag(cf));grid on
title('CF of a ratio of chi-squared and normal RVs ')

%% I.1.d EXAMPLE: PDF/CDF of a ratio of two independent RVs 
%  % with central symmetric triangular distribution on [-1,1] and
%  % rectangular distributions on [-1,1]

cfX1  = @(t) cfS_Triangular(t);
a     = -1;
b     = 1;
pdfX2 = @(x) 1/(b-a);
tol   = 1e-4;
cf    = @(t) cf_Ratio(t,cfX1,pdfX2,a,b,tol);
x     = linspace(-30,30,1001);
prob  = [0.9 0.95 0.99];

clear options;
options.N = 2^10;
figure;
result = cf2DistGP(cf,x,prob,options);

disp(result)

%% I.1.e EXAMPLE: PDF/CDF of a ratio of two independent RVs 
%  % with chi-squared (df = 5) distribution and the standard normal
%  % (Gaussian) distribution. 

cfX1  = @(t) cfX_ChiSquare(t,5);
pdfX2 = @(x) normpdf(x);
a     = -8;
b     = 8;
tol   = 1e-4;
cf    = @(t) cf_Prod(t,cfX1,pdfX2,a,b,tol);
prob  = [0.9 0.95 0.99];
x = linspace(-30,30,201);

clear options;
options.N = 2^12;
result = cf2DistGP(cf,x,prob,options);

disp(result)

%% I.1.f EXAMPLE: PDF/CDF of a ratio of two independent RVs 
%  % with the standard normal (Gaussian) distribution and the chi-squared
%  (df = 5) distribution.

cfX1  = @(t) cfS_Gaussian(t);
pdfX2 = @(x) chi2pdf(x,5);
a     = 0;
b     = Inf;
tol   = 1e-4;
cf    = @(t) cf_Ratio(t,cfX1,pdfX2,a,b,tol);
x     = linspace(-3,3,501);
prob  = [0.9 0.95 0.99];

clear options;
options.N = 2^10;
figure;
result = cf2DistGP(cf,x,prob,options);

disp(result)