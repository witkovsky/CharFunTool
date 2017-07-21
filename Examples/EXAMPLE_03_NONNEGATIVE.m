%% CharFunTool / EXAMPLES
%  CHARACTERISTIC FUNCTIONS OF SELECTED NONNEGATIVE CONTINUOUS DISTRIBUTIONS
%
%  (c) 2017 Viktor Witkovsky, witkovsky@gmail.com
%  Version: 10-May-2017

clear
close all

%% *EXPONENTIAL DISTIBUTION*
% |cfX_Exponential(t,lambda)| evaluates the characteristic function cf(t)
% of the Exponential distribution with the parameter lambda (rate, lambda >
% 0), i.e.   
%
% |cf(t) = cfX_Exponential(t,lambda) = lambda / (lambda - 1i*t)|
%
% For more details see WIKIPEDIA: 
% <https://en.wikipedia.org/wiki/Exponential_distribution>
%
% SYNTAX:
%
% * |cf = cfX_Exponential(t,lambda)|

clear 
close all

%% III.1.a EXAMPLE: CF of the Exponential distribution, lambda = 5
%
 lambda = 5;  
 t = linspace(-50,50,501);
 cf = cfX_Exponential(t,lambda);
 figure; plot(t,real(cf),t,imag(cf)),grid
 title('CF of the Exponential distribution with lambda = 5')

%% III.1.b EXAMPLE: PDF/CDF of the Exponential distribution, lambda = 5
%
 lambda = 5;  
 cf = @(t) cfX_Exponential(t,lambda);
 x  = linspace(0,1.5,101);
 clear options
 options.xMin = 0;
 options.SixSigmaRule = 8;
 result = cf2DistGP(cf,x,[],options);
 disp(result)
 
%% III.1.c EXAMPLE: PDF/CDF of the compound Binomial-Exponential distribution
%
 n = 25;  
 p = 0.3;
 lambda = 5;
 cfX  = @(t) cfX_Exponential(t,lambda);
 cf   = @(t) cfN_Binomial(t,n,p,cfX);
 x    = linspace(0,5,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)

%% *FISHER-SNEDECOR F DISTIBUTION*
% |cfX_FisherSnedecor(t,df1,df2,tol)| evaluates the characteristic
% function cf(t) of the Fisher-Snedecor F-distribution with the parameters
% df1 and df2 (degrees of freedom, df1>=1, df2>=1) computed for real vector
% argument t, i.e. 
%
% |cf(t) = cfX_FisherSnedecor(t,df1,df2)  
%        = (gamma(df1/2+df2/2)/gamma(df2/2)) * U(df1/2, 1-df2/2, -1i*(df2/df1)*t)|,
%
%  where U(a,b,z) denotes the confluent hypergeometric function of the
%  second kind.  For more details see WIKIPEDIA: 
%  <https://en.wikipedia.org/wiki/F-distribution> 
%
% SYNTAX:
%
% * |cf = cfX_FisherSnedecor(df1,df2,t,tol)|

%% III.2.a EXAMPLE: CF of the F-distribution with df1 = 3, df2 = 5
%
 df1 = 3;
 df2 = 5;
 t   = linspace(-30,30,2^10+1)';
 cf  = cfX_FisherSnedecor(t,df1,df2);
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the F-distribution')

%% III.2.b EXAMPLE: PDF/CDF of the F-distribution with df1 = 3, df2 = 5
%
 df1 = 3;
 df2 = 5;
 x = linspace(0,25,101);
 prob = [0.9 0.95 0.99];
 cf = @(t) cfX_FisherSnedecor(t,df1,df2);
 clear options
 options.xMin = 0;
 options.xMax = 500;
 options.N  = 2^15;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)

%% III.2.c EXAMPLE: PDF/CDF of the compound Binomial-Fisher-Snedecor distribution
%
 n = 25;  
 p = 0.3;
 df1 = 3;
 df2 = 5;
 cfX = @(t) cfX_FisherSnedecor(t,df1,df2);
 cf = @(t) cfN_Binomial(t,n,p,cfX);
 x = linspace(0,80,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% *GAMMA DISTRIBUTION*
% |cfX_Gamma(t,alpha,beta)| evaluates the characteristic function cf(t) of 
% the GAMMA distribution with the parameters alpha (shape, alpha > 0) and
% beta (rate, beta > 0), i.e.
%
% |cf(t) = cfX_Gamma(t,alpha,beta) = (1 - 1i*t / beta)^(-alpha)|;
%
% For more details see WIKIPEDIA: 
% <https://en.wikipedia.org/wiki/Gamma_distribution>
%
% SYNTAX:
%
% * |cf = cfX_Gamma(t,alpha,beta)|

%% III.3.a EXAMPLE: CF of the Gamma distribution, alpha = 2, beta = 2
%
 alpha = 2;
 beta = 2;
 t = linspace(-20,20,501);
 cf = cfX_Gamma(t,alpha,beta);
 figure; plot(t,real(cf),t,imag(cf)),grid
 title('CF of the Gamma distribution with alpha = 2, beta = 2')

%% III.3.b EXAMPLE: PDF/CDF of the Gamma distribution, alpha = 2, beta = 2
%
 alpha = 2;
 beta = 2;
 x = linspace(0,5,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.xMin = 0;
 options.N = 2^14;
 cf = @(t) cfX_Gamma(t,alpha,beta);
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% III.3.c EXAMPLE: PDF/CDF of the compound Binomial-Gamma distribution
%
 n = 25;  
 p = 0.3;
 alpha = 2;
 beta = 2;
 cfX = @(t)  cfX_Gamma(t,alpha,beta);
 cf = @(t) cfN_Binomial(t,n,p,cfX);
 x = linspace(0,25,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% *GENERALIZED PARETO DISTRIBUTION*
% |cfX_GeneralizedPareto(t,xi,sigma,theta)| Computes the characteristic
% function cf(t) of the Generalized Pareto distribution with parameters xi
% (shape, here xi >= 0), sigma (scale, sigma > 0), and theta (threshold,
% theta >= 0), for real (vector) argument t, i.e.  
%
% |cf(t) = cfX_GeneralizedPareto(t,xi,sigma,theta)|;
%
% The closed-form analytic expression of the characteristic function of the
% Generalized Pareto distribution is unknown. Thus, the characteristic
% function is numerically evalueted from its definition. 
% For more details see also WIKIPEDIA:    
% <https://en.wikipedia.org/wiki/Generalized_Pareto_distribution>
%  
% SYNTAX:
%
% * |cf = cfX_GeneralizedPareto(t,xi,sigma,theta)|
% * |cf = cfX_GeneralizedPareto(t,xi,sigma,theta,tol)|

%% III.4.a EXAMPLE: CF of the GeneralizedPareto distribution
%
 xi = 1;
 sigma = 1;
 theta = 1;
 t = linspace(-20,20,2^10+1)';
 cf = cfX_GeneralizedPareto(t,xi,sigma,theta);
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the GeneralizedPareto distribution')

%% III.4.b EXAMPLE: CDF/PDF of the  GeneralizedPareto distribution
%
 xi = 1;
 sigma = 1;
 theta = 0;
 x = linspace(0,300,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.N = 2^15;
 options.SixSigmaRule = 15;
 cf = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
 result = cf2DistGP(cf,x,prob,options);
 disp(result)

%% III.4.c EXAMPLE: PDF/CDF of the compound Poisson-GeneralizedPareto distribution
%
 xi = 1;
 sigma = 1;
 theta = 0;
 lambda = 10;
 cfX = @(t) cfX_GeneralizedPareto(t,xi,sigma,theta);
 cf = @(t) cfN_Poisson(t,lambda,cfX);
 x = linspace(0,5000,101);
 prob = [0.9 0.95 0.99];
 clear options
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
  
%% *CHISQUARED DISTRIBUTION*
% |cfX_ChiSquare(t,df)| evaluates the characteristic function cf(t) of 
% the CHI-SUQARED distribution with the parameter df (degrees of freedom,
% df >= 1), i.e. 
%
% |cf(t) = cfX_ChiSquare(t,df) = (1 - 2i*t)^(-df/2) =
%         = cfX_Gamma(t,alpha,beta) = cfX_Gamma(t,df/2,1/2)|;
% For more details see WIKIPEDIA: 
% <https://en.wikipedia.org/wiki/Chi-squared_distribution>
%
% SYNTAX:
%
% * |cf = cfX_ChiSquare(t,df)|

%% III.5.a EXAMPLE: CF of the ChiSquared distribution, df = 1
%
 df = 1;
 t = linspace(-50,50,501);
 cf = cfX_ChiSquare(t,df);
 figure; plot(t,real(cf),t,imag(cf)),grid
 title('CF of the ChiSquared distribution with df = 1')

%% III.5.b EXAMPLE: PDF/CDF of the ChiSquared distribution, df = 3
%
 df = 3;
 x = linspace(0,15,101);
 prob = [0.9 0.95 0.99];
 cf = @(t) cfX_ChiSquare(t,df);
 clear options
 options.xMin = 0;
 options.xMax = 22;
 options.N = 2^14;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)

%% III.5.c EXAMPLE: PDF/CDF of the compound Binomial-ChiSquared distribution
%
 n = 25;  
 p = 0.3;
 df = 3;
 cfX = @(t) cfX_ChiSquare(t,df);
 cf = @(t) cfN_Binomial(t,n,p,cfX);
 x = linspace(0,80,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% *LOG-LOGISTIC DISTRIBUTION*
% |cfX_LogLogistic(t,alpha,beta)| computes the characteristic function
% cf(t) of the LogLogistic distribution with parameters alpha (scale,
% alpha > 0) and beta (shape, beta > 0), for real (vector) argument t, i.e.
% 
% |cf(t) = cfX_LogLogistic(t,alpha,beta)|;
%
% The closed-form analytic expression of the characteristic function of the
% LogLogistic distribution is unknown. Thus, the characteristic function is
% numerically evalueted from its definition. 
% For more details see WIKIPEDIA:     
% <https://en.wikipedia.org/wiki/Log-logistic_distribution>
%  
% SYNTAX:
%
% * |cf = cfX_LogLogistic(t,alpha,beta)|
% * |cf = cfX_LogLogistic(t,alpha,beta,tol)|

%% III.6.a EXAMPLE: CF of the Loglogistic distribution, alpha=1, beta=1
%
 alpha = 1;
 beta = 1;
 t = linspace(-20,20,2^10+1)';
 cf = cfX_LogLogistic(t,alpha,beta);
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the LogLogistic distribution')

%% III.6.b EXAMPLE: CDF/PDF of the  Loglogistic distribution, alpha=1, beta=1
%
 alpha = 1;
 beta = 1;
 x = linspace(0,700,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.xMin = 0;
 options.xMax = 1000;
 options.N = 2^14;
 cf = @(t) cfX_LogLogistic(t,alpha,beta);
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% III.6.c EXAMPLE: PDF/CDF of the compound Poisson-Loglogistic distribution
%
 alpha = 1;
 beta = 1;
 lambda = 10;
 cfX = @(t)cfX_LogLogistic(t,alpha,beta);
 cf = @(t) cfN_Poisson(t,lambda,cfX);
 x = linspace(0,3500,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
  disp(result)
  
%% *LOG-NORMAL DISTRIBUTION*
% |cfX_LogNormal(t,mu,sigma)| computes the characteristic function cf(t) of
% the LogNormal distribution with parameters mu (real) and sigma > 0,
% computed for real (vector) argument t, i.e.
%
% |cf(t) = cfX_LogNormal(t,mu,sigma)|;
% For more details see WIKIPEDIA:
% <https://en.wikipedia.org/wiki/Log-normal_distribution>
%  
% SYNTAX:
%
% * |cf = cfX_LogNormal(t,mu,sigma)|
% * |cf = cfX_LogNormal(t,mu,sigma,tol)|

%% III.7.a EXAMPLE: CF of the LogNormal distribution, mu = 0, sigma = 1
%
 mu = 0;
 sigma = 1;
 t = linspace(-20,20,2^10+1)';
 cf = cfX_LogNormal(t,mu,sigma);
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the LogNormal distribution')

%% III.7.b EXAMPLE: CDF/PDF of the LogNormal distribution, mu = 0,sigma = 1
%
 mu = 0;
 sigma = 1;
 x = linspace(0,15,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.xMin = 0;
 options.N = 2^13;
 options.SixSigmaRule = 8;
 cf = @(t) cfX_LogNormal(t,mu,sigma);
 result = cf2DistGP(cf,x,prob,options);
 disp(result)

%% III.7.c EXAMPLE: PDF/CDF of the compound Poisson-LogNormal distribution
%
 mu = 0;
 sigma = 1;
 lambda = 10;
 cfX = @(t)cfX_LogNormal(t,mu,sigma);
 cf = @(t) cfN_Poisson(t,lambda,cfX);
 x = linspace(0,70,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)

%% *PARETO DISTRIBUTION*
% |cfX_Pareto(t,alpha,sigma,type,tol)| computes the characteristic function
% of the Pareto distribution with parameters alpha > 0 (shape), sigma > 0
% (scale, sometimes denoted by x_m - which is the minimum value of the
% distribution support, i.e. x >= x_m), and the parameter type ('typeI' or
% 'typeII'), computed for real vector argument t, i.e.
%
% |cf(t) = cfX_Pareto(t,alpha,sigma,type)|;
%
% For more details see WIKIPEDIA:
% <https://en.wikipedia.org/wiki/Pareto_distribution>
%
% SYNTAX:
%
% * |cf = cfX_Pareto(t,alpha)|
% * |cf = cfX_Pareto(t,alpha,sigma,type)|
%
% FURTHER DETAILS: 
%
% The Pareto Type I distribution (also known as European Pareto) is
% characterized by a scale parameter sigma and a shape parameter alpha,
% which is known as the tail index. The pdf of X is defined by pdf_X(x) =
% alpha * sigma^alpha / x^(alpha+1), for x >= sigma, otherwise it is zero.
% When this distribution is used to model the distribution of wealth, then
% the parameter alpha is called the Pareto index.
%
% The Type II Pareto distribution (also known as Lomax or American Pareto
% distribution, type = 'typeII' or type = 'Amer') is shifted to zero.
% The default type is 'TypeI' (type = 'typeI' or type = 'Eur').
% The CFs of the Pareto type distributions are defined by
% cf_{typeI}(t)  = alpha * exp(1i*sigma*t) * U(1,1-alpha,-1i*sigma*t)
% cf_{typeII}(t) = alpha * U(1,1-alpha,-1i*sigma*t)
%
% A special case of Pareto TypeI distribution is the Pareto 1 distribution
% with sigma = 1, cf_{1}(t) = alpha * exp(1i*t) * U(1,1-alpha,-1i*t), where
% U(a,b,z) denotes the confluent hypergeometric function of the second
% kind. 

%% III.8.a EXAMPLE: CF of the Pareto 1 distribution, alpha = 3/2
%
 alpha = 3/2;  
 t = linspace(-10,10,1001);
 cf = cfX_Pareto(t,alpha);
 figure; plot(t,real(cf),t,imag(cf)),grid
 title('CF of the Pareto 1 distribution with alpha = 3/2')

%% III.8.b EXAMPLE: CF of the Pareto Type I/EUR distribution
%
 alpha = 3/2;  
 sigma  = 2;
 type = 'eur';
 t = linspace(-10,10,1001);
 cf = cfX_Pareto(t,alpha,sigma,type);
 figure; plot(t,real(cf),t,imag(cf)),grid
 title('CF of the Pareto EUR distribution with alpha = 3/2, sigma = 2')

%% III.8.c EXAMPLE: CF of the Pareto Type II/AMER distribution
%
 alpha = 3/2;  
 sigma  = 2;
 type = 'amer';
 t = linspace(-10,10,1001);
 cf = cfX_Pareto(t,alpha,sigma,type);
 figure; plot(t,real(cf),t,imag(cf)),grid
 title('CF of the Pareto AMER distribution with alpha = 3/2, sigma = 2')

%% III.8.d EXAMPLE: PDF/CDF of the compound Poisson-Pareto distribution
%
 alpha = 3/2;  
 sigma  = 2;
 type = 'eur';
 lambda = 10;
 cfX = @(t) cfX_Pareto(t,alpha,sigma,type);
 cf = @(t) cfN_Poisson(t,lambda,cfX);
 x = linspace(0,500,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)

%% *PEARSON TYPE VI DISTRIBUTION*
% |cfX_PearsonVI(t,alpha,beta)| omputes the characteristic function of the
% Pearson type VI distribution with parameters alpha > 0 (shape) and beta >
% 0 (scale), computed for real vector argument t, i.e. 
%
% |cf(t) = cfX_PearsonVI(t,alpha,beta)
%         = (gamma(alpha+beta)./gamma(beta)) .* U(alpha,1-beta,-1i*t)|
%
% where U(a,b,z) denotes the confluent hypergeometric function of the
% second kind. For more detail see WIKIPEDIA:
% <https://en.wikipedia.org/wiki/Pearson_distribution>
%
% SYNTAX:
%
% * |cf = cfX_PearsonVI(t,alpha,beta)|

%% III.9.a EXAMPLE: CF of the Pearson Type VI distribution
%
 alpha = 3/2;
 beta = 2/3;
 t = linspace(-10,10,1001);
 cf = cfX_PearsonVI(t,alpha,beta);
 figure; plot(t,real(cf),t,imag(cf)),grid
 title('CF of the Pearson Type VI distribution with alpha = 3/2, beta = 2/3')

%% III.9.b EXAMPLE: PDF/CDF of the Pearson Type VI distribution
%
 alpha = 3/2;
 beta = 2/3;
 x = linspace(0,200,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.xMin = 0;
 options.N = 2^14;
 cf = @(t) cfX_PearsonVI(t,alpha,beta);
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% III.9.c EXAMPLE: PDF/CDF of the compound Binomial-Pearson Type VI distribution
%
 n = 25;  
 p = 0.3;
 alpha = 3/2;
 beta = 2/3;
 cfX = @(t) cfX_PearsonVI(t,alpha,beta);
 cf = @(t) cfN_Binomial(t,n,p,cfX);
 x = linspace(0,10000,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% *WEIBULL DISTRIBUTION*
% |cfX_Weibull(t,alpha,beta)| computes the characteristic function cf(t)
% of the Weibull distribution with the parameters alpha(scale parameter,
% elsewhere denoted also as lambda, alpha > 0) and beta (shape parameter,
% elsewhere denoted as k, beta > 0), for real (vector) argument t, i.e.
%
% |cf(t) = cfX_Weibull(t,alpha,beta);|
%
% The closed-form analytic expression of the characteristic function of the
% Weibull distribution is unknown. The series expansion of the CF, valid
% for abs(t)<1, is known and given, e.g. in WIKIPEDIA. 
% 
% Thus, for abs(t)<1 the CF is evalueted from its series expansion, and for
% abs(t)>=1 from its definition. For more details see WIKIPEDIA:  
% <https://en.wikipedia.org/wiki/Weibull_distribution>
%  
% SYNTAX:
%
% * |cf = cfX_Weibull(t,alpha,beta)|
% * |cf = cfX_Weibull(t,alpha,beta,tol)|
%
% REMARK (Caution Notice)
%
% This implementation of the algorithm is experimental. It DOES NOT WORK
% PROPERLY for all parameters!!!!   

%% III.10.a EXAMPLE CF of the Weibull distribution, alpha=1, beta=1
%
 alpha = 1;
 beta = 1;
 t = linspace(-20,20,2^10+1)';
 cf = cfX_Weibull(t,alpha,beta);
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the Weibull distribution')

%% III.10.b EXAMPLE: CDF/PDF of the  Weibull distribution, alpha=1, beta=1
%
 alpha = 1;
 beta = 1;
 x = linspace(0,7,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.xMin = 0;
 options.N = 2^10;
 cf = @(t) cfX_Weibull(t,alpha,beta);
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% III.10.c EXAMPLE: PDF/CDF of the compound Poisson-Weibull distribution
%
 alpha = 1;
 beta = 1;
 lambda = 10;
 cfX = @(t)cfX_Weibull(t,alpha,beta);
 cf = @(t) cfN_Poisson(t,lambda,cfX);
 x = linspace(0,35,101);
 prob = [0.9 0.95 0.99];
 clear options
 options.isCompound = true;
 result = cf2DistGP(cf,x,prob,options);
 disp(result)
 
%% *CHARACTERISTIC FUNCTION FROM PDF*
% |cfX_PDF(t,pdfFun,method,tol)| computes the characteristic function of
% the continuos distribution defined by its PDF function, pdfFun = @(x)
% pdf(x), here we assume x>=0, computed for real vector argument t.
%
% DEFINITION:
%  cfX_PDF is based on the standard integral representation of the
%  characteristic function of the continuous distribution defined by its 
%  PDF (here PDF is represented by the function handle pdfFun(x) defined 
%  for x >= 0). Then, 
%    CF(t) = Integral_0^inf exp(i*t*x) * pdfFun(x) dx.
%  Alternatively, by using the half-space Fourier Integral Transformation
%  (FIT), which could improve the highly oscillatory behaviour of the
%  integrand function, we get    
%    CF(t) = Integral_0^inf (i/t) * exp(-x) * pdfFun(i*x/t) dx.
%  If we define the integrand as 
%    funCF(t,x) = (i/t) * exp(-x) * pdfFun(i*x/t),
%  then by using a stabilizing transformation from [0,inf] to [0,1], we can
%  evaluate the CF by the following (well behaved) integral:
%    CF(t) = Integral_0^1 2x/(1-x)^3 * funCF(t,(x/(1-x))^2) dx.
%
%  Selection of the proper method (standard definition of FIT) depends on
%  the distribution and particular form of its PDF.  
%
%  cfX_PDF evaluates this integral by using the MATLAB built in function
%  'integral', with precission specified by tolerance tol (default value is 
%  tol = 1e-6).

%% III.11.a EXAMPLE CF from PDF of the Exponential distribution, mu = 1
%
 pdfFun = @(x) exp(-x);
 t = linspace(-20,20,2^10+1)';
 cf = cfX_PDF(t,pdfFun);
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the Exponential distribution')

%% III.11.b EXAMPLE: CF from PDF of the LogNormal distribution, mu = 0, sigma = 1
%
 mu = 0;
 sigma = 1;
 pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
 t = linspace(-20,20,2^10+1)';
 cf = cfX_PDF(t,pdfFun,'fit');
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the LogNormal distribution')

%% III.11.c EXAMPLE: CDF/PDF of the LogNormal distribution, mu = 0, sigma = 1
%
 mu = 0;
 sigma = 1; 
 pdfFun = @(x) exp(-0.5*((log(x)-mu)./sigma).^2)./(x.*sqrt(2*pi).*sigma);
 cf = @(t) cfX_PDF(t,pdfFun,'fit');
 clear options
 options.xMin = 0;
 result = cf2DistGP(cf,[],[],options);
 disp(result)
 
%% III.11.d EXAMPLE: CF from PDF of the Weibull distribution, a = 1.5, b = 1
%
 a = 1.5;
 b = 1;
 pdfFun = @(x) (x./a).^(b-1) .* exp(-((x./a).^b)) .* b ./ a;
 t = linspace(-20,20,2^10+1)';
 tol = 1e-10;
 cf = cfX_PDF(t,pdfFun,'fit',tol);
 plot(t,real(cf),t,imag(cf));grid
 title('Characteristic function of the Weibull distribution')