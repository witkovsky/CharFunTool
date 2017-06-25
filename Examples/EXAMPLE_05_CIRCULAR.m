%% CharFunTool / EXAMPLES
%  CHARACTERISTIC FUNCTIONS OF SELECTED CIRCULAR DISTRIBUTIONS
%
%  (c) 2017 Viktor Witkovsky, witkovsky@gmail.com
%  Version: 10-May-2017

clear
close all

%% *VON MISES DISTIBUTION*

% |cf_vonMises(t)| evaluates the characteristic function cf(t) of
% the von Mises distribution (circular normal distribution) with the
% parameters mu in (-pi,pi) and kappa > 0 (mu and 1/kappa are analogous to
% mu and sigma^2, the mean and variance in the normal distribution), on a
% circle e.g. the interval (-pi,pi), i.e.     
%   |cf_vonMises(t) = besseli(t,kappa)/besseli(0,kappa) .* exp(1i*t*mu)|.
% For more details see also WIKIPEDIA: 
% <https://en.wikipedia.org/wiki/Von_Mises_distribution>

%% I.a EXAMPLE: CF of the circular von Mises distribution on (-pi,pi)

 t = linspace(-10,10,501);
 cf = cf_vonMises(t);
 figure; plot(t,cf),grid
 title('CF of the uniform von Mises distribution on (-pi,pi)')

%% I.b EXAMPLE: CF of the mixture of the von Mises distribution on (-pi,pi)
 mu1 = 0;
 kappa1 = 5;
 mu2 = 1;
 kappa2 = 15;
 cf = @(t) 0.25*cf_vonMises(t,mu1,kappa1) + 0.75*cf_vonMises(t,mu2,kappa2);
 t = linspace(-10,10,501);
 figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
 title('CF of the von Mises distribution on (-pi,pi)')


%% I.c EXAMPLE: PDF/CDF of the von Mises distribution on (-pi,pi)
 mu = 0;
 kappa = 1;
 cf = @(t) cf_vonMises(t,mu,kappa);
 clear options
 options.xMin = -pi;
 options.xMax = pi;
 result = cf2DistGP(cf,[],[],options);
 disp(result);
 angle = result.x;
 radius = result.pdf;
 figure; polarplot(angle,radius);
 ax = gca; ax.ThetaAxisUnits = 'radians';

%% I.d EXAMPLE: PDF/CDF of the linear combinantion of 2 VM distributions

 mu1 = 0;
 kappa1 = 5;
 mu2 = 1;
 kappa2 = 15;
 c = [1 0.25];
 cf = @(t) cf_vonMises(c(1)*t,mu1,kappa1) .* cf_vonMises(c(2)*t,mu2,kappa2);
 clear options
 options.isCircular = true;
 result = cf2DistGP(cf,[],[],options);
 disp(result);
 angle = result.x;
 radius = result.pdf;
 figure; polarplot(angle,radius);
 ax = gca; ax.ThetaAxisUnits = 'radians';

%% I.e EXAMPLE: PDF/CDF of the mixture of 2 von Mises distribution

 mu1 = 0;
 kappa1 = 5;
 mu2 = 1;
 kappa2 = 15;
 mu3 = pi;
 kappa3 = 10;
 cf = @(t) 0.25*cf_vonMises(t,mu1,kappa1) + ...
      0.25*cf_vonMises(t,mu2,kappa2) + 0.5*cf_vonMises(t,mu3,kappa3);
 clear options
 options.isCircular = true;
 options.correctedCDF = true;
 x = linspace(0,2*pi);
 result = cf2DistGP(cf,x,[],options);
 disp(result);
 angle = result.x;
 radius = result.pdf;
 figure; polarplot(angle,radius);
 ax = gca; ax.ThetaAxisUnits = 'radians';
 
%% *COMBINATIONS OF CIRCULAR DISTIBUTION*
%% II.a EXAMPLE: PDF/CDF of the linear combinantion of wrapped distributions

 cf = @(t) cfS_Gaussian(0.25*t) .* cfS_Triangular(0.5*t) .* cfS_Arcsine(2*t);
 clear options
 options.isCircular = true;
 options.correctedCDF = true;
 result = cf2DistGP(cf,[],[],options);
 disp(result);
 angle = result.x;
 radius = result.pdf;
 figure; polarplot(angle,radius);
 ax = gca; ax.ThetaAxisUnits = 'radians';
 
 %% II.b EXAMPLE: PDF/CDF of the mixture of wrapped circular distributions
 df = 3;
 cf = @(t) 0.8*cfS_Student(t,df) + 0.2*cfS_Arcsine(pi*t/4);
 clear options
 options.isCircular = true;
 options.correctedCDF = true;
 options.N = 2^12;
  x = linspace(0,2*pi,501);
 result = cf2DistGP(cf,x,[],options);
 disp(result);
 angle = result.x;
 radius = result.pdf;
 figure; polarplot(angle,radius);
 ax = gca; ax.ThetaAxisUnits = 'radians';