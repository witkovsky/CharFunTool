%function cf = cfC_vonMises(t,mu,kappa)
%cfC_vonMises(t) evaluates the characteristic function cf(t) of
% the von Mises distribution (circular normal distribution) with the
% parameters mu in (-pi,pi) and kappa > 0 (mu and 1/kappa are analogous to
% mu and sigma^2, the mean and variance in the normal distribution), on a
% circle e.g. the interval (-pi,pi), i.e.     
%   cf(t) = cfC_vonMises(t,mu,kappa) 
%         = besseli(t,kappa)/besseli(0,kappa) .* exp(1i*t*mu).
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Von_Mises_distribution
%
% SYNTAX
%  cf = cfC_vonMises(t,mu,kappa) 
%
% EXAMPLE1 (CF of the uniform von Mises distribution on (-pi,pi))
%  t = linspace(-10,10,501);
%  cf = cfC_vonMises(t);
%  figure; plot(t,cf),grid
%  title('CF of the uniform von Mises distribution on (-pi,pi)')
%
% EXAMPLE2 (CF of the mixture of the von Mises distribution on (-pi,pi))
%  mu1 = 0;
%  kappa1 = 5;
%  mu2 = 1;
%  kappa2 = 15;
%  cf = @(t) 0.25*cfC_vonMises(t,mu1,kappa1) + 0.75*cfC_vonMises(t,mu2,kappa2);
%  t = linspace(-10,10,501);
%  figure; plot(t,real(cf(t)),t,imag(cf(t))),grid
%  title('CF of the von Mises distribution on (-pi,pi)')
%
% EXAMPLE3 (PDF/CDF of the von Mises distribution on (-pi,pi))
%  mu = 0;
%  kappa = 1;
%  cf = @(t) cfC_vonMises(t,mu,kappa);
%  clear options
%  options.xMin = -pi;
%  options.xMax = pi;
%  result = cf2DistGP(cf,[],[],options)
%  angle = result.x;
%  radius = result.pdf;
%  figure; polarplot(angle,radius);
%  ax = gca; ax.ThetaAxisUnits = 'radians';
%
% EXAMPLE4 (PDF/CDF of the linear combinantion of 2 von Mises distribution on (-pi,pi))
%  mu1 = 0;
%  kappa1 = 5;
%  mu2 = 1;
%  kappa2 = 15;
%  c = [1 0.25];
%  cf = @(t) cfC_vonMises(c(1)*t,mu1,kappa1) .* cfC_vonMises(c(2)*t,mu2,kappa2);
%  clear options
%  options.isCircular = true;
%  result = cf2DistGP(cf,[],[],options)
%  angle = result.x;
%  radius = result.pdf;
%  figure; polarplot(angle,radius);
%  ax = gca; ax.ThetaAxisUnits = 'radians';
%
% EXAMPLE5 (PDF/CDF of the mixture of 2 von Mises distribution on (0,2*pi))
%  mu1 = 0;
%  kappa1 = 5;
%  mu2 = 1;
%  kappa2 = 15;
%  mu3 = pi;
%  kappa3 = 10;
%  cf = @(t) 0.25*cfC_vonMises(t,mu1,kappa1) + ...
%       0.25*cfC_vonMises(t,mu2,kappa2) + 0.5*cfC_vonMises(t,mu3,kappa3);
%  clear options
%  options.isCircular = true;
%  options.correctedCDF = true;
%  x = linspace(0,2*pi);
%  result = cf2DistGP(cf,x,[],options)
%  angle = result.x;
%  radius = result.pdf;
%  figure; polarplot(angle,radius);
%  ax = gca; ax.ThetaAxisUnits = 'radians';
%
% REFERENCES:


% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfC_vonMises(t,mu,kappa)
