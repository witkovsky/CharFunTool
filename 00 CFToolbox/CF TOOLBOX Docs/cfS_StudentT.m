%function cf = cfS_StudentT(t,df)
%cfS_StudentT(t,df) evaluates the characteristic function cf(t) of
% the Student's t-distribution with the parameter df (degrees of freedom,
% df>=1) computed for real vector argument t, i.e.   
%   cf(t) = cfS_StudentT(t,df) 
%         = besselk(df/2,abs(t)*sqrt(df),1) * exp(-abs(t)*sqrt(df)) * ...
%          (sqrt(df)*abs(t))^(df/2) / 2^(df/2-1)/gamma(df/2);
% For more details see also WIKIPEDIA: 
% https://en.wikipedia.org/wiki/Student%27s_t-distribution
%
% SYNTAX
%  cf = cfS_StudentT(t,df)
%
% EXAMPLE1 (CF of the Student t-distribution with df = 3)
%  df = 3;
%  t = linspace(-5,5,501);
%  cf = cfS_StudentT(t,df);
%  figure; plot(t,cf),grid
%  title('CF of the Student t-distribution with df = 3)')
%
% EXAMPLE2 (PDF/CDF of the Student t-distribution with df = 3)
%  df = 3;
%  cf = @(t) cfS_StudentT(t,df);
%  x = linspace(-8,8,101);
%  clear options
%  options.SixSigmaRule = 8;
%  result = cf2DistGP(cf,x,[],options)
%
% REFERENCES:
% [1] WITKOVSKY V. (2016). Numerical inversion of a characteristic
%     function: An alternative tool to form the probability distribution of
%     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.  

% (c) 2016 Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 15-Nov-2016 13:36:26

%% ALGORITHM
%cf = cfS_StudentT(t,df);