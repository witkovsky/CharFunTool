function y = SqrtZX(z)
%SqrtZX 
%  Square root of a complex-valued vector z with using the unwrapped
%  arguments (angles), which is very important for correct evaluation of
%  functions (as e.g. the characteristic function) in the complex plane,
%  where the standard MATLAB SQRT function does not work as expected.
% 
%  In fact, for complex z the square root of z can be expressed as
%  y = sqrt(z) = sqrt(abs(z)).*exp(0.5i*angle(z)). Then, the corrected
%  verion is given by y = sqrt(abs(z)).*exp(0.5i*unwrap(angle(z))).
%
% SYNTAX:
%  y = SqrtZX(z)
% 
% INPUTS:
%  z     - vector or array of complex values, where the SQRT is evaluated.
%
% EXAMPLE 1
%  % CF of the asymptotic Anderson-Darling test statistic
%  t = linspace(-50,50,1000);
%  f = -SqrtZX((-2i*pi*t)./cos((pi/2).*SqrtZX(1+8i*t)));
%  figure;plot(t,real(f),t,imag(f));grid
%  title('CF of the asymptotic Anderson-Darling test statistic')
%
% EXAMPLE 2
%  % CF of the asymptotic Cramer-von Mises test statistic
%  t = linspace(-100,100,1000);
%  f = -SqrtZX(sqrt(2i*t)./sin(SqrtZX(2i*t)));
%  figure;plot(t,real(f),t,imag(f));grid
%  title('CF of the asymptotic Cramer-von Mises test statistic')

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 23-Oct-2018 20:52:43

%% ALGORITHM
%  cf = y = SqrtZX(z)

%%
sz = size(z);
z  = z(:);
y  = reshape(sqrt(abs(z)).*exp(0.5i*unwrap(angle(z))),sz);

end