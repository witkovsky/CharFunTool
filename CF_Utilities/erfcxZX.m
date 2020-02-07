function f = erfcxZX(z)
% erfcxZX
%  Computes the scaled complementary error function in complex argument z
%  by using the Faddeeva function.
%
%  The erfcxZX function is defined as
%   erfcxZX(z) = w(1i*z)
%  where w(z) is the Faddeeva function or the Kramp function, which is a
%  scaled complex complementary error function, in complex argument z.
%
%  SYNTAX:
%   f = erfcxZX(z)
%
%  INPUTS:
%   z      - possibly complex input argument (vector or array),
%
%  OUTPUTS:
%   f      - values of the erfcxZX function
%
%  EXAMPLE 1:
%  % erfcxZX the scaled complementary error function in complex argument z
%  z = 1i + linspace(-5,5);
%  f = erfcxZX(z);
%  plot(f)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Sep-2019 20:26:55
%% FUNCTION CALL
% f = erfcxZX(z)

%% ALGORITHM
if isreal(z)
    f = erfcx(z);
else
    f = Faddeeva(1i*z);
end
end