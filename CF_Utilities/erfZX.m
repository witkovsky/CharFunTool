function f = erfZX(z)
% erfZX
%  Computes the error function in complex argument z by using the
%  Faddeeva function.
%
%  The erfZX function is defined as
%   erfZX(z) = 1 - exp(-z^2) * w(1i*z)
%  where w(z) is the Faddeeva function or the Kramp function, which is a
%  scaled complex complementary error function, in complex argument z.
%
%  SYNTAX:
%   f = erfZX(z)
%
%  INPUTS:
%   z      - possibly complex input argument (vector or array),
%
%  OUTPUTS:
%   f      - values of the erfZX function
%
%  EXAMPLE 1:
%  % erfZX the error function in complex argument z
%  z = 1i + linspace(-5,5);
%  f = erfZX(z);
%  plot(f)

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 20-Sep-2019 20:26:55
%% FUNCTION CALL
% f = erfZX(z)

%% ALGORITHM
if isreal(z)
    f = erf(z);
else
    f = 1 - erfcZX(z);
end
end