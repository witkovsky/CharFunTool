function pts = ChebPoints(N,interval)
% ChebPoints Evaluate N chebpoints on the given interval [a,b]
%  
% SYNTAX:
%   pts = ChebPoints(N,interval)
%
% EXAMPLE1 (Barycentric interpolant of the Sine function on (-pi,pi))
%   x    = ChebPoints(32,[-pi,pi])';
%   f    = sin(x);
%   xNew = linspace(-pi,pi,201)';
%   fNew = InterpBarycentric(x,f,xNew);
%   plot(x,f,xNew,fNew,'.')
%   disp([xNew fNew sin(xNew)])

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 24-Jul-2017 10:06:48

%% FUNCTION
%  pts = ChebPoints(N,interval)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 2);
if nargin < 2, interval = []; end

if isempty(interval)
    interval = [-1,1];
end

%% ALGORITHM

a   = interval(1);
b   = interval(2);
pts = (a+b)/2 - cos(pi*(0:N)/N)*(b-a)/2;
end