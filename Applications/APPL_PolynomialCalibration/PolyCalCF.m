function cf = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,shift)
% PolyCalCF is an auxiliary function for the algorithm PolyCal which
% creates the function handle of the combined 'calibration' characteristic
% function of the random variable W, which is a linear function of the
% input varibles defined by W = wX'*XA + wX'*XB + sum(wX)*XB0 + wY'*YA +
% wY'*YB + sum(wY)*YB0 + shift. 
% That is, given the characteristic functions of the input variables XA,
% XB, XB0, YA, YB, YB0 (assumed that all their components are mulually
% independent random variables), and the vectors of coefficients wX  
% and wY, the characteristic function of W is defined as     
% cfW = @(t)(prod_{i=1}^m cfXA{i}(wX(i)*t)) .* ...
%           cfXB{i}(wX(i)*t)) .* cfXB0((sum_{i=1}^m wX(i))*t)) .* ... 
%           (prod_{i=1}^m cfYA{i}(wY(i)*t)) .* ...
%           cfYB{i}(wY(i)*t)) .* cfYB0((sum_{i=1}^m wY(i))*t)) .* ...
%           exp(1i*t*shift),  
% where wX and wY are m-dimensional vectors of coefficients, and cfXA =
% {@(t)cfXA{1}(t), ..., @(t)cfXA{m}(t)}, cfXB = {@(t)cfXB{1}(t), ...,
% @(t)cfXB{m}(t)}, cfYA = {@(t)cfYA{1}(t), ..., @(t)cfYA{m}(t)}, cfYB =
% {@(t)cfYB{1}(t), ..., @(t)cfYB{m}(t)}, are m-dimensional cell vectors of
% the characteristic functions of XA = (XA_1,...,XA_m), XB =
% (XB_1,...,XB_m), YA = (YA_1,...,YA_m), and YB = (YB_1,...,YB_m),
% respectively, and cfXB0 and cfYB0 are the characteristic functions of the
% 'common' random variables XB0 and YB0, respectively.
%
% SYNTAX
% cf = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,shift)
%
% INPUTS:
%  wX    - m-dimensional vector of coefficients which define the linear
%          combinations of the random vectors XA, XB, and with the random
%          variable XB0.  
%  wY    - m-dimensional vector of coefficients which define the linear
%          combinations of the random vectors YA, YB, and with the random
%          variable YB0.  
%  wX    - m-dimensional vector of coefficients which define the linear
%          combinations of the random vectors XA, XB, and with the random
%          variable XB0. 
% cfXA   - m-dimensional cell with function handles of the characteristic
%          functions of independent components of XA = (XA_1,...,XA_m).
% cfXB   - m-dimensional cell with function handles of the characteristic
%          functions of independent components of XB = (XB_1,...,XB_m).
% cfXB0  - function handle of the characteristic function of the random
%          variable XB0. 
% cfYA   - m-dimensional cell with function handles of the characteristic
%          functions of independent components of YA = (YA_1,...,YA_m).
% cfYB   - m-dimensional cell with function handles of the characteristic
%          functions of independent components of YB = (YB_1,...,YB_m).
% cfYB0  - function handle of the characteristic function of the random
%          variable YB0.
%
% OUTPUT:
% cf     - function handle of the combined characteristic function.
%
% EXAMPLE (CF of W=wX'*(XA+XB)+sum(wX)*XB0 + wY'*(YA+YB)+sum(wY)*YB0+shift)
%  wX = [1 2 3 4 5]/15;
%  wY = [1 1 1 1 1]/5;
%  shift = 1;
%  cfXA  = {@(t)cf_Normal(t), ...
%           @(t)cf_Normal(1.1*t), ...
%           @(t)cf_RectangularSymmetric(1.2*t), ...
%           @(t)cf_RectangularSymmetric(2*t), ...
%           @(t)cf_RectangularSymmetric(8*t)}';
%  cfXB  = {@(t)cf_RectangularSymmetric(t), ...
%           @(t)cf_RectangularSymmetric(1.1*t), ...
%           @(t)cf_RectangularSymmetric(1.2*t), ...
%           @(t)cf_RectangularSymmetric(2*t), ...
%           @(t)cf_RectangularSymmetric(8*t)}';
%  cfXB0 =  @(t)cf_ArcsineSymmetric(t);
%  cfYA  = {@(t)cf_Student(t,5), ...
%           @(t)cf_Normal(1.1*t), ...
%           @(t)cf_TriangularSymmetric(1.2*t), ...
%           @(t)cf_TriangularSymmetric(2*t), ...
%           @(t)cf_TriangularSymmetric(8*t)}';
%  cfYB  = {@(t)cf_RectangularSymmetric(t), ...
%           @(t)cf_RectangularSymmetric(1.1*t), ...
%           @(t)cf_RectangularSymmetric(1.2*t), ...
%           @(t)cf_RectangularSymmetric(2*t), ...
%           @(t)cf_RectangularSymmetric(8*t)}';
%  cfYB0 =  @(t)cf_ArcsineSymmetric(t);
%  cf    = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,shift);
%  figure
%  t = linspace(-5,5,201);
%  plot(t,real(cf(t)))
%  figure 
%  x      = linspace(-10,10,201);
%  prob   = [0.9 0.95 0.99];
%  result = cf2DistGP(cf,x,prob);
%  disp(result)
%
% REFERENCES:
% [1] WITKOVSKY V. and WIMMER G.: PolyCal - MATLAB algorithm for
% comparative polynomial calibration and its applications. AMCTM 2020. 

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 28-Sep-2020 16:46:23

%% Algorithm
narginchk(2, 9);
if nargin < 9, shift = []; end
if nargin < 8, cfYB0 = []; end
if nargin < 7, cfYB = []; end
if nargin < 6, cfYA = []; end
if nargin < 5, cfXB0 = []; end
if nargin < 4, cfXB = []; end
if nargin < 3, cfXA = []; end
if isempty(shift), shift = 0; end

m  = length(wX);
cf = @(t) cfXB0(sum(wX)*t) .* cfYB0(sum(wY)*t);
for i = 1:m
    cf = @(t) cf(t) .* cfXA{i}(wX(i)*t)  .* cfXB{i}(wX(i)*t) .* ...
        cfYA{i}(wY(i)*t) .* cfYB{i}(wY(i)*t);    
end

if ~isempty(shift)
    cf = @(t) cf(t) .* exp(1i*t*shift);
end
end