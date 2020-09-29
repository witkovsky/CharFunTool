function result = PolyCalFit(meanXNew,cfXNew,PolyCalResult,options)
% PolyCalFit computes the best estimate and the associated measurement
%  uncertainty together with the state-of-knowledge distribution of a
%  response variable YNew corresponding to a new stimulus variable XNew
%  based on the fitted polynomial calibration function resulted from the
%  comparative polynomial calibration experiment.
%
%  Here, we assume that the new stimulus XNew is a random variable fully
%  specified by its state-of-knowledge distribution which is represented by
%  its mean value meanXNew and the zero-mean probability distribution,
%  given by its characteristic function cfXNew.
%
%  Based on the state-of-knowledge distribution of XNew and the computed
%  results from the comparative polynomial calibration (PolyCalResult), the
%  algorithm computes the best estimate of the response variable YNew
%  (meanYNew) with the associated uncertainty (stdYNew) and the
%  state-of-knowledge distribution, specified by its characteristic
%  function (cfYNew).
%
% SYNTAX:
% result = PolyCalFit(meanXNew,cfXNew,PolyCalResult,options)
%
% INPUTS:
%  meanXNew      - n-dimensional vector of mean values of the new stimulus
%                  variables XNew. 
%  cfXNew        - n-dimensional cell with function handles of
%                  the characteristic functions of the zero-mean
%                  state-of-knowledge distributions of the components of
%                  the zero-mean shifted stimulus variable XNew-meanXNew.
%  PolyCalResult - is a structure with the calibration results generated by
%                  the algorithm PolyCal.
%  options       - structure with the optional parameters used for the
%                  numerical inversion algorithm cf2DistGP.
%
% OUTPUTS:
%  result       - structure with the following items:
%                 result.meanXNew,
%                 result.stdXNew,
%                 result.cfXNew,
%                 result.meanYNew,
%                 result.stdYNew,
%                 result.cfYNew,
%                 result.resultYNew,
%                 result.nLP,
%                 result.n,
%                 result.xPoints,
%                 result.xWeights.
%
% EXAMPLE 1
%  % Comparative Polynomial Calibration by PolyCal
%  x    = [4.7, 5.5, 6.5, 7.5, 8.4]';
%  y    = [5.8, 7.9, 10.4, 12.7, 15.2]';
%  xFit = [4.7, 6.5, 8.4]';
%  clear options
%  options.order = 1;
%  options.cfXA = {@(t)cf_Normal(0.01*t), ...
%          @(t)cf_Normal(0.1*t), ...
%          @(t)cf_Normal(0.1*t), ...
%          @(t)cf_Student(0.1*t,5), ...
%          @(t)cf_RectangularSymmetric(0.2*t)};
%  options.cfXB = {@(t)cf_RectangularSymmetric(0.1*t), ...
%          @(t)cf_RectangularSymmetric(0.1*t), ...
%          @(t)cf_TriangularSymmetric(0.1*t), ...
%          @(t)cf_TriangularSymmetric(0.1*t), ...
%          @(t)cf_RectangularSymmetric(0.2*t)};
%  options.cfXB0 = [];
%  options.cfYA = {@(t)cf_Normal(0.15*t), ...
%          @(t)cf_Normal(0.2*t), ...
%          @(t)cf_Normal(0.2*t), ...
%          @(t)cf_Normal(0.3*t), ...
%          @(t)cf_Normal(0.3*t)};
%  options.cfYB = {@(t)cf_RectangularSymmetric(0.2*t), ...
%          @(t)cf_RectangularSymmetric(0.2*t), ...
%          @(t)cf_TriangularSymmetric(0.3*t), ...
%          @(t)cf_TriangularSymmetric(0.3*t), ...
%          @(t)cf_RectangularSymmetric(0.3*t)};
%  options.cfXB0 = @(t)cf_ArcsineSymmetric(0.1*t);
%  options.tolDiff = 1e-2;
%  PolyCalResult = PolyCal(x,y,xFit,options);
%  % Set the inputs for the algorithm PolyCalFit
%  meanXNew = [5.5, 6, 6.5]';
%  cfXNew = { ...
%         @(t)cf_Normal(0.1*t) .* cf_RectangularSymmetric(0.2*t), ...
%         @(t)cf_Student(0.1*t,5) .* cf_RectangularSymmetric(0.2*t), ...
%         @(t)cf_Normal(0.2*t) .* cf_RectangularSymmetric(0.2*t)};
%  result = PolyCalFit(meanXNew,cfXNew,PolyCalResult,options);
%
% REFERENCES
%
% [1] ISO/TS 28037:2010. Determination and Use of Straight-Line Calibration
%     Functions. International StandardsOrganization, Geneva, September
%     (2010).
% [2] ISO/TS 28038:2018. Determination and Use of Polynomial Calibration
%     Functions. International Standards Organization, Geneva, December
%     (2018).
% [3] JCGM 100:2008. Evaluation of measurement data – Guide to the
%     expression of uncertainty in measurement (GUM 1995 with minor
%     corrections), ISO, BIPM, IEC, IFCC, ILAC, IUPAC, IUPAP and OIML,
%     (2008).
% [4] KUBÁČEK L. Foundations of Estimation Theory, Elsevier, Amsterdam
%     (1988).
% [5] WITKOVSKÝ V. Numerical inversion of a characteristic function: An
%     alternative tool to form the probability distribution of output
%     quantity in linear measurement models. ACTA IMEKO 5 (3), (2016)
%     32–44.
% [6] WITKOVSKÝ V. CharFunTool: The Characteristic Functions Toolbox
%     (MATLAB). https://github.com/witkovsky/CharFunTool, (2020).
% [7] WITKOVSKÝ V. and WIMMER G. Generalized polynomial comparative
%     calibration: Parameter estimation and applications. In: Advances in
%     Measurements and Instrumentation: Reviews. S.Y. Yurish (Ed.).
%     Barcelona, Spain, IFSA (International Frequency Sensor Association
%     Publishing) Publishing S.L., (2018) 15-52.
% [8] WITKOVSKÝ V. and WIMMER G. PolyCal - MATLAB algorithm for comparative
%     polynomial calibration and its applications. AMCTM 2020.

% Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 29-Sep-2020 14:03:52
%% Check the Inputs
narginchk(3, 4);

if nargin < 4, options = []; end

if ~isfield(options, 'nLP')
    options.nLP = [];
end

if ~isfield(options, 'isPlot')
    options.isPlot = true;
end

isPlot          = options.isPlot;
options.isPlot  = false;
nLP             = options.nLP;
order           = PolyCalResult.order;
betaEstimate    = PolyCalResult.betaEstimate;
LX              = PolyCalResult.LX;
KY              = PolyCalResult.KY;
cfXA            = PolyCalResult.cfXA;
cfXB            = PolyCalResult.cfXB;
cfXB0           = PolyCalResult.cfXB0;
cfYA            = PolyCalResult.cfYA;
cfYB            = PolyCalResult.cfYB;
cfYB0           = PolyCalResult.cfYB0;

if isempty(nLP)
    nLP = 101;
end

meanXNew = meanXNew(:);
n        = length(meanXNew);

[prob,xWeights] = LegendrePoints(nLP,0,1);
if n == 1
    [~,~,~,points] = cf2DistGP(cfXNew,[],prob,options);
    xPoints = meanXNew + points;
    w = ones(nLP,order+1);
    for k = 2:(order+1)
        w(:,k) = w(:,k-1).* xPoints(:);
    end
    cfy = @(t) 0;
    for j = 1:nLP
        y        = w(j,:)*betaEstimate;
        wX       = w(j,:)*LX;
        wY       = w(j,:)*KY;
        cfXpoint = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,y);
        cfy      = @(t) cfy(t) + xWeights(j) * cfXpoint(t);
    end
    cfYNew = cfy;
    resultYNew = cf2DistGP(cfYNew,[],[],options);
    if isPlot
        figure
        plot(resultYNew.x,resultYNew.pdf);grid
        xlabel('y')
        ylabel('pdf')
        title('State-of-Knowledge PDF of the Response Variable')
        figure
        plot(resultYNew.x,resultYNew.cdf);grid
        xlabel('y')
        ylabel('cdf')
        title('State-of-Knowledge CDF of the Response Variable')
    end
elseif n > 1
    resultYNew = cell(n,1);
    xPoints = zeros(n,nLP);
    for i = 1:n
        [~,~,~,points] = cf2DistGP(cfXNew{i},[],prob,options);
        xPoints(i,:) = meanXNew(i) + points;
    end
    cfYNew = cell(n,1);
    for i = 1:n
        w = ones(nLP,order+1);
        for k = 2:(order+1)
            w(:,k) = w(:,k-1).* xPoints(i,:)';
        end
        cfy = @(t) 0;
        for j = 1:nLP
            y        = w(j,:)*betaEstimate;
            wX       = w(j,:)*LX;
            wY       = w(j,:)*KY;
            cfXpoint = PolyCalCF(wX,wY,cfXA,cfXB,cfXB0,cfYA,cfYB,cfYB0,y);
            cfy      = @(t) cfy(t) + xWeights(j) * cfXpoint(t);
        end
        cfYNew{i} = cfy;
        resultYNew{i} = cf2DistGP(cfYNew{i},[],[],options);
        if isPlot
            figure
            plot(resultYNew{i}.x,resultYNew{i}.pdf);grid
            xlabel('y')
            ylabel('pdf')
            title('State-of-Knowledge PDF of the Response Variable')
            figure
            plot(resultYNew{i}.x,resultYNew{i}.cdf);grid
            xlabel('y')
            ylabel('cdf')
            title('State-of-Knowledge CDF of the Response Variable')
        end
    end
end

stdXNew            = PolyCalCF2Std(cfXNew);
[stdYNew,meanYNew] = PolyCalCF2Std(cfYNew);

result.meanXNew   = meanXNew;
result.stdXNew    = stdXNew;
result.cfXNew     = cfXNew;
result.meanYNew   = meanYNew;
result.stdYNew    = stdYNew;
result.cfYNew     = cfYNew;
result.resultYNew = resultYNew;
result.nLP        = nLP;
result.n          = n;
result.xPoints    = xPoints;
result.xWeights   = xWeights;

end