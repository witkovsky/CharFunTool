function [INT,ERR,SubIntervals,Q,Q1] = IntegralGK(fun,Intervals,division,isPlot)
%IntegralGK Evaluates the (G7,K15)-Gauss-Kronod quadrature of function fun
%  over all sub-intervals defined by the matrix Intervals, e.g., Intervals
%  = [-1 0;0 1] specifies integration over the interval [-1,1] subdividet
%  into two subintervals [-1,0] and [0,1]. 
%  The intervals specified by Intervals can be further  subdivided into
%  smaller sub-intervals in relative points specified by the index vector
%  division  (a vector with indices from 1 to 15, e.g. [3 5 8 11 13]). If
%  division = [], there is no further subdivision of the intervals
%  specified in Intervals.
%
% SYNTAX
%   [INT,ERR] = IntegralGK(fun,Intervals)
%   [INT,ERR,SubIntervals,Q,Q1] = IntegralGK(fun,Intervals,division,isPlot)
%
% EXAMPLE
%   [INT,ERR,SubInts] = IntegralGK(@(x) sin(x).*cos(5*x),[-1 0;0 1],[5 8 11])

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 30-Nov-2018 11:14:06


%% FUNCTION
%  [INT,ERR,SubIntervals,Q,Q1] = IntegralGK(fun,Intervals,division,isPlot)

%% CHECK/SET THE INPUT PARAMETERS
narginchk(2, 4);
if nargin < 4, isPlot   = []; end
if nargin < 3, division = []; end

if isempty(isPlot)
    isPlot = false;
end
%% ALGORITHM
[nodeGK,weightGK,WG,G] = GKnodes;
[SubIntervals,t,mids] = GetSubs(Intervals,nodeGK,division);
F   = fun(t);
Q1  = mids.*sum(bsxfun(@times,weightGK,F));
Q2  = mids.*sum(bsxfun(@times,WG,F(G,:)));
Q   = sum(reshape(Q1,length(division)+1,[]));
ERR = sum(abs(Q1-Q2),2);
INT = sum(Q);

if isPlot
    %LW = 'linewidth';
    figure
    plot(t,F,'.-')
    xlabel('x')
    ylabel('function')
    grid on
end

end
%% Function GetSubs
function [SubIntervals,nodes,mids] = GetSubs(Intervals,nodeGK,nodeIdx)
% GETSUBS Auxiliary function. Sub-division of the integration intervals for
%  adaptive Gauss-Kronod quadrature.
%
% EXAMPLE
%  nodeGK =  GKnodes;
%  % Set the indices from 1:15 of nodes used for division of subintervals
%  nodeIdx = [5 8 11];
%  Intervals = [0;1]; % Given as 2 x n  matrix of initial sub-intervals
%  [SubIntervals,nodes,mids] = GetSubs(Intervals,nodeGK,nodeIdx)

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 30-Nov-2018 11:14:06

%% ALGORITHM
mids = 0.5*(Intervals(2,:)-Intervals(1,:));
C = 0.5*(Intervals(2,:)+Intervals(1,:));
nodes = nodeGK*mids + ones(size(nodeGK))*C;

L = [Intervals(1,:); nodes(nodeIdx,:)];
U = [nodes(nodeIdx,:); Intervals(2,:)];
SubIntervals = [reshape(L, 1, []); reshape(U, 1, [])];

mids = 0.5*(SubIntervals(2,:)-SubIntervals(1,:));
C = 0.5*(SubIntervals(2,:)+SubIntervals(1,:));
nodes = nodeGK*mids + ones(size(nodeGK))*C;
end
%% Function GKnodes
function [XK,WK,WG,G] =  GKnodes
%GKNODES Auxiliary function. The (7-15) Gauss-Kronod nodes and weights.

% (c) Viktor Witkovsky (witkovsky@gmail.com)
% Ver.: 30-Nov-2018 11:14:06

%% ALGORITHM
nodes = [ ...
    0.2077849550078984676006894; 0.4058451513773971669066064; ...
    0.5860872354676911302941448; 0.7415311855993944398638648;...
    0.8648644233597690727897128; 0.9491079123427585245261897; ...
    0.9914553711208126392068547];
wt = [ ...
    0.2044329400752988924141620; 0.1903505780647854099132564; ...
    0.1690047266392679028265834; 0.1406532597155259187451896;
    0.1047900103222501838398763; 0.0630920926299785532907007;...
    0.0229353220105292249637320];
wt7 = [0.3818300505051189449503698; ...
    0.2797053914892766679014678; 0.1294849661688696932706114];

XK = [-nodes(end:-1:1); 0; nodes];
WK = [wt(end:-1:1); 0.2094821410847278280129992; wt];
WG = [wt7(end:-1:1); 0.4179591836734693877551020; wt7];
G = (2:2:15)';

end