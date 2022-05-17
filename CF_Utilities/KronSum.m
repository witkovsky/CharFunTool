function [S,id1,id2,Sfull] = KronSum(A,B,options)
%KronSum Kronecker tensor sum.
%   KronSum(X,Y) is the Kronecker tensor SUM of A and B. The result is a
%   large matrix formed by taking all possible sum combinations between the
%   elements of A and those of B. For example, if A is 2 by 3, then
%   KronSum(A,B) is
%
%      [ A(1,1)+B  A(1,2)+B  A(1,3)+B
%        A(2,1)+B  A(2,2)+B  A(2,3)+B ]
%
%   SYNTAX:
%     S = KronSum(A,B)
%     [S,id1,id2,Sfull] = KronSum(A,B,options)
%
%   INPUT:
%     A       - First input matrix of real numbers. Here, inputs must be
%               2-D matrices or vectors. If B = []
%     B       - Second input matrix of real numbers. Here, inputs must be
%               2-D matrices or vectors. If B = [] then by default we set B
%               = A. 
%     options - structure with the following default parameters:
%               - options.isSort = false,
%               - options.isUnique = false.
% 
%   OUTPUT:
%     S       - matrix with Kronecker sums
%     id1     - id vector of indices as a result from the applied sort of
%               unique algorithms
%     id2     - id vector of indices as a result from the applied unique
%               algorithm.
%
%   EXAMPLE 1
%   A = reshape(1:12,3,4);
%   B = (3:7)';
%   S = KronSum(A,B);
%   disp(S)
%
%   EXAMPLE 2
%   A = randn(3,4);
%   B = rand(2,3);
%   options.isUnique = true;
%   S = KronSum(0.1*A,-0.5*B,options);
%   disp(S)
%
%   The code is adapted from the MATLAB function kron.m.

%   Viktor Witkovsky (witkovsky@savba.sk)
%   Ver.: 17-May-2022 15:04:54

%% ALGORITHM
%[S,id1,id2,Sfull] = KronSum(A,B,options)

%% CHECK THE INPUT PARAMETERS
narginchk(1, 3);
if nargin < 2, B = []; end
if nargin < 3, options = []; end

if isempty(B)
    B = A;
end

if ~isfield(options, 'isUnique')
    options.isUnique = false;
end

if ~isfield(options, 'isSort')
    options.isSort = false;
end

if ~ismatrix(A) || ~ismatrix(B)
    error(message('MATLAB:kron:TwoDInput'));
end

id1 = [];
id2 = [];
ma = size(A, 1);
na = size(A, 2);
mb = size(B, 1);
nb = size(B, 2);
A = reshape(A, 1, ma, 1, na);
B = reshape(B, mb, 1, nb, 1);
Sfull = reshape(A + B, ma*mb, na*nb);
S = Sfull;

if options.isSort
    [S,id1] = sort(S(:));
end

if options.isUnique
    [S,id1,id2] = unique(S(:));
end

end