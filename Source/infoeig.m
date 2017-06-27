function [lambda, Q] = infoeig(F)
%INFOEIG Eigendecomposition an unstable information matrix.
%
%   [lambda, Q] = infoeig(F)
%
%   Symmetric matrix F is eigendecomposed into vector of eigenvalues lambda
%   and matrix of eigenvectors Q. All eigenvalues are sorted from smallest
%   to largest. F can be indefinite.
%
%   For many applications, it is necessary to decompose a Fisher
%   infomation matrix into its eigenvalues and eigenvectors. However,
%   certain matrices that are valid information matrices cannot be reliably
%   decomposed due to numerical instabilities. These problems include
%   having entirely zero columns and rows representing no information on a
%   parameter, or an infinity on the diagonal representing percect
%   information on a parameter. Both of these lead to singular matrices
%   that cannot be decomposed. This algorithm treats those two cases as
%   special and considers those rows and columns seperately from the rest
%   of the matrix.
%
%   A column and row of all zeros are transformed into an eigenvalue of
%   zero and an eigenvector with a 1 at the index of the column/row.
%
%   An infinity on the diagonal is transformed into an eigenvalue of inf
%   and an eigenvector with a 1 at the index of the infinity. Values on
%   the off-diagonal are irrelavent and are discarded.
%
%   Note that unlike the Matlab function eig, this function always supplies
%   the eigenvalues as the first argument and as a vector.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Size
nT = size(F,1);
diagF = diag(F);

% Find specials: zeros and infs
zeroInd = find(all(F == 0, 2));
n0 = numel(zeroInd);
negInfInd  = find(diagF == -inf);
posInfInd  = find(diagF == inf);
nNegInf = numel(negInfInd);
nPosInf = numel(posInfInd);

% Correct for empty cases
if isempty(zeroInd);   zeroInd   = zeros(0,1); end
if isempty(negInfInd); negInfInd = zeros(0,1); end
if isempty(posInfInd); posInfInd = zeros(0,1); end

% Eig non-special portion of F
invertableInd = true(nT,1);
invertableInd([zeroInd;negInfInd;posInfInd]) = false;

if nargout <= 1
    lambda = eig(F(invertableInd,invertableInd));
else
    [Qplain, lambda] = eig(F(invertableInd,invertableInd));
    lambda = diag(lambda);
end

% Recombine results
positiveLambda = (lambda >= 0);
nNegativeLambda = sum(~positiveLambda);
nPositiveLambda = sum(positiveLambda);
lambda = [inf(nNegInf,1)*-1;
          lambda(~positiveLambda);
          zeros(n0,1);
          lambda(positiveLambda);
          inf(nPosInf,1)];

if nargout >=2
    Q = zeros(nT,nT);
    Q(sub2ind([nT,nT], negInfInd,(1:nNegInf)')) = 1;% -inf
    Q(invertableInd, nNegInf+1:nNegInf+nNegativeLambda) = Qplain(:,~positiveLambda); % negative lambda
    Q(sub2ind([nT,nT], zeroInd,(nNegInf+nNegativeLambda+1:nNegInf+nNegativeLambda+n0)')) = 1; % 0
    Q(invertableInd,  nNegInf+nNegativeLambda+n0+1:nNegInf+nNegativeLambda+n0+nPositiveLambda) = Qplain(:,positiveLambda); % positive lambda
    Q(sub2ind([nT,nT], posInfInd,(nNegInf+nNegativeLambda+n0+nPositiveLambda+1:nT)')) = 1; % +inf
end
