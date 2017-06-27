function Finv = infoinv(F)
%INFOINV Invert an unstable information matrix.
%
%   Finv = infoinv(F)
%
%   This function inverts a Fisher information matrix into a covariance
%   matrix, but in a way that avoids zeros and infinites corrupting the
%   entire result. Instead, zeros on the diagonal are converted to
%   infinites on the diagonal and infinites on the diagonal are converted
%   to zeros on the diagonal. The rest of the matrix is forced to be
%   positive semidefinite and then inverted. This is done by
%   eigendecomposing the matrix, forcing all small eigenvalues to be above
%   a epsilon threshold, and then inverting.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Desparsify
if issparse(F)
    sparseF = true;
    F = full(F);
else
    sparseF = false;
end

% Size
nT = size(F,1);

% Find specials: zeros and infs
zeroInd = find(diag(F) == 0);
infInd  = find(diag(F) == inf);

% Invert non-special portion of F
invertableInd = true(nT,1);
invertableInd([zeroInd;infInd]) = false;
F = posdef(F(invertableInd,invertableInd));
F = inv(F);

% Create the diagonal of Finv and replace all previous zeros with inf
Finv = zeros(nT,1);
Finv(zeroInd) = inf;
Finv = diag(Finv); % Build the matrix from the diagonal

% Insert all the invertable portions
Finv(invertableInd,invertableInd) = F;

if sparseF
    Finv = sparse(Finv);
end
