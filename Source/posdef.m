function A = posdef(A, e)
%POSDEF Make symmetric matrix positive definite.
%
%   A = posdef(A, e)
%   Reconstitute a symmetric matrix as a symmetric positive definite
%   matrix. The eigenvectors of the new matrix will be identical and all
%   eigenvalues larger than e will be identical. All eigenvalues less than
%   e will now be equal to e.
%
%   A = posdef(A)
%   The default e is eight times the machine epsilon of the largest
%   eigenvalue.
%
%   This function is most useful in taking matrices that are supposed to be
%   positive definite, but for numerical reasons have eigenvalues that are
%   around zero. This function forces those essentially zero eigenvalues to
%   be the smallest detectable positive value, instead. Algorithms that
%   rely on the matrix being symmetric positive definite will be less
%   likely to crash.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Eigendecompose matrix
[lambda, Q] = infoeig(A);

% Default epsilon is machine epsilon of max eigenvalue
if nargin < 2 || isempty(e)
    e = 8*eps(max(lambda));
end

% Floor epsilon values
lambda(lambda < e) = e;

% Reconstitute matrix
A = symmat(Q*diag(lambda)*Q');
