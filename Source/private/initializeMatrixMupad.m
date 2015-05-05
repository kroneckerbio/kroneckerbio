function matout = initializeMatrixMupad(i,j,s,m,n)
% More efficient method for initializing sparse symbolic matrices.
% Utilizes a MuPAD function that generates a MuPAD table,
% then uses the table to initialize the matrix.

assert(length(i) == length(j), 'i and j must be of the same length')
assert(length(i) == length(s), 'i and s must be of the same length')
assert(isscalar(m), 'm must be a scalar')
assert(isscalar(n), 'n must be a scalar')
assert(isempty(i) || max(i) <= m, 'All elements of i must be less than or equal to m')
assert(isempty(j) || max(j) <= n, 'All elements of j must be less than or equal to n')

filename = which('sparse_mupad.mu');
read(symengine, filename);
matout = feval(symengine,'sparse_mupad',i, j, s, m, n);

end