function B = padarray(A, padSize, post)
% 
% function B = padarray(A, padSize)
% 
% This function is component of the full padarray function that merely pads
% the input array A with zeros in the appropriate dimensions, such that
% the final array size is B.

assert(strcmp(post, 'post'), '', 'Error in our implementation of padArray');
nDims = numel(padSize);

ind   = cell(1,nDims);
sizeB = zeros(1,nDims);
for k = 1:nDims
    M = size(A,k);
    ind{k}   = 1:M;
    sizeB(k) = M + padSize(k);
end

B         = zeros(sizeB);
B(ind{:}) = A;