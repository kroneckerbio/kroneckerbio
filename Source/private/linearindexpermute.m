function indexes = linearindexpermute(indexes, sizes, order)
%linearindexpermute Convert linear indexes so that they point to same
%   elements after a hypotethical permutation
%
%   new_indexes = linearindexpermute(indexes, sizes, order)
%
%   Indexes is some array of linear indexes into an array of
%   size sizes. If the array of size sizes were to be permuted, then the
%   elements to which indexes points will be at new_indexes.

% Extract subscript indexes corresponding to each linear index
[temp{1:numel(order)}] = ind2sub(sizes, indexes);

% Rearrange the indexes according to the permutation
temp = temp(order);

% Convert back to linear indexes
indexes = sub2ind(sizes(order), temp{:});
