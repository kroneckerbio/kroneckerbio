function C = mtimestall(A, B)
% mtimes has a bug in it that will allocate too much memory is A is very
% tall, yet very sparse. This is a workaround that avoids using mtimes and
% instead uses scalar multiplication.

if nnz(A) < size(A,1) && ~(isa(A,'sym') || isa(B,'sym'))
    % Allocate memory for output
    C = sparse([], [], [], size(A,1), size(B,2), nnz(A)*min(nnz(B),size(B,2)));
    
    % Find all nonzeros of B
    [i, j, s] = find(B);
    
    % Loop over nonzeros of B and multiply
    for index = 1:numel(i)
        C(:,j(index)) = C(:,j(index)) + A(:,i(index)) * s(index);
    end
else
    % This is handled by normal mtimes
    C = A * B;
end
