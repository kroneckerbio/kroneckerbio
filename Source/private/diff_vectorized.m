function D = diff_vectorized(nums, dens, simplifyMethod)

if nargin < 3
    simplifyMethod = '';
end

assert(length(nums) == length(dens), 'nums and dens should be of the same length')

if isempty(nums)
    D = sym(zeros(size(nums)));
    return
end

assert(size(nums,2) == 1 && size(dens,2) == 1, 'nums and dens should be column vectors')

switch simplifyMethod
    case 'simplifyFraction'
        filename = which('diff_vectorized_simplifyFraction.mu');
        procname = 'diff_vectorized_simplifyFraction';
    case 'simplify'
        filename = which('diff_vectorized_simplify.mu');
        procname = 'diff_vectorized_simplify';
    case ''
        filename = which('diff_vectorized.mu');
        procname = 'diff_vectorized';
    otherwise
        warning('Unrecognized expression simplification method. Defaulting to no simplification.')
        filename = which('diff_vectorized.mu');
        procname = 'diff_vectorized';
end
read(symengine, filename);
D = feval(symengine,procname,nums,dens);

end