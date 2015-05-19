function D = diff_vectorized(nums, dens)

assert(length(nums) == length(dens), 'nums and dens should be of the same length')

if isempty(nums)
    D = sym(zeros(size(nums)));
    return
end

assert(size(nums,2) == 1 && size(dens,2) == 1, 'nums and dens should be column vectors')

filename = which('diff_vectorized.mu');
read(symengine, filename);
D = feval(symengine,'diff_vectorized',nums,dens);

end