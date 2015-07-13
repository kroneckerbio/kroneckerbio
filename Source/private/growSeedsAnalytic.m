function seeds = growSeedsAnalytic(seeds, ns)

if nargin < 2
    ns = 0;
    if nargin < 1
        seeds = [];
    end
end

current = numel(seeds);
add = ns - current;
if add > 0 || isempty(seeds)
    % Double length
    add = max(current,add);
    seeds = [seeds; struct('Name', cell(add,1), 'ID', cell(add,1), 'Value', zeros(1))];
end