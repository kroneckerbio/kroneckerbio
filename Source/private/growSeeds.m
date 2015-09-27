function seeds = growSeeds(seeds, ns)

if nargin < 2
    ns = 0;
    if nargin < 1
        seeds = [];
    end
end

% Add more room in vector if necessary
current = numel(seeds);
add = ns - current;
if add > 0 || isempty(seeds)
    % Double length
    add = max(current,add);
    seeds = [seeds; struct('Name', cell(add,1), 'Value', cell(add,1))];
end
