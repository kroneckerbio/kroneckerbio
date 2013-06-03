function seeds = growSeeds(seeds, ns)

% Add more room in vector if necessary
current = numel(seeds);
add = ns - current;
if add > 0 || isempty(seeds)
    % Double length
    add = max(current,add);
    seeds = [seeds; struct('Name', cell(add,1), 'Value', cell(add,1))];
end
