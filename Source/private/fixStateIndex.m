function ind = fixStateIndex(m, ind)

if is(m, 'Integration')
    % Create a fake model
    [m.States(1:numel(m.x_names)).Name] = deal(m.x_names{:});
end

if ischar(ind) || iscellstr(ind)
    ind = nameToStateIndex(m, ind, false);
end
