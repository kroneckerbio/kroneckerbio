function ind = fixParameterIndex(m, ind)

if is(m, 'Integration')
    % Create a fake model
    [m.Parameters(1:numel(m.k_names)).Name] = deal(m.k_names{:});
end

if ischar(ind) || iscellstr(ind)
    ind = nameToParameterIndex(m, ind, false);
end
