function ind = fixOutputIndex(m, ind)

if is(m, 'Integration')
    % Create a fake model
    [m.Outputs(1:numel(m.y_names)).Name] = deal(m.y_names{:});
end

if ischar(ind) || iscellstr(ind)
    ind = nameToOutputIndex(m, ind, false);
end
