function ind = fixInputIndex(m, ind)

if is(m, 'Integration')
    % Create a fake model
    [m.Inputs(1:numel(m.u_names)).Name] = deal(m.u_names{:});
end

if ischar(ind) || iscellstr(ind)
    ind = nameToInputIndex(m, ind, false);
end
