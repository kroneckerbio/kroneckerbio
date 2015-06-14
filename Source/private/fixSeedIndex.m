function ind = fixSeedIndex(m, ind)

if is(m, 'Integration')
    % Create a fake model
    [m.Seeds(1:numel(m.s_names)).Name] = deal(m.s_names{:});
end

if ischar(ind) || iscellstr(ind)
    ind = nameToSeedIndex(m, ind, false);
end
