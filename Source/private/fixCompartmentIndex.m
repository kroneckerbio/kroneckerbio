function ind = fixCompartmentIndex(m, ind)

if ischar(ind) || iscellstr(ind)
    ind = nameToCompartmentIndex(m, ind, false);
end
