function states = emptySpecies(nx)

if isscalar(nx)
    nx = [nx,1];
end

states = struct('Name', cell(nx), 'Compartment', cell(nx), 'IsInput', cell(nx), 'Value', cell(nx));