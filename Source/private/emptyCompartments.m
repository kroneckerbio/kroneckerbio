function compartments = emptyCompartments(nv)

compartments = struct('Name', cell(nv,1), 'Dimension', cell(nv,1), 'Size', cell(nv,1));