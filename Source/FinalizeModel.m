function m = FinalizeModel(m)
%FinalizeModel Update the mathematical components of the model to reflect
%   changes made to the model
%
%   m = FinalizeModel(m)
%
%   Some components of Kronecker models are cyclically dependent
%   (compartment volumes depend on species and species are placed into
%   compartments). Because of this, it is not possible to design Kronecker
%   in such a way that the model is ready to use after the addition of any
%   component. Furthermore, rebuilding the mathematical components of
%   Kronecker takes a non-trivial amount of time. Since most users add
%   components in groups, simply storing the components until they are all
%   added saves time. Unfortunately, this requires the user to remember to
%   call this function on the model when he is done.
%
%   Forgetting to call this function and using a model which has Ready set
%   to false will result in undefined behavior.
%
%   Inputs
%   m: [ model struct scalar ]
%
%   Outputs
%   m: [ model struct scalar ]

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
assert(nargin >= 1, 'KroneckerBio:FinalizeModel:TooFewInputs', 'FinalizeModel requires at least 1 input argument')
assert(isscalar(m), 'KroneckerBio:FinalizeModel:MoreThanOneModel', 'The model structure must be scalar')

%% Place compartments
nvNew = m.add.nv;

% Check if item by this name already exists
for iv = 1:nvNew
    name = m.add.Compartments(iv).Name;
    
    if any(strcmp(name, [{m.Compartments.Name}, {m.add.Compartments(1:iv-1).Name}]))
        error('KroneckerBio:FinalizeModel:RepeatCompartment', 'There is already a compartment with the name %s (new compartment #%i)', name, iv)
    end
end

% Append new items
m.Compartments = [m.Compartments; m.add.Compartments(1:nvNew)];

% Update count
nv = numel(m.Compartments);
m.nv = nv;
vNames = vec({m.Compartments.Name});

%% Place seeds
nsNew = m.add.ns;

% Check if item by this name already exists
for is = 1:nsNew
    name = m.add.Seeds(is).Name;
    
    if any(strcmp(name, [{m.Seeds.Name}, {m.add.Seeds(1:is-1).Name}]))
        error('KroneckerBio:FinalizeModel:RepeatSeed', 'There is already a seed with the name %s (new seed #%i)', name, is)
    end
end

% Append new items
m.Seeds = [m.Seeds; m.add.Seeds(1:nsNew)];

% Update count
ns = numel(m.Seeds);
m.ns = ns;
sNames = {m.Seeds.Name};

%% Place inputs
nuNew = m.add.nu;

% Assemble full names of states
existing_full_names = strcat(vec({m.Inputs.Compartment}), '.', vec({m.Inputs.Name}));
new_full_names = strcat(vec({m.add.Inputs.Compartment}), '.', vec({m.add.Inputs.Name}));

% Check if a state by this name already exists
for iu = 1:nuNew
    compartment = m.add.Inputs(iu).Compartment;
    name = m.add.Inputs(iu).Name;
    full_name = [compartment '.' name];
    
    if any(strcmp(full_name, [existing_full_names, new_full_names(1:iu-1)]))
        error('KroneckerBio:FinalizeModel:RepeatInput', 'There is already an input with the name %s in compartment %s (new input #%i)', name, compartment, iu)
    end
end

% Append new items
m.Inputs = [m.Inputs; m.add.Inputs(1:nuNew)];

% Update count
nu = numel(m.Inputs);
m.nu = nu;
uNamesFull = [existing_full_names; new_full_names];
m.vuInd = lookup(vec({m.Inputs.Compartment}), vNames);

%% Place states
nxNew = m.add.nx;

% Assemble full names of states
existing_full_names = strcat(vec({m.States.Compartment}), '.', vec({m.States.Name}));
new_full_names = strcat(vec({m.add.States.Compartment}), '.', vec({m.add.States.Name}));

% Check if a state by this name already exists
for ix = 1:nxNew
    compartment = m.add.States(ix).Compartment;
    name = m.add.States(ix).Name;
    full_name = [compartment '.' name];
    
    if any(strcmp(full_name, [existing_full_names; new_full_names(1:ix-1)]))
        error('KroneckerBio:FinalizeModel:RepeatState', 'There is already a state with the name %s in compartment %s (new state #%i)', name, compartment, ix)
    end
end

% Append new item
m.States = [m.States; m.add.States(1:nxNew)];

% Update count
nx = numel(m.States);
xNamesFull = [existing_full_names; new_full_names];
m.nx = nx;
m.vxInd = lookup(vec({m.States.Compartment}), vNames);

%% Place outputs
nyNew = m.add.ny;

% Check if item by this name already exists
for iy = 1:nyNew
    name = m.add.Outputs(iy).Name;
    
    if any(strcmp(name, [{m.Outputs.Name}, {m.add.Outputs(1:iy-1).Name}]))
        error('KroneckerBio:FinalizeModel:RepeatOutput', 'There is already an output with the name %s (new output #%i)', name, iy)
    end
end

% Append new items
m.Outputs = [m.Outputs; m.add.Outputs(1:nyNew)];

% Update count
ny = numel(m.Outputs);
m.ny = ny;

%% Place parameters
nkNew = m.add.nk;

% Check if item by this name already exists
for ik = 1:nkNew
    name = m.add.Parameters(ik).Name;
    
    if any(strcmp(name, [{m.Parameters.Name}, {m.add.Parameters(1:ik-1).Name}]))
        error('KroneckerBio:FinalizeModel:RepeatParameter', 'There is already a parameter with the name %s (new parameter #%i)', name, ik)
    end
end

% Append new items
m.Parameters = [m.Parameters; m.add.Parameters(1:nkNew)];

% Update count
nk = numel(m.Parameters);
m.nk = nk;
kNames = vec({m.Parameters.Name});

%% Place reactions
xuNamesFull = [uNamesFull; xNamesFull];

nrNew = m.add.nr; % Number of reaction specifications
new_reactions = emptyReactions(nrNew);

% Loop over added reaction specifications
for ir = 1:nrNew
    % Possible compartments as a cell vector of strings
    if isempty(m.add.Reactions(ir).Compartment)
        possible_comp = vNames;
    else
        possible_comp = m.add.Reactions(ir).Compartment;
    end
    
    nReac = numel(m.add.Reactions(ir).Reactants);
    nProd = numel(m.add.Reactions(ir).Products);
    
    % Find full name for each reactant
    full_reactant_names = cell(nReac,1);
    for iReac = 1:nReac
        incomplete_name = m.add.Reactions(ir).Reactants{iReac};
        if any(incomplete_name == '.')
            assert(nnz(strcmp(incomplete_name, xuNamesFull)) == 1, 'KroneckerBio:FinalizeModel:ReactantNotFound', 'Reaction %s (#%i) has reactant %s (#%i), which was not found as a species', m.add.Reactions(ir).Name, ir, incomplete_name, iReac)

            full_reactant_names{iReac} = incomplete_name;
        else
            potential_full_names = strcat(possible_comp, '.', incomplete_name);
            species_index = lookup(xuNamesFull, potential_full_names);
            
            if nnz(species_index) < 1
                error('KroneckerBio:FinalizeModel:ReactantNotFound', 'Reaction %s (#%i) has reactant %s (#%i), which was not found as a species', m.add.Reactions(ir).Name, ir, incomplete_name, iReac)
            elseif nnz(species_index) > 1
                error('KroneckerBio:FinalizeModel:AmbiguousReactant', 'Reaction %s (#%i) has reactant %s (#%i), which is a species in several compartments so it is ambiguous as to which one should be applied', m.add.Reactions(ir).Name, ir, incomplete_name, iReac)
            end
            
            full_reactant_names{iReac} = xuNamesFull{logical(species_index)};
        end
    end

    % Find full name for each product
    full_product_names = cell(nProd,1);
    for iProd = 1:nProd
        incomplete_name = m.add.Reactions(ir).Products{iProd};
        if any(incomplete_name == '.')
            assert(nnz(strcmp(incomplete_name, xuNamesFull)) == 1, 'KroneckerBio:FinalizeModel:ProductNotFound', 'Reaction %s (#%i) has product %s (#%i), which was not found as a species', m.add.Reactions(ir).Name, ir, incomplete_name, iProd)

            full_product_names{iProd} = incomplete_name;
        else
            potential_full_names = strcat(possible_comp, '.', incomplete_name);
            species_index = lookup(xuNamesFull, potential_full_names);
            
            if nnz(species_index) < 1
                error('KroneckerBio:FinalizeModel:ProductNotFound', 'Reaction %s (#%i) has product %s (#%i), which was not found as a species', m.add.Reactions(ir).Name, ir, incomplete_name, iProd)
            elseif nnz(species_index) > 1
                error('KroneckerBio:FinalizeModel:AmbiguousProduct', 'Reaction %s (#%i) has product %s (#%i), which is a species in several compartments so it is ambiguous as to which one should be applied', m.add.Reactions(ir).Name, ir, incomplete_name, iProd)
            end
            
            full_product_names{iProd} = xuNamesFull{logical(species_index)};
        end
    end
    
    % Add new reaction to model
    new_reactions(ir).Name = m.add.Reactions(ir).Name;
    new_reactions(ir).Reactants = full_reactant_names;
    new_reactions(ir).Products = full_product_names;
    new_reactions(ir).Parameter = m.add.Reactions(ir).Parameter;
end

% Append new items
m.Reactions = [m.Reactions; new_reactions];

% Update count
nr = numel(m.Reactions);
m.nr = nr;

%% Process compartments
% Dimensions
m.dv = vec([m.Compartments.Dimension]);

% Size
m.v = vec([m.Compartments.Size]);

%% Process seeds
m.s = vec([m.Seeds.Value]);

%% Process states
% Set up map from seeds to initial conditions
state_seed_values = vec({m.States.InitialValue});

constant_state_IC = zeros(nx,1);

ndx0dsEntries = 0;
dx0dsEntries = zeros(nx,2);
dx0dsValues = zeros(nx,1);

for ix = 1:nx
    seed_value = state_seed_values{ix};
    
    % Constant entry first
    constant_entry = strcmp('', seed_value(:,1));
    if any(constant_entry)
        constant_state_IC(ix) = seed_value{constant_entry,2};
        seed_value(constant_entry,:) = [];
    end
    
    % Add seed entries
    for i = 1:size(seed_value,1)
        seed_index = find(strcmp(seed_value{i,1}, sNames));
        assert(numel(seed_index) == 1, 'KroneckerBio:FinalizeModel:InvalidSeedName', ['Species ' xNamesFull{ix} ' has an invalid seed ' seed_value])
        
        ndx0dsEntries = ndx0dsEntries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(dx0dsEntries,1);
        if ndx0dsEntries > currentLength
            addlength = max(currentLength, 1);
            dx0dsEntries = [dx0dsEntries; zeros(addlength,2)];
            dx0dsValues  = [dx0dsValues;  zeros(addlength,1)];
        end
        
        dx0dsEntries(ndx0dsEntries,:) = [ix, seed_index];
        dx0dsValues(ndx0dsEntries) = 1;
    end
end

% Assemble seed to state map
dx0dsEntries = dx0dsEntries(1:ndx0dsEntries,:);
dx0dsValues = dx0dsValues(1:ndx0dsEntries);

m.dx0ds = sparse(dx0dsEntries(:,1), dx0dsEntries(:,2), dx0dsValues, m.nx, m.ns);
m.x0c = constant_state_IC;

%% Process inputs
% Condense time varying inputs into u(t) function
m.u = completeInputFunction(m.Inputs);

% Put input parameters into q vector
m.nqu = zeros(nu,1);
for iu = 1:nu
    m.nqu(iu) = numel(m.Inputs(iu).Parameters);
end
m.q = cat(1, m.Inputs.Parameters, zeros(0,1));
nq = numel(m.q);
m.nq = nq;

%% Process outputs
% Entries in each sparse matrix for compartment conversion
nC1Entries = 0;
nC2Entries = 0;
ncEntries = 0;

C1Entries = zeros(0,2);
C1Values  = zeros(0,1);
C2Entries = zeros(0,2);
C2Values  = zeros(0,1);
cEntries  = zeros(0,2);
cValues   = zeros(0,1);

for iy = 1:ny
    nExpr = size(m.Outputs(iy).Expressions,1);
    for iExpr = 1:nExpr
        % Find states that match the expression
        match = find(~cellfun(@isempty, regexp(xNamesFull, m.Outputs(iy).Expressions{iExpr,1}, 'once')));
        nAdd = numel(match);
        nC1Entries = nC1Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(C1Entries,1);
        if nC1Entries > currentLength
            addlength = max(currentLength, nAdd);
            C1Entries = [C1Entries; zeros(addlength,2)];
            C1Values  = [C1Values;  zeros(addlength,1)];
        end
        
        % Add entries
        C1Entries(nC1Entries-nAdd+1:nC1Entries,1) = iy;
        C1Entries(nC1Entries-nAdd+1:nC1Entries,2) = match;
        C1Values(nC1Entries-nAdd+1:nC1Entries) = m.Outputs(iy).Expressions{iExpr,2};
        
        % Find inputs that match the expression
        match = find(~cellfun(@isempty, regexp(uNamesFull, m.Outputs(iy).Expressions{iExpr,1}, 'once')));
        nAdd = numel(match);
        nC2Entries = nC2Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(C2Entries,1);
        if nC2Entries > currentLength
            addlength = max(currentLength, nAdd);
            C2Entries = [C2Entries; zeros(addlength,2)];
            C2Values  = [C2Values; zeros(addlength,1)];
        end
        
        % Add entries
        C2Entries(nC2Entries-nAdd+1:nC2Entries,1) = iy;
        C2Entries(nC2Entries-nAdd+1:nC2Entries,2) = match;
        C2Values(nC2Entries-nAdd+1:nC2Entries) = m.Outputs(iy).Expressions{iExpr,2};

        % Find empty expressions, which are constants
        if isempty(m.Outputs(iy).Expressions{iExpr,1})
            ncEntries = ncEntries + 1;
            
            % Add more room in vector if necessary
            currentLength = size(cEntries,1);
            if ncEntries > currentLength
                addlength = max(currentLength, 1);
                cEntries = [cEntries; zeros(addlength,2)];
                cValues = [cValues; zeros(addlength,1)];
            end
            
            % Add entries
            cEntries(ncEntries,1) = iy;
            cEntries(ncEntries,2) = 1;
            cValues(ncEntries) = m.Outputs(iy).Expressions{iExpr,2};
        end
    end
end

% Remove duplicate entries
[C1Entries, ind] = unique(C1Entries(1:nC1Entries,:), 'rows');
C1Values = C1Values(ind);

[C2Entries, ind] = unique(C2Entries(1:nC2Entries,:), 'rows');
C2Values = C2Values(ind);

[cEntries, ind] = unique(cEntries(1:ncEntries,:), 'rows');
cValues = cValues(ind);

% Construct matrices
m.C1 = sparse(C1Entries(:,1), C1Entries(:,2), C1Values, m.ny, m.nx);
m.C2 = sparse(C2Entries(:,1), C2Entries(:,2), C2Values, m.ny, m.nu);
m.c  = sparse(cEntries(:,1),  cEntries(:,2),  cValues,  m.ny, 1);

%% Process parameters
% Put rate parameter into k vector
m.k = vec([m.Parameters.Value]);

%% Process Reactions
m.rOrder = zeros(nr,1);
m.krInd  = zeros(nr,1);

nSEntries  = 0;
nD1Entries = 0;
nD2Entries = 0;
nD3Entries = 0;
nD4Entries = 0;
nD5Entries = 0;
nD6Entries = 0;
ndEntries  = 0;

SEntries  = zeros(0,2);
SValues   = zeros(0,1);
dD1dkEntries = zeros(0,2);
dD1dkValues  = zeros(0,1);
dD2dkEntries = zeros(0,2);
dD2dkValues  = zeros(0,1);
dD3dkEntries = zeros(0,2);
dD3dkValues  = zeros(0,1);
dD4dkEntries = zeros(0,2);
dD4dkValues  = zeros(0,1);
dD5dkEntries = zeros(0,2);
dD5dkValues  = zeros(0,1);
dD6dkEntries = zeros(0,2);
dD6dkValues  = zeros(0,1);
dddkEntries  = zeros(0,2);
dddkValues   = zeros(0,1);

for ir =1:nr
    % Determine reaction type and species indexes
    nReac = numel(m.Reactions(ir).Reactants);
    nProd = numel(m.Reactions(ir).Products);
    reactants          = zeros(nReac,1);
    reactantsExist     = false(nReac,1);
    reactantsAreStates = false(nReac,1);
    products           = zeros(nProd,1);
    productsAreStates  = false(nProd,1);
    for iReac = 1:nReac
        if isempty(m.Reactions(ir).Reactants{iReac})
            reactantsExist(iReac) = false;
        else
            reactantsExist(iReac) = true;
            reactantsAreStates(iReac) = true;
            temp = find(strcmp(m.Reactions(ir).Reactants(iReac), xNamesFull), 1);
            if ~isempty(temp)
                % Store index to state species
                reactants(iReac) = temp;
                reactantsAreStates(iReac) = true;
            else
                % Store index to input species
                reactants(iReac) = find(strcmp(m.Reactions(ir).Reactants(iReac), uNamesFull), 1);
                reactantsAreStates(iReac) = false;
            end
        end
    end
    
    m.rOrder(ir) = nnz(reactantsExist); % Store order
    
    for iProd = 1:nProd
        temp = find(strcmp(m.Reactions(ir).Products{iProd}, xNamesFull), 1);
        if ~isempty(temp)
            % Store index to product only if it is a state
            products(iProd) = temp;
            productsAreStates(iProd) = true;
        else%it's an input
            % Redundant assignment
            products(iProd) = 0;
            productsAreStates(iProd) = false;
        end
    end
    
    % Extract parameter
    parameter = find(strcmp(m.Reactions(ir).Parameter{1}, kNames));
    assert(~isempty(parameter), 'KroneckerBio:FinalizeModel:MissingReactionParameter', 'Reaction %s (#%i) requires parameter %s, but no parameter by that name was found', m.Reactions(ir).Name, ir, m.Reactions(ir).Parameter{1})
    m.krInd(ir) = parameter; % Store index
    
    % The multiplier that is applied to the rate parameter
    modifier = m.Reactions(ir).Parameter{2};
    
    % Switch on reactant state
    if m.rOrder(ir) == 1 && reactantsAreStates(1)
        % D1/A1 reaction
        % Add S entries
        nAdd = 1 + nnz(productsAreStates);
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Subtract reactant
        SEntries(nSEntries-nAdd+1,1) = reactants(1);
        SEntries(nSEntries-nAdd+1,2) = ir;
        SValues(nSEntries-nAdd+1)    = -1;
        
        % Add products
        productsAddedSoFar = cumsum(productsAreStates);
        for iProd = 1:nProd
            if productsAreStates(iProd)
                Srow = nSEntries-nAdd+1+productsAddedSoFar(iProd);
                SEntries(Srow,1) = products(iProd);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            end
        end
        
        % Add D1 entry
        nD1Entries = nD1Entries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(dD1dkEntries,1);
        if nD1Entries > currentLength
            addlength = max(currentLength, 1);
            dD1dkEntries = [dD1dkEntries; zeros(addlength,2)];
            dD1dkValues  = [dD1dkValues;  zeros(addlength,1)];
        end
        
        dD1dkEntries(nD1Entries,1) = sub2ind([nr,nx], ir, reactants(1));
        dD1dkEntries(nD1Entries,2) = parameter;
        dD1dkValues(nD1Entries)    = modifier;

    elseif m.rOrder(ir) == 2 && reactantsAreStates(1) && reactantsAreStates(2)
        % D2/A2 reaction
        % Order reactants so that freest species is second
        if m.dv(m.vxInd(reactants(1))) > m.dv(m.vxInd(reactants(2))) || (m.dv(m.vxInd(reactants(1))) == m.dv(m.vxInd(reactants(2))) && reactants(1) > reactants(2))
            [reactants(1), reactants(2)] = deal(reactants(2), reactants(1)); % Swap
        end
        
        % Add S entries
        nAdd = 2 + nnz(productsAreStates);
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Subtract reactant 1
        SEntries(nSEntries-nAdd+1,1) = reactants(1);
        SEntries(nSEntries-nAdd+1,2) = ir;
        SValues(nSEntries-nAdd+1)    = -1;
        
        % Subtract reactant 2
        SEntries(nSEntries-nAdd+2,1) = reactants(2);
        SEntries(nSEntries-nAdd+2,2) = ir;
        SValues(nSEntries-nAdd+2)    = -1;

        % Add products
        productsAddedSoFar = cumsum(productsAreStates);
        for iProd = 1:nProd
            if productsAreStates(iProd)
                Srow = nSEntries-nAdd+2+productsAddedSoFar(iProd);
                SEntries(Srow,1) = products(iProd);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            end
        end
        
        % Add D2 entry
        nD2Entries = nD2Entries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(dD2dkEntries,1);
        if nD2Entries > currentLength
            addlength = max(currentLength, 1);
            dD2dkEntries = [dD2dkEntries; zeros(addlength,2)];
            dD2dkValues  = [dD2dkValues;  zeros(addlength,1)];
        end
        
        dD2dkEntries(nD2Entries,1) = sub2ind([nr,nx,nx], ir, reactants(2), reactants(1));
        dD2dkEntries(nD2Entries,2) = parameter;
        dD2dkValues(nD2Entries)    = modifier;
        
    elseif m.rOrder(ir) == 2 && xor(reactantsAreStates(1), reactantsAreStates(2))
        % Order reactants so that freest species is second
        if reactantsAreStates(1)
            if m.dv(m.vxInd(reactants(1))) > m.dv(m.vuInd(reactants(2)))
                [reactants(1), reactants(2)] = deal(reactants(2), reactants(1)); % Swap
                [reactantsAreStates(1), reactantsAreStates(2)] = deal(reactantsAreStates(2), reactantsAreStates(1));
            end
        else%reactantsAreStates(2)
            if m.dv(m.vuInd(reactants(1))) >= m.dv(m.vxInd(reactants(2)))
                [reactants(1), reactants(2)] = deal(reactants(2), reactants(1)); % Swap
                [reactantsAreStates(1), reactantsAreStates(2)] = deal(reactantsAreStates(2), reactantsAreStates(1));
            end
        end
        
        if reactantsAreStates(2)
            % D3/A3 reaction
            % Add S entries
            nAdd = 1 + nnz(productsAreStates);
            nSEntries = nSEntries + nAdd;
            
            % Add more room in vector if necessary
            currentLength = size(nSEntries,1);
            if nSEntries > currentLength
                addlength = max(currentLength, nAdd);
                SEntries = [SEntries; zeros(addlength,2)];
                SValues  = [SValues;  zeros(addlength,1)];
            end
            
            % Subtract reactant 2
            SEntries(nSEntries-nAdd+1,1) = reactants(2);
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = -1;
            
            % Add products
            productsAddedSoFar = cumsum(productsAreStates);
            for iProd = 1:nProd
                if productsAreStates(iProd)
                    Srow = nSEntries-nAdd+1+productsAddedSoFar(iProd);
                    SEntries(Srow,1) = products(iProd);
                    SEntries(Srow,2) = ir;
                    SValues(Srow)    = 1;
                end
            end
        
            % Add D3 entry
            nD3Entries = nD3Entries + 1;
            
            % Add more room in vector if necessary
            currentLength = size(dD3dkEntries,1);
            if nD3Entries > currentLength
                addlength = max(currentLength, 1);
                dD3dkEntries = [dD3dkEntries; zeros(addlength,2)];
                dD3dkValues  = [dD3dkValues;  zeros(addlength,1)];
            end
            
            dD3dkEntries(nD3Entries,1) = sub2ind([nr,nx,nu], ir, reactants(2), reactants(1));
            dD3dkEntries(nD3Entries,2) = parameter;
            dD3dkValues(nD3Entries)    = modifier;
            
        else%reactantsAreStates(1)
            % D4/A4 reaction
            % Add S entries
            nAdd = 1 + nnz(productsAreStates);
            nSEntries = nSEntries + nAdd;
            
            % Add more room in vector if necessary
            currentLength = size(nSEntries,1);
            if nSEntries > currentLength
                addlength = max(currentLength, nAdd);
                SEntries = [SEntries; zeros(addlength,2)];
                SValues  = [SValues;  zeros(addlength,1)];
            end
            
            % Subtract reactant 1
            SEntries(nSEntries-nAdd+1,1) = reactants(1);
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = -1;
            
            % Add products
            productsAddedSoFar = cumsum(productsAreStates);
            for iProd = 1:nProd
                if productsAreStates(iProd)
                    Srow = nSEntries-nAdd+1+productsAddedSoFar(iProd);
                    SEntries(Srow,1) = products(iProd);
                    SEntries(Srow,2) = ir;
                    SValues(Srow)    = 1;
                end
            end
        
            % Add D4 entry
            nD4Entries = nD4Entries + 1;
        
            % Add more room in vector if necessary
            currentLength = size(dD4dkEntries,1);
            if nD4Entries > currentLength
                addlength = max(currentLength, 1);
                dD4dkEntries = [dD4dkEntries; zeros(addlength,2)];
                dD4dkValues  = [dD4dkValues;  zeros(addlength,1)];
            end
        
            dD4dkEntries(nD4Entries,1) = sub2ind([nr,nu,nx], ir, reactants(2), reactants(1));
            dD4dkEntries(nD4Entries,2) = parameter;
            dD4dkValues(nD4Entries)    = modifier;
            
        end
        
    elseif m.rOrder(ir) == 2 && ~reactantsAreStates(1) && ~reactantsAreStates(2)
        % D5/A5 reaction
        % Order reactants so that freest species is second
        if m.dv(m.vuInd(reactants(1))) > m.dv(m.vuInd(reactants(2))) || (m.dv(m.vuInd(reactants(1))) == m.dv(m.vuInd(reactants(2))) && reactants(1) > reactants(2))
            [reactants(1), reactants(2)] = deal(reactants(2), reactants(1)); % Swap
        end
        
        % Add S entries
        nAdd = nnz(productsAreStates);
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Add products
        productsAddedSoFar = cumsum(productsAreStates);
        for iProd = 1:nProd
            if productsAreStates(iProd)
                Srow = nSEntries-nAdd+productsAddedSoFar(iProd);
                SEntries(Srow,1) = products(iProd);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            end
        end
        
        % Add D5 entry
        nD5Entries = nD5Entries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, 1);
            dD5dkEntries = [dD5dkEntries; zeros(addlength,2)];
            dD5dkValues  = [dD5dkValues;  zeros(addlength,1)];
        end
        
        dD5dkEntries(nD5Entries,1) = sub2ind([nr,nu,nu], ir, reactants(2), reactants(1));
        dD5dkEntries(nD5Entries,2) = parameter;
        dD5dkValues(nD5Entries)    = modifier;
        
    elseif m.rOrder(ir) == 1 && ~reactantsAreStates(1)
        % D6/A6 reaction
        % Add S entries
        nAdd = nnz(productsAreStates);
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Add products
        productsAddedSoFar = cumsum(productsAreStates);
        for iProd = 1:nProd
            if productsAreStates(iProd)
                Srow = nSEntries-nAdd+productsAddedSoFar(iProd);
                SEntries(Srow,1) = products(iProd);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            end
        end
        
        % Add D6 entry
        nD6Entries = nD6Entries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(dD1dkEntries,1);
        if nD1Entries > currentLength
            addlength = max(currentLength, 1);
            dD6dkEntries = [dD6dkEntries; zeros(addlength,2)];
            dD6dkValues  = [dD6dkValues;  zeros(addlength,1)];
        end
        
        dD6dkEntries(nD6Entries,1) = sub2ind([nr,nu], ir, reactants(1));
        dD6dkEntries(nD6Entries,2) = parameter;
        dD6dkValues(nD6Entries)    = modifier;
        
    elseif m.rOrder(ir) == 0
        % d/a reaction
        % Add S entries
        nAdd = nnz(productsAreStates);
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Add products
        productsAddedSoFar = cumsum(productsAreStates);
        for iProd = 1:nProd
            if productsAreStates(iProd)
                Srow = nSEntries-nAdd+productsAddedSoFar(iProd);
                SEntries(Srow,1) = products(iProd);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            end
        end
        
        % Add d entry
        ndEntries = ndEntries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(dddkEntries,1);
        if ndEntries > currentLength
            addlength = max(currentLength, 1);
            dddkEntries = [dddkEntries; zeros(addlength,2)];
            dddkValues  = [dddkValues;  zeros(addlength,1)];
        end
        
        dddkEntries(ndEntries,1) = ir;
        dddkEntries(ndEntries,2) = parameter;
        dddkValues(ndEntries)    = modifier;
        
    else
        error('Unknown error occured while processing reactions')
    end
end

%% Construct matrices of model
% Construct stoichiometry matrix
m.S  = sparse(SEntries(1:nSEntries,1), SEntries(1:nSEntries,2), SValues(1:nSEntries), nx, nr);

% Construct D matrices from the sparse entries
m.dD1dk = sparse(dD1dkEntries(1:nD1Entries,1), dD1dkEntries(1:nD1Entries,2), dD1dkValues(1:nD1Entries), nr*nx,    nk);
m.dD2dk = sparse(dD2dkEntries(1:nD2Entries,1), dD2dkEntries(1:nD2Entries,2), dD2dkValues(1:nD2Entries), nr*nx*nx, nk);
m.dD3dk = sparse(dD3dkEntries(1:nD3Entries,1), dD3dkEntries(1:nD3Entries,2), dD3dkValues(1:nD3Entries), nr*nu*nx, nk);
m.dD4dk = sparse(dD4dkEntries(1:nD4Entries,1), dD4dkEntries(1:nD4Entries,2), dD4dkValues(1:nD4Entries), nr*nx*nu, nk);
m.dD5dk = sparse(dD5dkEntries(1:nD5Entries,1), dD5dkEntries(1:nD5Entries,2), dD5dkValues(1:nD5Entries), nr*nu*nu, nk);
m.dD6dk = sparse(dD6dkEntries(1:nD6Entries,1), dD6dkEntries(1:nD6Entries,2), dD6dkValues(1:nD6Entries), nr*nu,    nk);
m.dddk  = sparse(dddkEntries(1:ndEntries,1),   dddkEntries(1:ndEntries,2),   dddkValues(1:ndEntries),   nr,       nk);

% Compute A matrices
m.dA1dk = reshape(m.S * reshape(m.dD1dk, nr, nx*nk   ), nx*nx,    nk);
m.dA2dk = reshape(m.S * reshape(m.dD2dk, nr, nx*nx*nk), nx*nx*nx, nk);
m.dA3dk = reshape(m.S * reshape(m.dD3dk, nr, nu*nx*nk), nx*nu*nx, nk);
m.dA4dk = reshape(m.S * reshape(m.dD4dk, nr, nx*nu*nk), nx*nx*nu, nk);
m.dA5dk = reshape(m.S * reshape(m.dD5dk, nr, nu*nu*nk), nx*nu*nu, nk);
m.dA6dk = reshape(m.S * reshape(m.dD6dk, nr, nu*nk   ), nx*nu,    nk);
m.dadk  = m.S * m.dddk;

% Alternative shape for easier computation of derivatives
m.dA1dk_fk_x  = spermute132(m.dA1dk, [nx,nx,nk],    [nx*nk,nx]);
m.dA2dk_fk_xx = spermute132(m.dA2dk, [nx,nx*nx,nk], [nx*nk,nx*nx]);
m.dA3dk_fk_ux = spermute132(m.dA3dk, [nx,nu*nx,nk], [nx*nk,nu*nx]);
m.dA4dk_fk_xu = spermute132(m.dA4dk, [nx,nx*nu,nk], [nx*nk,nx*nu]);
m.dA5dk_fk_uu = spermute132(m.dA5dk, [nx,nu*nu,nk], [nx*nk,nu*nu]);
m.dA6dk_fk_u  = spermute132(m.dA6dk, [nx,nu,nk],    [nx*nk,nu]);

m.dD1dk_rk_x  = spermute132(m.dD1dk, [nr,nx,nk],    [nr*nk,nx]);
m.dD2dk_rk_xx = spermute132(m.dD2dk, [nr,nx*nx,nk], [nr*nk,nx*nx]);
m.dD3dk_rk_ux = spermute132(m.dD3dk, [nr,nu*nx,nk], [nr*nk,nu*nx]);
m.dD4dk_rk_xu = spermute132(m.dD4dk, [nr,nx*nu,nk], [nr*nk,nx*nu]);
m.dD5dk_rk_uu = spermute132(m.dD5dk, [nr,nu*nu,nk], [nr*nk,nu*nu]);
m.dD6dk_rk_u  = spermute132(m.dD6dk, [nr,nu,nk],    [nr*nk,nu]);

% Construct fullest possible matrices to determine which columns are used
kRand = sparse(rand(nk,1));
m.A1 = reshape(m.dA1dk * kRand, nx,nx);
m.A2 = reshape(m.dA2dk * kRand, nx,nx*nx);
m.A3 = reshape(m.dA3dk * kRand, nx,nu*nx);
m.A4 = reshape(m.dA4dk * kRand, nx,nx*nu);
m.A5 = reshape(m.dA5dk * kRand, nx,nu*nu);
m.A6 = reshape(m.dA6dk * kRand, nx,nu);
m.a  = m.dadk * kRand;

m.D1 = reshape(m.dD1dk * kRand, nr,nx);
%m.D2 = reshape(m.dD2dk * kRand, nr,nx*nx);
m.D2 = reshape(mtimestall(m.dD2dk, kRand), nr,nx*nx);
m.D3 = reshape(m.dD3dk * kRand, nr,nu*nx);
m.D4 = reshape(m.dD4dk * kRand, nr,nx*nu);
m.D5 = reshape(m.dD5dk * kRand, nr,nu*nu);
m.D6 = reshape(m.dD6dk * kRand, nr,nu);
m.d  = m.dddk * kRand;

% Determine used columns of bimolecular matrices
[unused, D2UsedColumns] = find(m.D2);
[unused, D3UsedColumns] = find(m.D3);
[unused, D4UsedColumns] = find(m.D4);
[unused, D5UsedColumns] = find(m.D5);

D2UsedColumns = unique(D2UsedColumns);
D3UsedColumns = unique(D3UsedColumns);
D4UsedColumns = unique(D4UsedColumns);
D5UsedColumns = unique(D5UsedColumns);

[D2UsedSpecies2, D2UsedSpecies1] = ind2sub([nx,nx], D2UsedColumns);
[D3UsedSpecies2, D3UsedSpecies1] = ind2sub([nx,nu], D3UsedColumns);
[D4UsedSpecies2, D4UsedSpecies1] = ind2sub([nu,nx], D4UsedColumns);
[D5UsedSpecies2, D5UsedSpecies1] = ind2sub([nu,nu], D5UsedColumns);

%% Empty out added items
m.add.nv  = 0;
m.add.nu  = 0;
m.add.ns  = 0;
m.add.nx  = 0;
m.add.ny  = 0;
m.add.nk  = 0;
m.add.nr  = 0;

m.add.Compartments = growCompartments([], 0);
m.add.Inputs       = growInputs([], 0);
m.add.Seeds        = growSeeds([], 0);
m.add.States       = growStates([], 0);
m.add.Outputs      = growOutputs([], 0);
m.add.Parameters   = growParameters([], 0);
m.add.Reactions    = growReactions([], 0);

%% Final build of model
m = final(m, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2);

end

function m = final(m, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2)

% Constants
nx = m.nx;
nu = m.nu;
nr = m.nr;

% Build kronecker matrices
sparsek = sparse(m.k);
m.A1 = reshape(m.dA1dk * sparsek, nx,nx);
m.A2 = reshape(m.dA2dk * sparsek, nx,nx*nx);
m.A3 = reshape(m.dA3dk * sparsek, nx,nu*nx);
m.A4 = reshape(m.dA4dk * sparsek, nx,nx*nu);
m.A5 = reshape(m.dA5dk * sparsek, nx,nu*nu);
m.A6 = reshape(m.dA6dk * sparsek, nx,nu);
m.a  = m.dadk * m.k;

m.D1 = reshape(m.dD1dk * sparsek, nr,nx);
%m.D2 = reshape(m.dD2dk * sparsek, nr,nx*nx);
m.D2 = reshape(mtimestall(m.dD2dk, sparsek), nr,nx*nx);
m.D3 = reshape(m.dD3dk * sparsek, nr,nu*nx);
m.D4 = reshape(m.dD4dk * sparsek, nr,nx*nu);
m.D5 = reshape(m.dD5dk * sparsek, nr,nu*nu);
m.D6 = reshape(m.dD6dk * sparsek, nr,nu);
m.d  = reshape(m.dddk  * sparsek, nr,1);

% Handles
m.f = fHidden(m.A1, m.A2, m.A3, m.A4, m.A5, m.A6, m.a, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.dfdx = dfdxHidden(m.A1, m.A2, m.A3, m.A4, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.dfdu = dfduHidden(m.A3, m.A4, m.A5, m.A6, m.v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.dfdk = dfdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.dadk, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.d2fdx2  = d2fdx2Hidden(m.A2, m.v, D2UsedColumns, D2UsedSpecies2, m.vxInd);
m.d2fdu2  = d2fdu2Hidden(m.A5, m.v, D5UsedColumns, D5UsedSpecies2, m.vuInd);
m.d2fdk2  = d2fdk2Hidden(m.nx, m.nk);
m.d2fdudx = d2fdudxHidden(m.A3, m.A4, m.v, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdxdu = d2fdxduHidden(m.A3, m.A4, m.v, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdkdx = d2fdkdxHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdxdk = d2fdxdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdkdu = d2fdkduHidden(m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdudk = d2fdudkHidden(m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.r = rHidden(m.D1, m.D2, m.D3, m.D4, m.D5, m.D6, m.d, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.drdx = drdxHidden(m.D1, m.D2, m.D3, m.D4, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.drdu = drduHidden(m.D3, m.D4, m.D5, m.D6, m.v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.drdk = drdkHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.dD5dk_rk_uu, m.dD6dk_rk_u, m.dddk, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.d2rdx2  = d2rdx2Hidden(m.D2, m.v, D2UsedColumns, D2UsedSpecies2, m.vxInd);
m.d2rdk2  = d2rdk2Hidden(m.nr, m.nk);
m.d2rdxdk = d2rdxdkHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2rdkdx = d2rdkdxHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd, m.nr);

m.Ready = true;
m.Update = @Update;

    function mout = Update(k, s, q)
        % Copy existing model
        mout = m;
        
        % Apply changes
        mout.k = k;
        mout.s = s;
        mout.q = q;
        
        % Distribute values
        if m.nk >= 1
            k = num2cell(k);
            [mout.Parameters.Value] = k{:};
        end
        
        if m.ns >= 1
            sr = num2cell(s(1:numel(m.Seeds)));
            [mout.Seeds.Value] = sr{:};
        end
        
        qIndex = 0;
        for iu = 1:nu
            mout.Inputs(iu).Parameters = q(qIndex+1:qIndex+m.nqu(iu));
            qIndex = qIndex + m.nqu(iu);
        end
        
        % Rebuild model
        mout = final(mout, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2);
    end
end

function handle = fHidden(A1, A2, A3, A4, A5, A6, a, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @f;

    function val = f(t, x, u)
%   f = A1*x + A2*(x kron x/vx) + A3*(u kron x/vx) + A4*(x kron u/vu) + A5*(u kron u/vu) + A6*u + a
        % Compartment column
        xvx = x ./ v(vxInd);
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvx(D2UsedSpecies2), nx*nx,1);
        ukronx = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvx(D3UsedSpecies2), nx*nu,1);
        xkronu = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvu(D4UsedSpecies2), nu*nx,1);
        ukronu = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvu(D5UsedSpecies2), nu*nu,1);
        
        val = A1 * x + A2 * xkronx + A3 * ukronx + A4 * xkronu + A5 * ukronu + A6 * u + a; % f_
    end
end

function handle = dfdxHidden(A1, A2, A3, A4, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdx;

    function val = dfdx(t, x, u)
%   dfdx = A1 + A2*(Ix kron x/vx) + A2*(x kron diag(1/vx)) + A3*(u kron diag(1/vx))
        % Compartment column
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vx(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vx(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nx*nu,nx);
        
        val = A1 + A2 * (Ixkronxvx + xkron1vx) + A3 * ukron1vx + A4 * Ixkronuvu; % f_x
    end
end

function handle = dfduHidden(A3, A4, A5, A6, v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdu;

    function val = dfdu(t, x, u)
%   dfdu = A3*(Iu kron x/vx) + A4*(x kron diag(1/vu)) + A5*(Iu kron u/vu) + A5*(u kron diag(1/vu)) + A6
        % Compartment column
        xvx = x ./ v(vxInd);
        vu = 1 ./ v(vuInd);
        uvu = u .* vu;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(D3UsedColumns, D3UsedSpecies1, xvx(D3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vu(D4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(D5UsedColumns, D5UsedSpecies1, uvu(D5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vu(D5UsedSpecies2), nu*nu,nu);
        
        val = A3 * Iukronxvx + A4 * xkron1vu + A5 * (Iukronuvu + ukron1vu) + A6; % f_u
    end
end

function handle = dfdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, dadk, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @dfdk;

    function val = dfdk(t, x, u)
%   dfdk = dA1dk*x + dA2dk*(x kron x/vx) + dA3dk*(u kron x/vx) + dA4dk*(x kron u/vu) + dA5dk*(u kron u/vu) + dA6dk*u + dadk
        % Compartment column
        xvx = x ./ v(vxInd);
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvx(D2UsedSpecies2), nx*nx,1);
        ukronx = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvx(D3UsedSpecies2), nx*nu,1);
        xkronu = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvu(D4UsedSpecies2), nu*nx,1);
        ukronu = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvu(D5UsedSpecies2), nu*nu,1);
        
        val = dA1dk_fk_x * sparse(x) + dA2dk_fk_xx * xkronx + dA3dk_fk_ux * ukronx + dA4dk_fk_xu * xkronu + dA5dk_fk_uu * ukronu + dA6dk_fk_u * sparse(u); % fk_
        val = reshape(val, nx,nk) + dadk; % f_k
    end
end

function handle = d2fdx2Hidden(A2, v, D2UsedColumns, D2UsedSpecies2, vxInd)
% Constants
nx = numel(vxInd);

% Reverse the internal subscripts in the UsedColumns linear index
[x1ind,x2ind] = ind2sub([nx,nx], D2UsedColumns);
D2UsedColumnsReverse = sub2ind([nx, nx], x2ind, x1ind);

% Return handle
handle = @d2fdx2;

    function val = d2fdx2(t, x, u)
%   d2fdx2 = 2*A2*(Ix kron diag(1/vx))
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        
        % Second-order derivative of (x kron x/vx)
        Ixkron1vx = sparse([D2UsedColumns; D2UsedColumns], [D2UsedColumns; D2UsedColumnsReverse], [vxinv(D2UsedSpecies2); vxinv(D2UsedSpecies2)], nx*nx,nx*nx);
        
        val = A2 * Ixkron1vx; % f_xx
        val = reshape(val, nx*nx, nx); % fx_x
    end
end

function handle = d2fdu2Hidden(A5, v, D5UsedColumns, D5UsedSpecies2, vuInd)
% Constants
nx = size(A5,1);
nu = numel(vuInd);

% Reverse the internal subscripts in the UsedColumns linear index
[u1ind,u2ind] = ind2sub([nu,nu], D5UsedColumns);
D5UsedColumnsReverse = sub2ind([nu, nu], u2ind, u1ind);

% Return handle
handle = @d2fdu2;

    function val = d2fdu2(t, x, u)
%   d2fdu2 = 2*A5*(Iu kron diag(1/vu))
        % Compartment column
        vuinv = 1 ./ v(vuInd);
        
        % Second-order derivative of (u kron u/vu)
        Iukron1vu = sparse([D5UsedColumns; D5UsedColumns], [D5UsedColumns; D5UsedColumnsReverse], [vuinv(D5UsedSpecies2); vuinv(D5UsedSpecies2)], nu*nu,nu*nu);
        
        val = A5 * Iukron1vu; % f_uu
        val = reshape(val, nx*nu, nu); % fu_u
    end
end

function handle = d2fdk2Hidden(nx, nk)
% Return handle
handle = @d2fdk2;

    function val = d2fdk2(t, x, u)
%   d2fdk2 = 0
        val = sparse(nx*nk, nk); % fk_k
    end
end

function handle = d2fdudxHidden(A3, A4, v, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdudx;

    function val = d2fdudx(t, x, u)
%   d2fdudx = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1))
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        
        % Second order derivative of (u kron x/vx) and (x kron u/vu)
        Iukron1vx = sparse(D3UsedColumns, D3UsedColumns, vxinv(D3UsedSpecies2), nx*nu,nx*nu); % xu_xu
        Ixkron1vu = sparse(D4UsedColumns, D4UsedColumns, vuinv(D4UsedSpecies2), nu*nx,nu*nx); % ux_ux
        
        val = A3 * Iukron1vx + spermute132(A4 * Ixkron1vu, [nx,nu,nx],[nx,nx*nu]); % f_xu
        val = reshape(val, nx*nx,nu); % fx_u
    end
end

function handle = d2fdxduHidden(A3, A4, v, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdxdu;

    function val = d2fdxdu(t, x, u)
%   d2fdxdu = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1 - x *{xxx} dvxdx ^ -2)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1 - u *{uxu} dvudu ^ -2))
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        
        % Sparse and non-sparse kronecker multiplication
        Iukron1vx = sparse(D3UsedColumns, D3UsedColumns, vxinv(D3UsedSpecies2), nx*nu,nx*nu); % xu_xu
        Ixkron1vu = sparse(D4UsedColumns, D4UsedColumns, vuinv(D4UsedSpecies2), nu*nx,nu*nx); % ux_ux
        
        val = spermute132(A3 * Iukron1vx, [nx,nx,nu],[nx,nu*nx]) + A4 * Ixkron1vu; % f_ux
        val = reshape(val, nx*nu,nx); % fu_x
    end
end

function handle = d2fdkdxHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @d2fdkdx;

    function val = d2fdkdx(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        xvx = x .* vxinv;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dA1dk_fk_x + dA2dk_fk_xx * (Ixkronxvx + xkron1vx) + dA3dk_fk_ux * ukron1vx + dA4dk_fk_xu * Ixkronuvu; % fx_k
        val = spermute132(val, [nx,nk,nx], [nx*nx,nk]); % fk_x 
    end
end

function handle = d2fdxdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdxdk;

    function val = d2fdxdk(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        xvx = x .* vxinv;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dA1dk_fk_x + dA2dk_fk_xx * (Ixkronxvx + xkron1vx) + dA3dk_fk_ux * ukron1vx + dA4dk_fk_xu * Ixkronuvu; % fx_k
    end
end

function handle = d2fdkduHidden(dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA3dk_fk_ux, 1) / nx;

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        xvx = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvu = u .* vuinv;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(D3UsedColumns, D3UsedSpecies1, xvx(D3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(D5UsedColumns, D5UsedSpecies1, uvu(D5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        
        val = dA3dk_fk_ux * Iukronxvx + dA4dk_fk_xu * xkron1vu + dA5dk_fk_uu * (Iukronuvu + ukron1vu) + dA6dk_fk_u; % fk_u
        val = spermute132(val, [nx,nk,nu], [nx*nu,nk]); % fu_k
    end
end

function handle = d2fdudkHidden(dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        xvx = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvu = u .* vuinv;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(D3UsedColumns, D3UsedSpecies1, xvx(D3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(D5UsedColumns, D5UsedSpecies1, uvu(D5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        
        val = dA3dk_fk_ux * Iukronxvx + dA4dk_fk_xu * xkron1vu + dA5dk_fk_uu * (Iukronuvu + ukron1vu) + dA6dk_fk_u; % fk_u
    end
end

function handle = rHidden(D1, D2, D3, D4, D5, D6, d, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @r;

    function val = r(t, x, u)
%   r = D1*x + D2*(x kron x/vx) + D3*(u kron x/vx) + D4*(x kron u/vu) + D5*(u kron u/vu) + D6*u + d
        % Compartment column
        xvx = x ./ v(vxInd);
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvx(D2UsedSpecies2), nx*nx,1);
        ukronx = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvx(D3UsedSpecies2), nx*nu,1);
        xkronu = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvu(D4UsedSpecies2), nu*nx,1);
        ukronu = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvu(D5UsedSpecies2), nu*nu,1);
        
        val = D1 * x + D2 * xkronx + D3 * ukronx + D4 * xkronu + D5 * ukronu + D6 * u + d; % r_
    end
end

function handle = drdxHidden(D1, D2, D3, D4, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @drdx;

    function val = drdx(t, x, u)
%   drdx = D1 + D2*(Ix kron x/vx) + D2*(x kron diag(1/vx)) + D3*(u kron diag(1/vx))
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        xvx = x .* vxinv;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nx*nu,nx);
        
        val = D1 + D2 * (Ixkronxvx + xkron1vx) + D3 * ukron1vx + D4 * Ixkronuvu; % r_x
    end
end

function handle = drduHidden(D3, D4, D5, D6, v, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @drdu;

    function val = drdu(t, x, u)
%   drdu = D3*(Iu kron x/vx) + D4*(x kron diag(1/vu)) + D5*(Iu kron u/vu) + D5*(u kron diag(1/vu)) + D6
        % Compartment column
        xvx = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvu = u .* vuinv;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(D3UsedColumns, D3UsedSpecies1, xvx(D3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(D5UsedColumns, D5UsedSpecies1, uvu(D5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        
        val = D3 * Iukronxvx + D4 * xkron1vu + D5 * (Iukronuvu + ukron1vu) + D6; % r_u
    end
end

function handle = drdkHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, dD5dk_rk_uu, dD6dk_rk_u, dddk, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nr = size(dddk,1);
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dD1dk_rk_x, 1) / nr;

% Return handle
handle = @drdk;

    function val = drdk(t, x, u)
%   drdk = dD1dk*x + dD2dk*(x kron x/vx) + dD3dk*(u kron x/vx) + dD4dk*(x kron u/vu) + dD5dk*(u kron u/vu) + dD6dk*u + dddk
        % Compartment column
        xvx = x ./ v(vxInd);
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvx(D2UsedSpecies2), nx*nx,1);
        ukronx = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvx(D3UsedSpecies2), nx*nu,1);
        xkronu = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvu(D4UsedSpecies2), nu*nx,1);
        ukronu = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvu(D5UsedSpecies2), nu*nu,1);
        
        val = dD1dk_rk_x * sparse(x) + dD2dk_rk_xx * xkronx + dD3dk_rk_ux * ukronx + dD4dk_rk_xu * xkronu + dD5dk_rk_uu * ukronu + dD6dk_rk_u * sparse(u); % rk_
        val = reshape(val, nr,nk) + dddk; % r_k
    end
end

function handle = d2rdx2Hidden(D2, v, D2UsedColumns, D2UsedSpecies2, vxInd)
% Constants
nr = size(D2,1);
nx = numel(vxInd);

% Reverse the internal subscripts in the UsedColumns linear index
[x1ind,x2ind] = ind2sub([nx,nx], D2UsedColumns);
D2UsedColumnsReverse = sub2ind([nx, nx], x2ind, x1ind);

% Return handle
handle = @d2rdx2;

    function val = d2rdx2(t, x, u)
%   d2rdx2 = 2*D2*(Ix kron diag(1/vx))
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        
        % Second-order derivative of (x kron x/vx)
        Ixkron1vx = sparse([D2UsedColumns; D2UsedColumns], [D2UsedColumns; D2UsedColumnsReverse], [vxinv(D2UsedSpecies2); vxinv(D2UsedSpecies2)], nx*nx,nx*nx);
        
        val = D2 * Ixkron1vx; % r_xx
        val = reshape(val, nr*nx, nx); % rx_x
    end
end

function handle = d2rdk2Hidden(nr, nk)
% Return handle
handle = @d2rdk2;

    function val = d2rdk2(t, x, u)
%   d2rdk2 = 0
        val = sparse(nr*nk, nk); % rk_k
    end
end

function handle = d2rdxdkHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2rdxdk;

    function val = d2rdxdk(t, x, u)
%   d2rdxdk = dD1dk + dD2dk*(Ix kron x/vx) + dD2dk*(x kron diag(1/vx)) + dD3dk*(u kron diag(1/vx)) + dD4dk*(Ix kron u/vu)
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        xvx = x .* vxinv;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dD1dk_rk_x + dD2dk_rk_xx * (Ixkronxvx + xkron1vx) + dD3dk_rk_ux * ukron1vx + dD4dk_rk_xu * Ixkronuvu; % rx_k
    end
end

function handle = d2rdkdxHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd, nr)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dD1dk_rk_x, 1) / nr;

% Return handle
handle = @d2rdkdx;

    function val = d2rdkdx(t, x, u)
%   d2rdxdk = dD1dk + dD2dk*(Ix kron x/vx) + dD2dk*(x kron diag(1/vx)) + dD3dk*(u kron diag(1/vx)) + dD4dk*(Ix kron u/vu)
        % Compartment column
        vxinv = 1 ./ v(vxInd);
        xvx = x .* vxinv;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dD1dk_rk_x + dD2dk_rk_xx * (Ixkronxvx + xkron1vx) + dD3dk_rk_ux * ukron1vx + dD4dk_rk_xu * Ixkronuvu; % rx_k
        val = spermute132(val, [nr,nk,nx], [nr*nx,nk]); % rk_x 
    end
end
