function m = finalizeModelMassActionAmount(m)
%FinalizeModel Update the mathematical components of the model to reflect
%   changes made to a mass action amount model
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

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Extract model components
nv = numel(m.Compartments);
nk = numel(m.Parameters);
ns = numel(m.Seeds);
nu = numel(m.Inputs);
nx = numel(m.States);
nr = numel(m.Reactions);
nz = numel(m.Rules);
ny = numel(m.Outputs);
nxu = nx + nu;

v_names = vec({m.Compartments.Name});
k_names = vec({m.Parameters.Name});
s_names = vec({m.Seeds.Name});
u_names = vec({m.Inputs.Name});
x_names = vec({m.States.Name});
z_names = vec({m.Rules.Name});
r_names = vec({m.Reactions.Name});
y_names = vec({m.Outputs.Name});
xu_names = [x_names; u_names];

% Make list of all compartment.species in model
u_full_names = vec(strcat({m.Inputs.Compartment}, '.', {m.Inputs.Name}));
x_full_names = vec(strcat({m.States.Compartment}, '.', {m.States.Name}));
xu_full_names = [x_full_names; u_full_names];

% Make logical vector of which species names are unique
unique_xu_names = false(nu+nx,1);
for ixu = 1:nu+nx
    unique_xu_names(ixu) = ~ismember(xu_names{ixu}, [xu_names(1:ixu-1); xu_names(ixu+1:end)]);
end
unique_x_names = unique_xu_names(1:nx);
unique_u_names = unique_xu_names(nx+(1:nu));

%% Warn on detected rules
% Rules aren't currently supported in massaction models
if nz > 0
    warning('KroneckerBio:FinalizeModel:MassActionRule', 'Rules not supported in massaction models. Ignoring.')
end

%% Process compartments
% Dimensions
m.dv = vec([m.Compartments.Dimension]);

% Size
% Entries in each sparse matrix for compartment conversion
nB1Entries = 0;
nB2Entries = 0;
nbEntries = 0;

B1Entries = zeros(0,2);
B1Values  = zeros(0,1);
B2Entries = zeros(0,2);
B2Values  = zeros(0,1);
bEntries  = zeros(0,2);
bValues   = zeros(0,1);

for iv = 1:nv
    if isnumeric(m.Compartments(iv).Size)
        % Handle compartment size given as a number
        nbEntries = nbEntries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(bEntries,1);
        if nbEntries > currentLength
            addlength = max(currentLength, 1);
            bEntries = [bEntries; zeros(addlength,2)];
            bValues = [bValues; zeros(addlength,1)];
        end
        
        % Add entries
        bEntries(nbEntries,1) = iv;
        bEntries(nbEntries,2) = 1;
        bValues(nbEntries) = m.Compartments(iv).Size;
    elseif iscell(m.Compartments(iv).Size)
        % Handle compartment size given as a list of strings
        nExpr = size(m.Compartments(iv).Size,1);
        for iExpr = 1:nExpr
            % Find states that match the expression
            match = find(strcmp(m.Compartments(iv).Size{iExpr,1}, x_full_names));
            nAdd = numel(match);
            nB1Entries = nB1Entries + nAdd;
            
            % Add more room in vector if string
            currentLength = size(B1Entries,1);
            if nB1Entries > currentLength
                addlength = max(currentLength, nAdd);
                B1Entries = [B1Entries; zeros(addlength,2)];
                B1Values  = [B1Values;  zeros(addlength,1)];
            end
            
            % Add entries
            B1Entries(nB1Entries-nAdd+1:nB1Entries,1) = iv;
            B1Entries(nB1Entries-nAdd+1:nB1Entries,2) = match;
            B1Values(nB1Entries-nAdd+1:nB1Entries) = m.Compartments(iv).Size{iExpr,2};
            
            % Find inputs that match the string
            match = find(strcmp(m.Compartments(iv).Size{iExpr,1}, u_full_names));
            nAdd = numel(match);
            nB2Entries = nB2Entries + nAdd;
            
            % Add more room in vector if necessary
            currentLength = size(B2Entries,1);
            if nB2Entries > currentLength
                addlength = max(currentLength, nAdd);
                B2Entries = [B2Entries; zeros(addlength,2)];
                B2Values  = [B2Values; zeros(addlength,1)];
            end
            
            % Add entries
            B2Entries(nB2Entries-nAdd+1:nB2Entries,1) = iv;
            B2Entries(nB2Entries-nAdd+1:nB2Entries,2) = match;
            B2Values(nB2Entries-nAdd+1:nB2Entries) = m.Compartments(iv).Size{iExpr,2};
            
            % Find empty expressions, which are constants
            if isempty(m.Compartments(iv).Size{iExpr,1})
                nbEntries = nbEntries + 1;
                
                % Add more room in vector if necessary
                currentLength = size(bEntries,1);
                if nbEntries > currentLength
                    addlength = max(currentLength, 1);
                    bEntries = [bEntries; zeros(addlength,2)];
                    bValues = [bValues; zeros(addlength,1)];
                end
                
                % Add entries
                bEntries(nbEntries,1) = iv;
                bEntries(nbEntries,2) = 1;
                bValues(nbEntries) = m.Compartments(iv).Size{iExpr,2};
            end
        end
    else
        error('KroneckerBio:Compartment:Size', 'Invalid compartment size type encountered')
    end
end

% Remove duplicate entries
[B1Entries, ind] = unique(B1Entries(1:nB1Entries,:), 'rows');
B1Values = B1Values(ind);

[B2Entries, ind] = unique(B2Entries(1:nB2Entries,:), 'rows');
B2Values = B2Values(ind);

[bEntries, ind] = unique(bEntries(1:nbEntries,:), 'rows');
bValues = bValues(ind);

% Construct matrices
B1 = sparse(B1Entries(:,1), B1Entries(:,2), B1Values, nv, nx);
B2 = sparse(B2Entries(:,1), B2Entries(:,2), B2Values, nv, nu);
b  = sparse(bEntries(:,1),  bEntries(:,2),  bValues,  nv, 1);

% Add to model
m.B1 = B1;
m.B2 = B2;
m.b = b;

m.v       = @v;
m.dvdx    = @dvdx;
m.dvdu    = @dvdu;
m.d2vdx2  = @d2vdx2;
m.d2vdu2  = @d2vdu2;
m.d2vdudx = @d2vdudx;
m.d2vdxdu = @d2vdxdu;

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
        seed_index = find(strcmp(seed_value{i,1}, s_names));
        assert(numel(seed_index) == 1, 'KroneckerBio:FinalizeModel:InvalidSeedName', ['Species ' x_full_names{ix} ' has an invalid seed ' seed_value{i,1}])
        
        ndx0dsEntries = ndx0dsEntries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(dx0dsEntries,1);
        if ndx0dsEntries > currentLength
            addlength = max(currentLength, 1);
            dx0dsEntries = [dx0dsEntries; zeros(addlength,2)];
            dx0dsValues  = [dx0dsValues;  zeros(addlength,1)];
        end
        
        dx0dsEntries(ndx0dsEntries,:) = [ix, seed_index];
        dx0dsValues(ndx0dsEntries)    = seed_value{i, 2};
    end
end

% Assemble seed to state map
dx0dsEntries = dx0dsEntries(1:ndx0dsEntries,:);
dx0dsValues = dx0dsValues(1:ndx0dsEntries);

dx0ds_val = sparse(dx0dsEntries(:,1), dx0dsEntries(:,2), dx0dsValues, m.nx, m.ns);
x0c = constant_state_IC;

m.x0 = @x0;
m.dx0ds = @dx0ds;
m.dx0dk = @dx0dk;
m.d2x0ds2 = @d2x0ds2;
m.d2x0dk2 = @d2x0dk2;
m.d2x0dkds = @d2x0dkds;
m.d2x0dsdk = @d2x0dsdk;

%% Process inputs
% Condense time varying inputs into u(t) function
m.u = vec([m.Inputs.DefaultValue]);

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
    nExpr = size(m.Outputs(iy).Expression,1);
    for iExpr = 1:nExpr
        % Find states that match the string
        match = find(strcmp(m.Outputs(iy).Expression{iExpr,1}, x_full_names));
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
        C1Values(nC1Entries-nAdd+1:nC1Entries) = m.Outputs(iy).Expression{iExpr,2};
        
        % Find inputs that match the string
        match = find(strcmp(m.Outputs(iy).Expression{iExpr,1}, u_full_names));
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
        C2Values(nC2Entries-nAdd+1:nC2Entries) = m.Outputs(iy).Expression{iExpr,2};

        % Find empty expressions, which are constants
        if isempty(m.Outputs(iy).Expression{iExpr,1})
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
            cValues(ncEntries) = m.Outputs(iy).Expression{iExpr,2};
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
C1 = sparse(C1Entries(:,1), C1Entries(:,2), C1Values, m.ny, m.nx);
C2 = sparse(C2Entries(:,1), C2Entries(:,2), C2Values, m.ny, m.nu);
c  = sparse(cEntries(:,1),  cEntries(:,2),  cValues,  m.ny, 1);

% Add to model
m.C1 = C1;
m.C2 = C2;
m.c  = c;

m.y       = @y;
m.dydx    = @dydx;
m.dydu    = @dydu;
m.dydk    = @dydk;
m.d2ydx2  = @d2ydx2;
m.d2ydu2  = @d2ydu2;
m.d2ydk2  = @d2ydk2;
m.d2ydudx = @d2ydudx;
m.d2ydxdu = @d2ydxdu;
m.d2ydkdx = @d2ydkdx;
m.d2ydxdk = @d2ydxdk;
m.d2ydkdu = @d2ydkdu;
m.d2ydudk = @d2ydudk;

%% Process parameters
% Put rate parameter into k vector
m.k = vec([m.Parameters.Value]);

%% Process Reactions
m.rOrder = zeros(nr,1);
m.krInd  = zeros(nr,1);

nSEntries  = 0;
nSuEntries = 0; % Stoichiometry matrix for inputs used only for duplicate reaction checking
nD1Entries = 0;
nD2Entries = 0;
nD3Entries = 0;
nD4Entries = 0;
nD5Entries = 0;
nD6Entries = 0;
ndEntries  = 0;

SEntries  = zeros(0,2);
SValues   = zeros(0,1);
SuEntries  = zeros(0,2);
SuValues   = zeros(0,1);
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

for ir = 1:nr
    % Determine reaction type and species indexes
    n_reac = numel(m.Reactions(ir).Reactants);
    n_prod = numel(m.Reactions(ir).Products);
    reactants          = zeros(n_reac,1);
    reactantsExist     = false(n_reac,1);
    reactantsAreStates = false(n_reac,1);
    products           = zeros(n_prod,1);
    productsAreStates  = false(n_prod,1);
    for i_reac = 1:n_reac
        if isempty(m.Reactions(ir).Reactants{i_reac})
            reactantsExist(i_reac) = false;
        else
            reactantsExist(i_reac) = true;
            reactantsAreStates(i_reac) = true;
            temp = find(strcmp(m.Reactions(ir).Reactants(i_reac), x_full_names), 1);
            if ~isempty(temp)
                % Store index to state species
                reactants(i_reac) = temp;
                reactantsAreStates(i_reac) = true;
            else
                % Store index to input species
                reactants(i_reac) = find(strcmp(m.Reactions(ir).Reactants(i_reac), u_full_names), 1);
                reactantsAreStates(i_reac) = false;
            end
        end
    end
    
    m.rOrder(ir) = nnz(reactantsExist); % Store order
    
    for i_prod = 1:n_prod
        temp = find(strcmp(m.Reactions(ir).Products{i_prod}, x_full_names), 1);
        if ~isempty(temp)
            % It's a state
            products(i_prod) = temp;
            productsAreStates(i_prod) = true;
        else
            % It's an input
            temp = find(strcmp(m.Reactions(ir).Products{i_prod}, u_full_names), 1);
            products(i_prod) = temp;
            productsAreStates(i_prod) = false;
        end
    end
    
    % Extract parameter
    parameter = find(strcmp(m.Reactions(ir).Parameter{1}, k_names));
    assert(~isempty(parameter), 'KroneckerBio:FinalizeModel:MissingReactionParameter', 'Reaction %s (#%i) requires parameter %s, but no parameter by that name was found', m.Reactions(ir).Name, ir, m.Reactions(ir).Parameter{1})
    m.krInd(ir) = parameter; % Store index
    
    % The multiplier that is applied to the rate parameter
    modifier = m.Reactions(ir).Parameter{2};
    
    % Switch on reactant state
    if m.rOrder(ir) == 1 && reactantsAreStates(1)
        % D1/A1 reaction
        % Add S entries
        nAdd = 1 + nnz(productsAreStates);
        nuAdd = nnz(~productsAreStates);
        nSEntries = nSEntries + nAdd;
        nSuEntries = nSuEntries + nuAdd;
        
        % Add more room in vector if necessary
        SEntries = make_room(SEntries, nSEntries);
        SValues  = make_room(SValues, nSEntries);
        SuEntries = make_room(SuEntries, nSuEntries);
        SuValues = make_room(SuValues, nSuEntries);
        
        % Subtract reactant
        SEntries(nSEntries-nAdd+1,1) = reactants(1);
        SEntries(nSEntries-nAdd+1,2) = ir;
        SValues(nSEntries-nAdd+1)    = -1;
        
        % Add products
        states_added_so_far = cumsum(productsAreStates);
        inputs_added_so_far = cumsum(~productsAreStates);
        for i_prod = 1:n_prod
            if productsAreStates(i_prod)
                Srow = nSEntries-nAdd+1+states_added_so_far(i_prod);
                SEntries(Srow,1) = products(i_prod);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            else
                Srow = nSuEntries-nuAdd+inputs_added_so_far(i_prod);
                SuEntries(Srow,1) = products(i_prod);
                SuEntries(Srow,2) = ir;
                SuValues(Srow)    = 1;
            end
        end
        
        % Add D1 entry
        nD1Entries = nD1Entries + 1;
        
        % Add more room in vector if necessary
        dD1dkEntries = make_room(dD1dkEntries, nD1Entries);
        dD1dkValues  = make_room(dD1dkValues, nD1Entries);
        
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
        nuAdd = nnz(~productsAreStates);
        nSEntries = nSEntries + nAdd;
        nSuEntries = nSuEntries + nuAdd;
        
        % Add more room in vector if necessary
        SEntries = make_room(SEntries, nSEntries);
        SValues  = make_room(SValues, nSEntries);
        SuEntries = make_room(SuEntries, nSuEntries);
        SuValues = make_room(SuValues, nSuEntries);
        
        % Subtract reactant 1
        SEntries(nSEntries-nAdd+1,1) = reactants(1);
        SEntries(nSEntries-nAdd+1,2) = ir;
        SValues(nSEntries-nAdd+1)    = -1;
        
        % Subtract reactant 2
        SEntries(nSEntries-nAdd+2,1) = reactants(2);
        SEntries(nSEntries-nAdd+2,2) = ir;
        SValues(nSEntries-nAdd+2)    = -1;

        % Add products
        states_added_so_far = cumsum(productsAreStates);
        inputs_added_so_far = cumsum(~productsAreStates);
        for i_prod = 1:n_prod
            if productsAreStates(i_prod)
                Srow = nSEntries-nAdd+2+states_added_so_far(i_prod);
                SEntries(Srow,1) = products(i_prod);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            else
                Srow = nSuEntries-nuAdd+inputs_added_so_far(i_prod);
                SuEntries(Srow,1) = products(i_prod);
                SuEntries(Srow,2) = ir;
                SuValues(Srow)    = 1;
            end
        end
        
        % Add D2 entry
        nD2Entries = nD2Entries + 1;
        
        % Add more room in vector if necessary
        dD2dkEntries = make_room(dD2dkEntries, nD2Entries);
        dD2dkValues  = make_room(dD2dkValues, nD2Entries);
        
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
            nuAdd = 1 + nnz(~productsAreStates);
            nSEntries = nSEntries + nAdd;
            nSuEntries = nSuEntries + nuAdd;
            
            % Add more room in vector if necessary
            SEntries = make_room(SEntries, nSEntries);
            SValues  = make_room(SValues, nSEntries);
            SuEntries = make_room(SuEntries, nSuEntries);
            SuValues = make_room(SuValues, nSuEntries);
            
            % Subtract reactant 1
            SuEntries(nSuEntries-nuAdd+1,1) = reactants(1);
            SuEntries(nSuEntries-nuAdd+1,2) = ir;
            SuValues(nSuEntries-nuAdd+1)    = -1;

            % Subtract reactant 2
            SEntries(nSEntries-nAdd+1,1) = reactants(2);
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = -1;
            
            % Add products
            states_added_so_far = cumsum(productsAreStates);
            inputs_added_so_far = cumsum(~productsAreStates);
            for i_prod = 1:n_prod
                if productsAreStates(i_prod)
                    Srow = nSEntries-nAdd+1+states_added_so_far(i_prod);
                    SEntries(Srow,1) = products(i_prod);
                    SEntries(Srow,2) = ir;
                    SValues(Srow)    = 1;
                else
                    Srow = nSuEntries-nuAdd+1+inputs_added_so_far(i_prod);
                    SuEntries(Srow,1) = products(i_prod);
                    SuEntries(Srow,2) = ir;
                    SuValues(Srow)    = 1;
                end
            end
        
            % Add D3 entry
            nD3Entries = nD3Entries + 1;
            
            % Add more room in vector if necessary
            dD3dkEntries = make_room(dD3dkEntries, nD3Entries);
            dD3dkValues  = make_room(dD3dkValues, nD3Entries);
            
            dD3dkEntries(nD3Entries,1) = sub2ind([nr,nx,nu], ir, reactants(2), reactants(1));
            dD3dkEntries(nD3Entries,2) = parameter;
            dD3dkValues(nD3Entries)    = modifier;
            
        else%reactantsAreStates(1)
            % D4/A4 reaction
            % Add S entries
            nAdd = 1 + nnz(productsAreStates);
            nuAdd = 1 + nnz(~productsAreStates);
            nSEntries = nSEntries + nAdd;
            nSuEntries = nSuEntries + nuAdd;
            
            % Add more room in vector if necessary
            SEntries = make_room(SEntries, nSEntries);
            SValues  = make_room(SValues, nSEntries);
            SuEntries = make_room(SuEntries, nSuEntries);
            SuValues = make_room(SuValues, nSuEntries);
            
            % Subtract reactant 1
            SEntries(nSEntries-nAdd+1,1) = reactants(1);
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = -1;
            
            % Subtract reactant 2
            SuEntries(nSuEntries-nuAdd+1,1) = reactants(2);
            SuEntries(nSuEntries-nuAdd+1,2) = ir;
            SuValues(nSuEntries-nuAdd+1)    = -1;

            % Add products
            states_added_so_far = cumsum(productsAreStates);
            inputs_added_so_far = cumsum(~productsAreStates);
            for i_prod = 1:n_prod
                if productsAreStates(i_prod)
                    Srow = nSEntries-nAdd+1+states_added_so_far(i_prod);
                    SEntries(Srow,1) = products(i_prod);
                    SEntries(Srow,2) = ir;
                    SValues(Srow)    = 1;
                else
                    Srow = nSuEntries-nuAdd+1+inputs_added_so_far(i_prod);
                    SuEntries(Srow,1) = products(i_prod);
                    SuEntries(Srow,2) = ir;
                    SuValues(Srow)    = 1;
                end
            end
        
            % Add D4 entry
            nD4Entries = nD4Entries + 1;
        
            % Add more room in vector if necessary
            dD4dkEntries = make_room(dD4dkEntries, nD4Entries);
            dD4dkValues  = make_room(dD4dkValues, nD4Entries);
        
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
        nuAdd = 2 + nnz(~productsAreStates);
        nSEntries = nSEntries + nAdd;
        nSuEntries = nSuEntries + nuAdd;
        
        % Add more room in vector if necessary
        SEntries = make_room(SEntries, nSEntries);
        SValues  = make_room(SValues, nSEntries);
        SuEntries = make_room(SuEntries, nSuEntries);
        SuValues = make_room(SuValues, nSuEntries);
        
        % Subtract reactant 1
        SuEntries(nSuEntries-nuAdd+1,1) = reactants(1);
        SuEntries(nSuEntries-nuAdd+1,2) = ir;
        SuValues(nSuEntries-nuAdd+1)    = -1;

        % Subtract reactant 2
        SuEntries(nSuEntries-nuAdd+2,1) = reactants(2);
        SuEntries(nSuEntries-nuAdd+2,2) = ir;
        SuValues(nSuEntries-nuAdd+2)    = -1;

        % Add products
        states_added_so_far = cumsum(productsAreStates);
        inputs_added_so_far = cumsum(~productsAreStates);
        for i_prod = 1:n_prod
            if productsAreStates(i_prod)
                Srow = nSEntries-nAdd+states_added_so_far(i_prod);
                SEntries(Srow,1) = products(i_prod);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            else
                Srow = nSuEntries-nuAdd+2+inputs_added_so_far(i_prod);
                SuEntries(Srow,1) = products(i_prod);
                SuEntries(Srow,2) = ir;
                SuValues(Srow)    = 1;
            end
        end
        
        % Add D5 entry
        nD5Entries = nD5Entries + 1;
        
        % Add more room in vector if necessary
        dD5dkEntries = make_room(dD5dkEntries, nD5Entries);
        dD5dkValues  = make_room(dD5dkValues, nD5Entries);
        
        dD5dkEntries(nD5Entries,1) = sub2ind([nr,nu,nu], ir, reactants(2), reactants(1));
        dD5dkEntries(nD5Entries,2) = parameter;
        dD5dkValues(nD5Entries)    = modifier;
        
    elseif m.rOrder(ir) == 1 && ~reactantsAreStates(1)
        % D6/A6 reaction
        % Add S entries
        nAdd = nnz(productsAreStates);
        nuAdd = 1 + nnz(~productsAreStates);
        nSEntries = nSEntries + nAdd;
        nSuEntries = nSuEntries + nuAdd;
        
        % Add more room in vector if necessary
        SEntries = make_room(SEntries, nSEntries);
        SValues  = make_room(SValues, nSEntries);
        SuEntries = make_room(SuEntries, nSuEntries);
        SuValues = make_room(SuValues, nSuEntries);
        
        % Subtract reactant 1
        SuEntries(nSuEntries-nuAdd+1,1) = reactants(1);
        SuEntries(nSuEntries-nuAdd+1,2) = ir;
        SuValues(nSuEntries-nuAdd+1)    = -1;

        % Add products
        states_added_so_far = cumsum(productsAreStates);
        inputs_added_so_far = cumsum(~productsAreStates);
        for i_prod = 1:n_prod
            if productsAreStates(i_prod)
                Srow = nSEntries-nAdd+states_added_so_far(i_prod);
                SEntries(Srow,1) = products(i_prod);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            else
                Srow = nSuEntries-nuAdd+1+inputs_added_so_far(i_prod);
                SuEntries(Srow,1) = products(i_prod);
                SuEntries(Srow,2) = ir;
                SuValues(Srow)    = 1;
            end
        end
        
        % Add D6 entry
        nD6Entries = nD6Entries + 1;
        
        % Add more room in vector if necessary
        dD6dkEntries = make_room(dD6dkEntries, nD6Entries);
        dD6dkValues  = make_room(dD6dkValues, nD6Entries);
        
        dD6dkEntries(nD6Entries,1) = sub2ind([nr,nu], ir, reactants(1));
        dD6dkEntries(nD6Entries,2) = parameter;
        dD6dkValues(nD6Entries)    = modifier;
        
    elseif m.rOrder(ir) == 0
        % d/a reaction
        % Add S entries
        nAdd = nnz(productsAreStates);
        nuAdd = nnz(~productsAreStates);
        nSEntries = nSEntries + nAdd;
        nSuEntries = nSuEntries + nuAdd;
        
        % Add more room in vector if necessary
        SEntries = make_room(SEntries, nSEntries);
        SValues  = make_room(SValues, nSEntries);
        SuEntries = make_room(SuEntries, nSuEntries);
        SuValues = make_room(SuValues, nSuEntries);
        
        % Add products
        states_added_so_far = cumsum(productsAreStates);
        inputs_added_so_far = cumsum(~productsAreStates);
        for i_prod = 1:n_prod
            if productsAreStates(i_prod)
                Srow = nSEntries-nAdd+states_added_so_far(i_prod);
                SEntries(Srow,1) = products(i_prod);
                SEntries(Srow,2) = ir;
                SValues(Srow)    = 1;
            else
                Srow = nSuEntries-nuAdd+inputs_added_so_far(i_prod);
                SuEntries(Srow,1) = products(i_prod);
                SuEntries(Srow,2) = ir;
                SuValues(Srow)    = 1;
            end
        end
        
        % Add d entry
        ndEntries = ndEntries + 1;
        
        % Add more room in vector if necessary
        dddkEntries = make_room(dddkEntries, ndEntries);
        dddkValues  = make_room(dddkValues, ndEntries);
        
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
Su = sparse(SuEntries(1:nSuEntries,1), SuEntries(1:nSuEntries,2), SuValues(1:nSuEntries), nu, nr);

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
%m.A2 = reshape(m.dA2dk * kRand, nx,nx*nx);
m.A2 = reshape(mtimestall(m.dA2dk, kRand), nx,nx*nx);
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
[~, D2UsedColumns] = find(m.D2);
[~, D3UsedColumns] = find(m.D3);
[~, D4UsedColumns] = find(m.D4);
[~, D5UsedColumns] = find(m.D5);

D2UsedColumns = vec(unique(D2UsedColumns)); % vec compensates for Matlab bug when nr = 1
D3UsedColumns = vec(unique(D3UsedColumns));
D4UsedColumns = vec(unique(D4UsedColumns));
D5UsedColumns = vec(unique(D5UsedColumns));

% *1 and *2 refer to the order in (1 kron 2)
[D2UsedSpecies2, D2UsedSpecies1] = ind2sub([nx,nx], D2UsedColumns);
[D3UsedSpecies2, D3UsedSpecies1] = ind2sub([nx,nu], D3UsedColumns);
[D4UsedSpecies2, D4UsedSpecies1] = ind2sub([nu,nx], D4UsedColumns);
[D5UsedSpecies2, D5UsedSpecies1] = ind2sub([nu,nu], D5UsedColumns);

%% Warn on duplicate reactions
% Encode reaction names as numbers
[~, ~, i_reaction_name] = unique(r_names);

% Make one index zero so the sparse drops it (usually the empty string)
i_reaction_name = sparse(i_reaction_name - 1);

% Assemble reactions into one big matrix
reaction_matrix = [sparse(m.krInd), m.S.', Su.', i_reaction_name];

% All rows that already appeared elsewhere
[~, i_unique_reactions] = unique(reaction_matrix, 'rows');
i_duplicates = setdiff(vec(1:nr), i_unique_reactions);

% Remove duplicate duplicates
[~, i_unique_i_duplicates] = unique(reaction_matrix(i_duplicates, :), 'rows');
i_unique_duplicates = i_duplicates(i_unique_i_duplicates);

for i_dup = row(i_unique_duplicates)
    % Find all rows corresponding to this duplicate and warn
    dup_inds = find(ismember(reaction_matrix, reaction_matrix(i_dup,:), 'rows'));
    
    warning_string = 'The following reactions are identical: ';
    for dup_ind = row(dup_inds)
        if isempty(r_names{dup_ind})
            warning_string = [warning_string, ' #', num2str(dup_ind), ','];
        else
            warning_string = [warning_string, ' #', num2str(dup_ind), ' "', r_names{dup_ind}, '",'];
        end
    end
    warning_string = warning_string(1:end-1); % Drop final comma
    
    warning('KroneckerBio:FinalizeModel:RepeatReactions', warning_string)
end

%% Final build of model
m = final(m, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Function handles not dependent on parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function val = v(t, x, u)
        val = B1 * x + B2 * u + b;
    end

    function val = dvdx(t, x, u)
        val = B1;
    end

    function val = dvdu(t, x, u)
        val = B2;
    end

    function val = d2vdx2(t, x, u)
        val = sparse(nv*nx, nx);
    end

    function val = d2vdu2(t, x, u)
        val = sparse(nv*nu, nu);
    end

    function val = d2vdudx(t, x, u)
        val = sparse(nv*nx, nu);
    end

    function val = d2vdxdu(t, x, u)
        val = sparse(nv*nu, nx);
    end

    function val = y(t, x, u)
        val = C1 * x + C2 * u + repmat(c, [1,numel(t)]);
    end

    function val = dydx(t, x, u)
        val = C1;
    end

    function val = dydu(t, x, u)
        val = C2;
    end

    function val = dydk(t, x, u)
        val = sparse(ny, nk);
    end

    function val = d2ydx2(t, x, u)
        val = sparse(ny*nx, nx);
    end

    function val = d2ydu2(t, x, u)
        val = sparse(ny*nu, nu);
    end

    function val = d2ydk2(t, x, u)
        val = sparse(ny*nk, nk);
    end

    function val = d2ydudx(t, x, u)
        val = sparse(ny*nx, nu);
    end

    function val = d2ydxdu(t, x, u)
        val = sparse(ny*nu, nx);
    end

    function val = d2ydkdx(t, x, u)
        val = sparse(ny*nx, nk);
    end

    function val = d2ydxdk(t, x, u)
        val = sparse(ny*nk, nx);
    end

    function val = d2ydkdu(t, x, u)
        val = sparse(ny*nu, nk);
    end

    function val = d2ydudk(t, x, u)
        val = sparse(ny*nk, nu);
    end

    function val = x0(s)
        val = dx0ds_val*s + x0c;
    end

    function val = dx0ds(s)
        val = dx0ds_val;
    end

    function val = dx0dk(s)
        val = sparse(nx,nk);
    end

    function val = d2x0ds2(s)
        val = sparse(nx*ns, ns);
    end

    function val = d2x0dk2(s)
        val = sparse(nx*nk, nk);
    end

    function val = d2x0dkds(s)
        val = sparse(nx*ns, nk);
    end

    function val = d2x0dsdk(s)
        val = sparse(nx*nk, ns);
    end
end

function m = final(m, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2)

% Constants
nx = m.nx;
nu = m.nu;
nr = m.nr;

% Build kronecker matrices
sparsek = sparse(m.k);
m.A1 = reshape(m.dA1dk * sparsek, nx,nx);
%m.A2 = reshape(m.dA2dk * sparsek, nx,nx*nx);
m.A2 = reshape(mtimestall(m.dA2dk, sparsek), nx,nx*nx);
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

m.dfdx = dfdxHidden(m.A1, m.A2, m.A3, m.A4, m.A5, m.v, m.dvdx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.dfdu = dfduHidden(m.A2, m.A3, m.A4, m.A5, m.A6, m.v, m.dvdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.dfdk = dfdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.dadk, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.d2fdx2  = d2fdx2Hidden(m.A2, m.A3, m.A4, m.A5, m.v, m.dvdx, m.d2vdx2, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdu2  = d2fdu2Hidden(m.A2, m.A3, m.A4, m.A5, m.v, m.dvdu, m.d2vdu2, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdk2  = d2fdk2Hidden(m.nx, m.nk);
m.d2fdudx = d2fdudxHidden(m.A2, m.A3, m.A4, m.A5, m.v, m.dvdx, m.dvdu, m.d2vdudx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdxdu = d2fdxduHidden(m.A2, m.A3, m.A4, m.A5, m.v, m.dvdx, m.dvdu, m.d2vdxdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdkdx = d2fdkdxHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.v, m.dvdx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdxdk = d2fdxdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.v, m.dvdx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdkdu = d2fdkduHidden(m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.v, m.dvdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdudk = d2fdudkHidden(m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.v, m.dvdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.r = rHidden(m.D1, m.D2, m.D3, m.D4, m.D5, m.D6, m.d, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.drdx = drdxHidden(m.D1, m.D2, m.D3, m.D4, m.D5, m.v, m.dvdx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.drdu = drduHidden(m.D2, m.D3, m.D4, m.D5, m.D6, m.v, m.dvdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.drdk = drdkHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.dD5dk_rk_uu, m.dD6dk_rk_u, m.dddk, m.v, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.d2rdx2  = d2rdx2Hidden(m.D2, m.D3, m.D4, m.D5, m.v, m.dvdx, m.d2vdx2, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2rdu2  = d2rdu2Hidden(m.D2, m.D3, m.D4, m.D5, m.v, m.dvdu, m.d2vdu2, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2rdk2  = d2rdk2Hidden(m.nx, m.nk);
m.d2rdudx = d2rdudxHidden(m.D2, m.D3, m.D4, m.D5, m.v, m.dvdx, m.dvdu, m.d2vdudx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2rdxdu = d2rdxduHidden(m.D2, m.D3, m.D4, m.D5, m.v, m.dvdx, m.dvdu, m.d2vdxdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2rdkdx = d2rdkdxHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.dD5dk_rk_uu, m.v, m.dvdx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd, m.nr);
m.d2rdxdk = d2rdxdkHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.dD5dk_rk_uu, m.v, m.dvdx, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2rdkdu = d2rdkduHidden(m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.dD5dk_rk_uu, m.dD6dk_rk_u, m.v, m.dvdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd, m.nr);
m.d2rdudk = d2rdudkHidden(m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.dD5dk_rk_uu, m.dD6dk_rk_u, m.v, m.dvdu, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.Ready = true;
m.Update = @Update;

    function mout = Update(k)
        % Copy existing model
        mout = m;
        
        % Apply changes
        mout.k = k;
        
        % Distribute values
        if m.nk >= 1
            k = num2cell(k);
            [mout.Parameters.Value] = k{:};
        end
        
        % Rebuild model
        mout = final(mout, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2);
    end
end

function handle = fHidden(A1, A2, A3, A4, A5, A6, a, v_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @f;

    function val = f(t, x, u)
%   f = A1*x + A2*(x kron x/vx) + A3*(u kron x/vx) + A4*(x kron u/vu) + A5*(u kron u/vu) + A6*u + a
        % Compartment column
        v = v_fun(t,x,u);
        xvxinv = x ./ v(vxInd);
        uvuinv = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvxinv(D2UsedSpecies2), nx*nx,1);
        ukronx = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvxinv(D3UsedSpecies2), nx*nu,1);
        xkronu = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvuinv(D4UsedSpecies2), nu*nx,1);
        ukronu = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvuinv(D5UsedSpecies2), nu*nu,1);
        
        val = A1 * x + A2 * xkronx + A3 * ukronx + A4 * xkronu + A5 * ukronu + A6 * u + a; % f_
    end
end

function handle = dfdxHidden(A1, A2, A3, A4, A5, v_fun, dvdx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdx;

    function val = dfdx(t, x, u)
%   dfdx = A1 + A2*(Ix kron x/vx) + A2*(x kron diag(1/vx)) + A3*(u kron diag(1/vx))
        % Compartment column
        v = v_fun(t,x,u);
        dvdx = dvdx_fun(t,x,u);
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1));
        vxinv = 1 ./ v(vxInd);
        xvxinv = x .* vxinv;
        uvuinv = u ./ v(vuInd);
        dvxinvdx = dvinvdx(vxInd,:);
        dvuinvdx = dvinvdx(vuInd,:);
        
        % A2
        Ixkronxvxinv   = sparse(D2UsedColumns, D2UsedSpecies1, xvxinv(D2UsedSpecies2), nx*nx,nx);
        xkronIxvxinv   = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        xkronxdvxinvdx = sparse(nx*nx,nx);
        xkronxdvxinvdx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:));
        
        % A3
        ukronIxvxinv   = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        ukronxdvxinvdx = sparse(nu*nx,nx);
        ukronxdvxinvdx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:));
        
        % A4
        Ixkronuvuinv   = sparse(D4UsedColumns, D4UsedSpecies1, uvuinv(D4UsedSpecies2), nx*nu,nx);
        xkronudvuinvdx = sparse(nx*nu,nx);
        xkronudvuinvdx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:));
        
        % A5
        ukronudvuinvdx = sparse(nu*nu,nx);
        ukronudvuinvdx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:));
        
        val = A1 + A2 * (Ixkronxvxinv + xkronIxvxinv + xkronxdvxinvdx) + A3 * (ukronIxvxinv + ukronxdvxinvdx) + A4 * (Ixkronuvuinv + xkronudvuinvdx) + A5 * ukronudvuinvdx; % f_x
    end
end

function handle = dfduHidden(A2, A3, A4, A5, A6, v_fun, dvdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdu;

    function val = dfdu(t, x, u)
%   dfdu = A3*(Iu kron x/vx) + A4*(x kron diag(1/vu)) + A5*(Iu kron u/vu) + A5*(u kron diag(1/vu)) + A6
        % Compartment column
        v = v_fun(t,x,u);
        dvdu = dvdu_fun(t,x,u);
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1));
        xvxinv = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvuinv = u .* vuinv;
        dvxinvdu = dvinvdu(vxInd,:);
        dvuinvdu = dvinvdu(vuInd,:);
        
        % A2
        xkronxdvxinvdu = sparse(nx*nx,nu);
        xkronxdvxinvdu(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:));

        % A3
        Iukronxvxinv   = sparse(D3UsedColumns, D3UsedSpecies1, xvxinv(D3UsedSpecies2), nx*nu,nu);
        ukronxdvxinvdu = sparse(nu*nx,nu);
        ukronxdvxinvdu(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:));
        
        % A4
        xkronIuvuinv   = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        xkronudvuinvdu = sparse(nx*nu,nu);
        xkronudvuinvdu(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:));

        % A5
        Iukronuvuinv   = sparse(D5UsedColumns, D5UsedSpecies1, uvuinv(D5UsedSpecies2), nu*nu,nu);
        ukronIuvuinv   = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        ukronudvuinvdu = sparse(nu*nu,nu);
        ukronudvuinvdu(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:));
        
        val = A2 * xkronxdvxinvdu + A3 * (Iukronxvxinv + ukronxdvxinvdu) + A4 * (xkronIuvuinv + xkronudvuinvdu) + A5 * (Iukronuvuinv + ukronIuvuinv + ukronudvuinvdu) + A6; % f_u
    end
end

function handle = dfdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, dadk, v_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @dfdk;

    function val = dfdk(t, x, u)
%   dfdk = dA1dk*x + dA2dk*(x kron x/vx) + dA3dk*(u kron x/vx) + dA4dk*(x kron u/vu) + dA5dk*(u kron u/vu) + dA6dk*u + dadk
        % Compartment column
        v = v_fun(t,x,u);
        xvxinv = x ./ v(vxInd);
        uvuinv = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronxinv = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvxinv(D2UsedSpecies2), nx*nx,1);
        ukronxinv = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvxinv(D3UsedSpecies2), nx*nu,1);
        xkronuinv = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvuinv(D4UsedSpecies2), nu*nx,1);
        ukronuinv = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvuinv(D5UsedSpecies2), nu*nu,1);
        
        val = dA1dk_fk_x * sparse(x) + dA2dk_fk_xx * xkronxinv + dA3dk_fk_ux * ukronxinv + dA4dk_fk_xu * xkronuinv + dA5dk_fk_uu * ukronuinv + dA6dk_fk_u * sparse(u); % fk_
        val = reshape(val, nx,nk) + dadk; % f_k
    end
end

% Notes on second-order derivatives:
%   When a derivative has the form Dd * (Ia kron Ib / vbinv)
%   with dimensions nb*na by nb*na:
%   IakronIbvbinv = sparse(DdUsedColumns, DdUsedColumns, vbinv(DdUsedSpecies2), nb*na,nb*na);
%
%   When a derivative has the form Dd * (Ia kron Ib / vbinv)
%   with dimensions nb*na by na*nb:
%   DdUsedColumnsReverse = linearindexpermute(DdUsedColumns, [nb,na], [2,1]);
%   IakronIbvbinv = sparse(DdUsedColumns, DdUsedColumnsReverse, vbinv(DdUsedSpecies2), nb*na,na*nb);
%
%   When a derivative has the form Dd * (Ia kron b / dvbinvdc)
%   with dimensions nb*na by nc*na:
%   rows = repmat(DdUsedColumns, [nc,1]);
%   cols = vec(bsxfun(@plus, 1:nc, (DdUsedSpecies1-1)*nc));
%   Iakronbdvbinvdc = sparse(rows, cols, bsxfun(@times, b(DdUsedSpecies2), dvbinvdc(DdUsedSpecies2,:)), nb*na,nc*nb);
%
%   When a derivative has the form Dd * (Ia kron b / dvbinvdc)
%   with dimensions nb*na by na*nc:
%   rows = repmat(DdUsedColumns, [nc,1]);
%   cols = lineaerindexpermute(vec(bsxfun(@plus, 1:nc, (DdUsedSpecies1-1)*nc)), [nc,na], [2,1]);
%   Iakronbdvbinvdc = sparse(rows, cols, bsxfun(@times, b(DdUsedSpecies2), dvbinvdc(DdUsedSpecies2,:)), nb*na,na*nc);
%
%   When a derivative has the form Dd * (a kron Ib / dvbinvdc)
%   with dimensions nb*na by nb*nc:
%   rows = repmat(DdUsedColumns, [nc,1]);
%   cols = vec(bsxfun(@plus, 1:nb:nb*nc, DdUsedSpecies2-1));
%   akronIbdvbinvdc = sparse(rows, cols, bsxfun(@times, a(DdUsedSpecies1), dvbinvdc(DdUsedSpecies2,:)), nb*na,nb*nc);
%
%   When a derivative has the form Dd * (a kron Ib / dvbinvdc)
%   with dimensions nb*na by nc*nb:
%   rows = repmat(DdUsedColumns, [nc,1]);
%   cols = linearindexpermute(vec(bsxfun(@plus, 1:nb:nb*nc, DdUsedSpecies2-1)), [nb,nc], [2,1]);
%   akronIbdvbinvdc = sparse(rows, cols, bsxfun(@times, a(DdUsedSpecies1), dvbinvdc(DdUsedSpecies2,:)), nb*na,nc*nb);

function handle = d2fdx2Hidden(A2, A3, A4, A5, v_fun, dvdx_fun, d2vdx2_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% A2
D2UsedColumnsReverse = linearindexpermute(D2UsedColumns, [nx,nx], [2,1]);

Ixkronxdvxinvdx_rows = repmat(D2UsedColumns,[nx,1]);
Ixkronxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx, (D2UsedSpecies1-1)*nx));
 
xkronIxdvxinvdx_rows = Ixkronxdvxinvdx_rows;
xkronIxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx:nx*nx, D2UsedSpecies2-1));

% A3
ukronIxdvxinvdx_rows = repmat(D3UsedColumns,[nx,1]);
ukronIxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx:nx*nx, D3UsedSpecies2-1));

% A4
Ixkronudvuinvdx_rows = repmat(D4UsedColumns,[nx,1]);
Ixkronudvuinvdx_cols = vec(bsxfun(@plus, 1:nx, (D4UsedSpecies1-1)*nx));

% Return handle
handle = @d2fdx2;

    function val = d2fdx2(t, x, u)
%   d2fdx2 = 2*A2*(Ix kron diag(1/vx))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdx = dvdx_fun(t,x,u); % v_x
        d2vdx2 = d2vdx2_fun(t,x,u); % vx_x
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1)); % v_x
        d2vinvdx2 = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdx .* -2)), nv,1,nx), full(dvdx) .* -1), nv,nx*nx) + bsxfun(@times, (v .^ -2), reshape(d2vdx2, nv,nx*nx) .* -1)); % v_xx
        vxinv = 1 ./ v(vxInd); % x_
        dvxinvdx = dvinvdx(vxInd,:); % x_x
        d2vxinvdx2 = d2vinvdx2(vxInd,:); % x_xx
        dvuinvdx = dvinvdx(vuInd,:); % u_x
        d2vuinvdx2 = d2vinvdx2(vuInd,:); % u_xx
        
        % A2
        %IxkronIxvxinv = sparse([D2UsedColumns; D2UsedColumns], [D2UsedColumns; D2UsedColumnsReverse], [vxinv(D2UsedSpecies2); vxinv(D2UsedSpecies2)], nx*nx,nx*nx);
        IxkronIxvxinv_1 = sparse(D2UsedColumns, D2UsedColumnsReverse, vxinv(D2UsedSpecies2), nx*nx,nx*nx);
        IxkronIxvxinv_2 = sparse(D2UsedColumns, D2UsedColumns, vxinv(D2UsedSpecies2), nx*nx,nx*nx);
        
        %Ixkronxdvxinvdx = kron(eye(nx), diag(x) * dvxinvdx);
        Ixkronxdvxinvdx_1 = sparse(Ixkronxdvxinvdx_rows, Ixkronxdvxinvdx_cols, bsxfun(@times, x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:)), nx*nx,nx*nx);
        Ixkronxdvxinvdx_2 = Ixkronxdvxinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        %xkronIxdvxinvdx = reshape(kron(x, reshape(bsxfun(@times, eye(nx), reshape(full(dvxinvdx),[nx,1,nx])), [nx,nx*nx])), [nx*nx,nx*nx]);
        xkronIxdvxinvdx_1 = sparse(xkronIxdvxinvdx_rows, xkronIxdvxinvdx_cols, bsxfun(@times, x(D2UsedSpecies1), dvxinvdx(D2UsedSpecies2,:)), nx*nx,nx*nx);
        xkronIxdvxinvdx_2 = xkronIxdvxinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        %xkronxd2vxinvdx2 = kron(x, diag(x) * d2vxinvdx2);
        xkronxd2vxinvdx2 = sparse(nx*nx,nx*nx);
        xkronxd2vxinvdx2(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdx2(D2UsedSpecies2,:));

        % A3
        ukronIxdvxinvdx_1 = sparse(ukronIxdvxinvdx_rows, ukronIxdvxinvdx_cols, bsxfun(@times, u(D3UsedSpecies1), dvxinvdx(D3UsedSpecies2,:)), nx*nu,nx*nx);
        ukronIxdvxinvdx_2 = ukronIxdvxinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        ukronxd2vxinvdx2 = sparse(nx*nu,nx*nx);
        ukronxd2vxinvdx2(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdx2(D3UsedSpecies2,:));

        % A4
        Ixkronudvuinvdx_1 = sparse(Ixkronudvuinvdx_rows, Ixkronudvuinvdx_cols, bsxfun(@times, u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:)), nu*nx,nx*nx);
        Ixkronudvuinvdx_2 = Ixkronudvuinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        xkronud2vuinvdx2 = sparse(nu*nx,nx*nx);
        xkronud2vuinvdx2(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdx2(D4UsedSpecies2,:));
        
        % A5
        ukronud2vuinvdx2 = sparse(nu*nu,nx*nx);
        ukronud2vuinvdx2(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdx2(D5UsedSpecies2,:));

        val = A2 * (IxkronIxvxinv_1 + Ixkronxdvxinvdx_1 + IxkronIxvxinv_2 + xkronIxdvxinvdx_1 + Ixkronxdvxinvdx_2 + xkronIxdvxinvdx_2 + xkronxd2vxinvdx2)...
              + A3 * (ukronIxdvxinvdx_1 + ukronIxdvxinvdx_2 + ukronxd2vxinvdx2)...
              + A4 * (Ixkronudvuinvdx_1 + Ixkronudvuinvdx_2 + xkronud2vuinvdx2)...
              + A5 * ukronud2vuinvdx2; % f_xx
        val = reshape(val, nx*nx, nx); % fx_x
    end
end

function handle = d2fdu2Hidden(A2, A3, A4, A5, v_fun, dvdu_fun, d2vdu2_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% A3
Iukronxdvxinvdu_rows = repmat(D3UsedColumns,[nu,1]);
Iukronxdvxinvdu_cols = vec(bsxfun(@plus, 1:nu, (D3UsedSpecies1-1)*nu));

% A4
xkronIudvuinvdu_rows = repmat(D4UsedColumns,[nu,1]);
xkronIudvuinvdu_cols = vec(bsxfun(@plus, 1:nu:nu*nu, D4UsedSpecies2-1));

% A5
[u1ind,u2ind] = ind2sub([nu,nu], D5UsedColumns);
D5UsedColumnsReverse = sub2ind([nu, nu], u2ind, u1ind);

Iukronudvuinvdu_rows = repmat(D5UsedColumns,[nu,1]);
Iukronudvuinvdu_cols = vec(bsxfun(@plus, 1:nu, (D5UsedSpecies1-1)*nu));
 
ukronIudvuinvdu_rows = Iukronudvuinvdu_rows;
ukronIudvuinvdu_cols = vec(bsxfun(@plus, 1:nu:nu*nu, D5UsedSpecies2-1));

% Return handle
handle = @d2fdu2;

    function val = d2fdu2(t, x, u)
%   d2fdu2 = 2*A5*(Iu kron diag(1/vu))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdu = dvdu_fun(t,x,u); % v_u
        d2vdu2 = d2vdu2_fun(t,x,u); % vu_u
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1)); % v_u
        d2vinvdu2 = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdu .* -2)), nv,1,nu), full(dvdu) .* -1), nv,nu*nu) + bsxfun(@times, (v .^ -2), reshape(d2vdu2, nv,nu*nu) .* -1)); % v_uu
        dvxinvdu = dvinvdu(vxInd,:); % x_u
        d2vxinvdu2 = d2vinvdu2(vxInd,:); % x_uu
        vuinv = 1 ./ v(vuInd); % u_
        dvuinvdu = dvinvdu(vuInd,:); % u_u
        d2vuinvdu2 = d2vinvdu2(vuInd,:); % u_uu
        
        % A2
        xkronxd2vxinvdu2 = sparse(nx*nx,nu*nu);
        xkronxd2vxinvdu2(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdu2(D2UsedSpecies2,:));

        % A3
        Iukronxdvxinvdu_1 = sparse(Iukronxdvxinvdu_rows, Iukronxdvxinvdu_cols, bsxfun(@times, x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:)), nx*nu,nu*nu);
        Iukronxdvxinvdu_2 = Iukronxdvxinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        ukronxd2vxinvdu2 = sparse(nx*nu,nu*nu);
        ukronxd2vxinvdu2(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdu2(D3UsedSpecies2,:));
        
        % A4
        xkronIudvuinvdu_1 = sparse(xkronIudvuinvdu_rows, xkronIudvuinvdu_cols, bsxfun(@times, x(D4UsedSpecies1), dvuinvdu(D4UsedSpecies2,:)), nu*nx,nu*nu);
        xkronIudvuinvdu_2 = xkronIudvuinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        xkronud2vuinvdu2 = sparse(nu*nx,nu*nu);
        xkronud2vuinvdu2(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdu2(D4UsedSpecies2,:));

        % A5
        IukronIuvuinv_1 = sparse(D5UsedColumns, D5UsedColumnsReverse, vuinv(D5UsedSpecies2), nu*nu,nu*nu);
        IukronIuvuinv_2 = sparse(D5UsedColumns, D5UsedColumns, vuinv(D5UsedSpecies2), nu*nu,nu*nu);
        
        %Iukronudvuinvdu = kron(eye(nu), diag(u) * dvuinvdu);
        Iukronudvuinvdu_1 = sparse(Iukronudvuinvdu_rows, Iukronudvuinvdu_cols, bsxfun(@times, u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:)), nu*nu,nu*nu);
        Iukronudvuinvdu_2 = Iukronudvuinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        %ukronIudvuinvdu = reshape(kron(u, reshape(bsxfun(@times, eye(nu), reshape(full(dvuinvdu),[nu,1,nu])), [nu,nu*nu])), [nu*nu,nu*nu]);
        ukronIudvuinvdu_1 = sparse(ukronIudvuinvdu_rows, ukronIudvuinvdu_cols, bsxfun(@times, u(D5UsedSpecies1), dvuinvdu(D5UsedSpecies2,:)), nu*nu,nu*nu);
        ukronIudvuinvdu_2 = ukronIudvuinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        %ukronud2vuinvdu2 = kron(u, diag(u) * d2vuinvdu2);
        ukronud2vuinvdu2 = sparse(nu*nu,nu*nu);
        ukronud2vuinvdu2(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdu2(D5UsedSpecies2,:));

        val = A2 * xkronxd2vxinvdu2...
              + A3 * (Iukronxdvxinvdu_1 + Iukronxdvxinvdu_2 + ukronxd2vxinvdu2)...
              + A4 * (xkronIudvuinvdu_1 + xkronIudvuinvdu_2 + xkronud2vuinvdu2)...
              + A5 * (IukronIuvuinv_1 + Iukronudvuinvdu_1 + IukronIuvuinv_2 + ukronIudvuinvdu_1 + Iukronudvuinvdu_2 + ukronIudvuinvdu_2 + ukronud2vuinvdu2); % f_uu
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

function handle = d2fdudxHidden(A2, A3, A4, A5, v_fun, dvdx_fun, dvdu_fun, d2vdudx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% A2
Ixkronxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
Ixkronxdvxinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu, (D2UsedSpecies1-1)*nu)), [nu,nx], [2,1]);

xkronIxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
xkronIxdvxinvdu_cols = vec(bsxfun(@plus, 1:nx:nx*nu, D2UsedSpecies2-1));

% A3
Iukronxdvxinvdx_rows = repmat(D3UsedColumns,[nx,1]);
Iukronxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx, (D3UsedSpecies1-1)*nx));

ukronIxdvxinvdu_rows = repmat(D3UsedColumns,[nu,1]);
ukronIxdvxinvdu_cols = vec(bsxfun(@plus, 1:nx:nx*nu, D3UsedSpecies2-1));

% A4
D4UsedColumnsReverse = linearindexpermute(D4UsedColumns, [nu, nx], [2,1]);

Ixkronudvuinvdu_rows = repmat(D4UsedColumns,[nu,1]);
Ixkronudvuinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu, (D4UsedSpecies1-1)*nu)), [nu,nx], [2,1]);

xkronIudvuinvdx_rows = repmat(D4UsedColumns,[nx,1]);
xkronIudvuinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu:nu*nx, D4UsedSpecies2-1)), [nu,nx], [2,1]);

% A5
Iukronudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
Iukronudvuinvdx_cols = vec(bsxfun(@plus, 1:nx, (D5UsedSpecies1-1)*nx));

ukronIudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
ukronIudvuinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu:nu*nx, D5UsedSpecies2-1)), [nu,nx], [2,1]);

% Return handle
handle = @d2fdudx;

    function val = d2fdudx(t, x, u)
%   d2fdudx = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdx = dvdx_fun(t,x,u); % v_x
        dvdu = dvdu_fun(t,x,u); % v_u
        d2vdudx = d2vdudx_fun(t,x,u); % vx_x

        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1)); % v_x
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1)); % v_u
        d2vinvdudx = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdu .* -2)), nv,1,nu), full(dvdx) .* -1), nv,nx*nu) + bsxfun(@times, (v .^ -2), reshape(d2vdudx, nv,nx*nu) .* -1)); % v_xu

        vxinv = 1 ./ v(vxInd); % x_
        dvxinvdx = dvinvdx(vxInd,:); % x_x
        dvxinvdu = dvinvdu(vxInd,:); % x_u
        d2vxinvdudx = d2vinvdudx(vxInd,:); % x_xu
        
        vuinv = 1 ./ v(vuInd); % u_
        dvuinvdx = dvinvdx(vuInd,:); % u_x
        dvuinvdu = dvinvdu(vuInd,:); % u_u
        d2vuinvdudx = d2vinvdudx(vuInd,:); % u_uu
        
        % A2 derivatives
        Ixkronxdvxinvdu = sparse(Ixkronxdvxinvdu_rows, Ixkronxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nx*nu);
        xkronIxdvxinvdu = sparse(xkronIxdvxinvdu_rows, xkronIxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies1), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nx*nu);
        
        xkronxd2vxinvdudx = sparse(nx*nx,nx*nu);
        xkronxd2vxinvdudx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdudx(D2UsedSpecies2,:));
        
        % A3 derivatives
        IukronIxvxinv   = sparse(D3UsedColumns, D3UsedColumns, vxinv(D3UsedSpecies2), nx*nu,nx*nu);
        ukronIxdvxinvdu = sparse(ukronIxdvxinvdu_rows, ukronIxdvxinvdu_cols, bsxfun(@times, u(D3UsedSpecies1), dvxinvdu(D3UsedSpecies2,:)), nx*nu,nx*nu);
        Iukronxdvxinvdu = sparse(Iukronxdvxinvdx_rows, Iukronxdvxinvdx_cols, bsxfun(@times, x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:)), nx*nu,nx*nu);

        ukronxd2vxinvdudx = sparse(nx*nu,nx*nu);
        ukronxd2vxinvdudx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdudx(D3UsedSpecies2,:));

        % A4 derivatives
        IxkronIuvuinv_1   = sparse(D4UsedColumns, D4UsedColumnsReverse, vuinv(D4UsedSpecies2), nu*nx,nx*nu);
        Ixkronudvuinvdu_1 = sparse(Ixkronudvuinvdu_rows, Ixkronudvuinvdu_cols, bsxfun(@times, u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:)), nu*nx,nx*nu);
        xkronIudvuinvdx_1 = sparse(xkronIudvuinvdx_rows, xkronIudvuinvdx_cols, bsxfun(@times, x(D4UsedSpecies1), dvuinvdx(D4UsedSpecies2,:)), nu*nx,nx*nu);

        xkronud2vuinvdudx = sparse(nu*nx,nx*nu);
        xkronud2vuinvdudx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdudx(D4UsedSpecies2,:));
        
        % A5 derivatives
        Iukronudvuinvdx = sparse(Iukronudvuinvdx_rows, Iukronudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nx*nu);
        ukronIudvuinvdx = sparse(ukronIudvuinvdx_rows, ukronIudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies1), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nx*nu);
        
        ukronud2vuinvdudx = sparse(nu*nu,nx*nu);
        ukronud2vuinvdudx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdudx(D5UsedSpecies2,:));

        val = A2 * (Ixkronxdvxinvdu + xkronIxdvxinvdu + xkronxd2vxinvdudx) ...
              + A3 * (IukronIxvxinv + ukronIxdvxinvdu + Iukronxdvxinvdu + ukronxd2vxinvdudx) ...
              + A4 * (IxkronIuvuinv_1 + Ixkronudvuinvdu_1 + xkronIudvuinvdx_1 + xkronud2vuinvdudx) ...
              + A5 * (Iukronudvuinvdx + ukronIudvuinvdx + ukronud2vuinvdudx);
        val = reshape(val, nx*nx,nu); % fx_u
    end
end

function handle = d2fdxduHidden(A2, A3, A4, A5, v_fun, dvdx_fun, dvdu_fun, d2vdxdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% A2
Ixkronxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
Ixkronxdvxinvdu_cols = vec(bsxfun(@plus, 1:nu, (D2UsedSpecies1-1)*nu));

xkronIxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
xkronIxdvxinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx:nx*nu, D2UsedSpecies2-1)), [nx,nu], [2,1]);

% A3
D3UsedColumnsReverse = linearindexpermute(D3UsedColumns, [nx,nu], [2,1]);

Iukronxdvxinvdx_rows = repmat(D3UsedColumns,[nx,1]);
Iukronxdvxinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx, (D3UsedSpecies1-1)*nx)), [nx,nu], [2,1]);

ukronIxdvxinvdu_rows = repmat(D3UsedColumns,[nu,1]);
ukronIxdvxinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx:nx*nu, D3UsedSpecies2-1)), [nx,nu], [2,1]);

% A4
Ixkronudvuinvdu_rows = repmat(D4UsedColumns,[nu,1]);
Ixkronudvuinvdu_cols = vec(bsxfun(@plus, 1:nu, (D4UsedSpecies1-1)*nu));

xkronIudvuinvdx_rows = repmat(D4UsedColumns,[nx,1]);
xkronIudvuinvdx_cols = vec(bsxfun(@plus, 1:nu:nu*nx, D4UsedSpecies2-1));

% A5
Iukronudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
Iukronudvuinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx, (D5UsedSpecies1-1)*nx)), [nx,nu], [2,1]);

ukronIudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
ukronIudvuinvdx_cols = vec(bsxfun(@plus, 1:nu:nu*nx, D5UsedSpecies2-1));

% Return handle
handle = @d2fdxdu;

    function val = d2fdxdu(t, x, u)
%   d2fdudx = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdx = dvdx_fun(t,x,u); % v_x
        dvdu = dvdu_fun(t,x,u); % v_u
        d2vdxdu = d2vdxdu_fun(t,x,u); % vx_x

        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1)); % v_x
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1)); % v_u
        d2vinvdxdu = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdx .* -2)), nv,1,nx), full(dvdu) .* -1), nv,nu*nx) + bsxfun(@times, (v .^ -2), reshape(d2vdxdu, nv,nu*nx) .* -1)); % v_ux

        vxinv = 1 ./ v(vxInd); % x_
        dvxinvdx = dvinvdx(vxInd,:); % x_x
        dvxinvdu = dvinvdu(vxInd,:); % x_u
        d2vxinvdxdu = d2vinvdxdu(vxInd,:); % x_xu
        
        vuinv = 1 ./ v(vuInd); % u_
        dvuinvdx = dvinvdx(vuInd,:); % u_x
        dvuinvdu = dvinvdu(vuInd,:); % u_u
        d2vuinvdxdu = d2vinvdxdu(vuInd,:); % u_uu
        
        % A2
        Ixkronxdvxinvdu = sparse(Ixkronxdvxinvdu_rows, Ixkronxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nu*nx);
        xkronIxdvxinvdu = sparse(xkronIxdvxinvdu_rows, xkronIxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies1), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nu*nx);
        
        xkronxd2vxinvdudx = sparse(nx*nx,nu*nx);
        xkronxd2vxinvdudx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdxdu(D2UsedSpecies2,:));
        
        % A3
        IukronIxvxinv = sparse(D3UsedColumns, D3UsedColumnsReverse, vxinv(D3UsedSpecies2), nx*nu,nu*nx);
        ukronIxdvxinvdu = sparse(ukronIxdvxinvdu_rows, ukronIxdvxinvdu_cols, bsxfun(@times, u(D3UsedSpecies1), dvxinvdu(D3UsedSpecies2,:)), nx*nu,nu*nx);
        Iukronxdvxinvdu = sparse(Iukronxdvxinvdx_rows, Iukronxdvxinvdx_cols, bsxfun(@times, x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:)), nx*nu,nu*nx);

        ukronxd2vxinvdudx = sparse(nx*nu,nu*nx);
        ukronxd2vxinvdudx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdxdu(D3UsedSpecies2,:));

        % A4
        IxkronIuvuinv_1   = sparse(D4UsedColumns, D4UsedColumns, vuinv(D4UsedSpecies2), nu*nx,nu*nx);
        Ixkronudvuinvdu_1 = sparse(Ixkronudvuinvdu_rows, Ixkronudvuinvdu_cols, bsxfun(@times, u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:)), nu*nx,nu*nx);
        xkronIudvuinvdx_1 = sparse(xkronIudvuinvdx_rows, xkronIudvuinvdx_cols, bsxfun(@times, x(D4UsedSpecies1), dvuinvdx(D4UsedSpecies2,:)), nu*nx,nu*nx);

        xkronud2vuinvdudx = sparse(nu*nx,nu*nx);
        xkronud2vuinvdudx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdxdu(D4UsedSpecies2,:));
        
        % A5
        Iukronudvuinvdx = sparse(Iukronudvuinvdx_rows, Iukronudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nu*nx);
        ukronIudvuinvdx = sparse(ukronIudvuinvdx_rows, ukronIudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies1), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nu*nx);
        
        ukronud2vuinvdudx = sparse(nu*nu,nu*nx);
        ukronud2vuinvdudx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdxdu(D5UsedSpecies2,:));

        val = A2 * (Ixkronxdvxinvdu + xkronIxdvxinvdu + xkronxd2vxinvdudx) ...
              + A3 * (IukronIxvxinv + ukronIxdvxinvdu + Iukronxdvxinvdu + ukronxd2vxinvdudx) ...
              + A4 * (IxkronIuvuinv_1 + Ixkronudvuinvdu_1 + xkronIudvuinvdx_1 + xkronud2vuinvdudx) ...
              + A5 * (Iukronudvuinvdx + ukronIudvuinvdx + ukronud2vuinvdudx);
        val = reshape(val, nx*nu,nx); % fu_x
    end
end

function handle = d2fdkdxHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, v_fun, dvdx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @d2fdkdx;

    function val = d2fdkdx(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        v = v_fun(t,x,u);
        dvdx = dvdx_fun(t,x,u);
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1));
        vxinv = 1 ./ v(vxInd);
        xvxinv = x .* vxinv;
        uvuinv = u ./ v(vuInd);
        dvxinvdx = dvinvdx(vxInd,:);
        dvuinvdx = dvinvdx(vuInd,:);
        
        % A2
        Ixkronxvxinv = sparse(D2UsedColumns, D2UsedSpecies1, xvxinv(D2UsedSpecies2), nx*nx,nx);
        xkronIxvxinv  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        xkronxdvxinvdx = sparse(nx*nx,nx);
        xkronxdvxinvdx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:));

        % A3
        ukronIxvxinv  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        ukronxdvxinvdx = sparse(nu*nx,nx);
        ukronxdvxinvdx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:));

        % A4
        Ixkronuvuinv = sparse(D4UsedColumns, D4UsedSpecies1, uvuinv(D4UsedSpecies2), nu*nx,nx);
        xkronudvuinvdx = sparse(nx*nu,nx);
        xkronudvuinvdx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:));
        
        % A5
        ukronudvuinvdx = sparse(nu*nu,nx);
        ukronudvuinvdx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:));

        val = dA1dk_fk_x ...
              + dA2dk_fk_xx * (Ixkronxvxinv + xkronIxvxinv + xkronxdvxinvdx) ...
              + dA3dk_fk_ux * (ukronIxvxinv + ukronxdvxinvdx) ...
              + dA4dk_fk_xu * (Ixkronuvuinv + xkronudvuinvdx) ...
              + dA5dk_fk_uu * ukronudvuinvdx; % fk_x
        val = spermute132(val, [nx,nk,nx], [nx*nx,nk]); % fx_k
    end
end

function handle = d2fdxdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, v_fun, dvdx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @d2fdkdx;

    function val = d2fdkdx(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        v = v_fun(t,x,u);
        dvdx = dvdx_fun(t,x,u);
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1));
        vxinv = 1 ./ v(vxInd);
        xvxinv = x .* vxinv;
        uvuinv = u ./ v(vuInd);
        dvxinvdx = dvinvdx(vxInd,:);
        dvuinvdx = dvinvdx(vuInd,:);
        
        % A2
        Ixkronxvxinv = sparse(D2UsedColumns, D2UsedSpecies1, xvxinv(D2UsedSpecies2), nx*nx,nx);
        xkronIxvxinv  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        xkronxdvxinvdx = sparse(nx*nx,nx);
        xkronxdvxinvdx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:));

        % A3
        ukronIxvxinv  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        ukronxdvxinvdx = sparse(nu*nx,nx);
        ukronxdvxinvdx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:));

        % A4
        Ixkronuvuinv = sparse(D4UsedColumns, D4UsedSpecies1, uvuinv(D4UsedSpecies2), nu*nx,nx);
        xkronudvuinvdx = sparse(nx*nu,nx);
        xkronudvuinvdx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:));
        
        % A5
        ukronudvuinvdx = sparse(nu*nu,nx);
        ukronudvuinvdx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:));

        val = dA1dk_fk_x ...
              + dA2dk_fk_xx * (Ixkronxvxinv + xkronIxvxinv + xkronxdvxinvdx) ...
              + dA3dk_fk_ux * (ukronIxvxinv + ukronxdvxinvdx) ...
              + dA4dk_fk_xu * (Ixkronuvuinv + xkronudvuinvdx) ...
              + dA5dk_fk_uu * ukronudvuinvdx; % fk_x
    end
end

function handle = d2fdkduHidden(dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, v_fun, dvdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA3dk_fk_ux, 1) / nx;

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        v = v_fun(t,x,u);
        dvdu = dvdu_fun(t,x,u);
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1));
        xvxinv = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvuinv = u .* vuinv;
        dvxinvdu = dvinvdu(vxInd,:);
        dvuinvdu = dvinvdu(vuInd,:);
        
        % A2
        xkronxdvxinvdu = sparse(nx*nx,nu);
        xkronxdvxinvdu(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:));
        
        % A3
        Iukronxvxinv = sparse(D3UsedColumns, D3UsedSpecies1, xvxinv(D3UsedSpecies2), nx*nu,nu);
        ukronxdvxinvdu = sparse(nu*nx,nu);
        ukronxdvxinvdu(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:));
        
        % A4
        xkronIuvuinv  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        xkronudvuinvdu = sparse(nx*nu,nu);
        xkronudvuinvdu(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:));
        
        % A5
        Iukronuvuinv = sparse(D5UsedColumns, D5UsedSpecies1, uvuinv(D5UsedSpecies2), nu*nu,nu);
        ukronIuvuinv  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        ukronudvuinvdu = sparse(nu*nu,nu);
        ukronudvuinvdu(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:));
        
        val = dA2dk_fk_xx * xkronxdvxinvdu ...
              + dA3dk_fk_ux * (Iukronxvxinv + ukronxdvxinvdu) ...
              + dA4dk_fk_xu * (xkronIuvuinv + xkronudvuinvdu) ...
              + dA5dk_fk_uu * (Iukronuvuinv + ukronIuvuinv + ukronudvuinvdu) ...
              + dA6dk_fk_u; % fk_u
        val = spermute132(val, [nx,nk,nu], [nx*nu,nk]); % fu_k
    end
end

function handle = d2fdudkHidden(dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, v_fun, dvdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA3dk_fk_ux, 1) / nx;

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        v = v_fun(t,x,u);
        dvdu = dvdu_fun(t,x,u);
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1));
        xvxinv = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvuinv = u .* vuinv;
        dvxinvdu = dvinvdu(vxInd,:);
        dvuinvdu = dvinvdu(vuInd,:);
        
        % A2
        xkronxdvxinvdu = sparse(nx*nx,nu);
        xkronxdvxinvdu(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:));
        
        % A3
        Iukronxvxinv = sparse(D3UsedColumns, D3UsedSpecies1, xvxinv(D3UsedSpecies2), nx*nu,nu);
        ukronxdvxinvdu = sparse(nu*nx,nu);
        ukronxdvxinvdu(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:));
        
        % A4
        xkronIuvuinv  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        xkronudvuinvdu = sparse(nx*nu,nu);
        xkronudvuinvdu(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:));
        
        % A5
        Iukronuvuinv = sparse(D5UsedColumns, D5UsedSpecies1, uvuinv(D5UsedSpecies2), nu*nu,nu);
        ukronIuvuinv  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        ukronudvuinvdu = sparse(nu*nu,nu);
        ukronudvuinvdu(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:));
        
        val = dA2dk_fk_xx * xkronxdvxinvdu ...
              + dA3dk_fk_ux * (Iukronxvxinv + ukronxdvxinvdu) ...
              + dA4dk_fk_xu * (xkronIuvuinv + xkronudvuinvdu) ...
              + dA5dk_fk_uu * (Iukronuvuinv + ukronIuvuinv + ukronudvuinvdu) ...
              + dA6dk_fk_u; % fk_u
    end
end

function handle = rHidden(D1, D2, D3, D4, D5, D6, d, v_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @r;

    function val = r(t, x, u)
%   r = D1*x + D2*(x kron x/vx) + D3*(u kron x/vx) + D4*(x kron u/vu) + D5*(u kron u/vu) + D6*u + d
        % Compartment column
        v = v_fun(t,x,u);
        xvxinv = x ./ v(vxInd);
        uvuinv = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronx = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvxinv(D2UsedSpecies2), nx*nx,1);
        ukronx = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvxinv(D3UsedSpecies2), nx*nu,1);
        xkronu = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvuinv(D4UsedSpecies2), nu*nx,1);
        ukronu = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvuinv(D5UsedSpecies2), nu*nu,1);
        
        val = D1 * x + D2 * xkronx + D3 * ukronx + D4 * xkronu + D5 * ukronu + D6 * u + d; % r_
    end
end

function handle = drdxHidden(D1, D2, D3, D4, D5, v_fun, dvdx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @drdx;

    function val = drdx(t, x, u)
%   drdx = D1 + D2*(Ix kron x/vx) + D2*(x kron diag(1/vx)) + D3*(u kron diag(1/vx))
        % Compartment column
        v = v_fun(t,x,u);
        dvdx = dvdx_fun(t,x,u);
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1));
        vxinv = 1 ./ v(vxInd);
        xvxinv = x .* vxinv;
        uvuinv = u ./ v(vuInd);
        dvxinvdx = dvinvdx(vxInd,:);
        dvuinvdx = dvinvdx(vuInd,:);
        
        % A2
        Ixkronxvxinv   = sparse(D2UsedColumns, D2UsedSpecies1, xvxinv(D2UsedSpecies2), nx*nx,nx);
        xkronIxvxinv   = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        xkronxdvxinvdx = sparse(nx*nx,nx);
        xkronxdvxinvdx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:));
        
        % A3
        ukronIxvxinv   = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        ukronxdvxinvdx = sparse(nu*nx,nx);
        ukronxdvxinvdx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:));
        
        % A4
        Ixkronuvuinv   = sparse(D4UsedColumns, D4UsedSpecies1, uvuinv(D4UsedSpecies2), nx*nu,nx);
        xkronudvuinvdx = sparse(nx*nu,nx);
        xkronudvuinvdx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:));
        
        % A5
        ukronudvuinvdx = sparse(nu*nu,nx);
        ukronudvuinvdx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:));
        
        val = D1 + D2 * (Ixkronxvxinv + xkronIxvxinv + xkronxdvxinvdx) + D3 * (ukronIxvxinv + ukronxdvxinvdx) + D4 * (Ixkronuvuinv + xkronudvuinvdx) + D5 * ukronudvuinvdx; % f_x
    end
end

function handle = drduHidden(D2, D3, D4, D5, D6, v_fun, dvdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @drdu;

    function val = drdu(t, x, u)
%   drdu = D3*(Iu kron x/vx) + D4*(x kron diag(1/vu)) + D5*(Iu kron u/vu) + D5*(u kron diag(1/vu)) + D6
        % Compartment column
        v = v_fun(t,x,u);
        dvdu = dvdu_fun(t,x,u);
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1));
        xvxinv = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvuinv = u .* vuinv;
        dvxinvdu = dvinvdu(vxInd,:);
        dvuinvdu = dvinvdu(vuInd,:);
        
        % A2
        xkronxdvxinvdu = sparse(nx*nx,nu);
        xkronxdvxinvdu(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:));

        % A3
        Iukronxvxinv   = sparse(D3UsedColumns, D3UsedSpecies1, xvxinv(D3UsedSpecies2), nx*nu,nu);
        ukronxdvxinvdu = sparse(nu*nx,nu);
        ukronxdvxinvdu(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:));
        
        % A4
        xkronIuvuinv   = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        xkronudvuinvdu = sparse(nx*nu,nu);
        xkronudvuinvdu(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:));

        % A5
        Iukronuvuinv   = sparse(D5UsedColumns, D5UsedSpecies1, uvuinv(D5UsedSpecies2), nu*nu,nu);
        ukronIuvuinv   = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        ukronudvuinvdu = sparse(nu*nu,nu);
        ukronudvuinvdu(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:));
        
        val = D2 * xkronxdvxinvdu + D3 * (Iukronxvxinv + ukronxdvxinvdu) + D4 * (xkronIuvuinv + xkronudvuinvdu) + D5 * (Iukronuvuinv + ukronIuvuinv + ukronudvuinvdu) + D6; % f_u
    end
end

function handle = drdkHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, dD5dk_rk_uu, dD6dk_rk_u, dddk, v_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
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
        v = v_fun(t,x,u);
        xvxinv = x ./ v(vxInd);
        uvuinv = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        xkronxinv = sparse(D2UsedColumns, ones(numel(D2UsedColumns),1), x(D2UsedSpecies1) .* xvxinv(D2UsedSpecies2), nx*nx,1);
        ukronxinv = sparse(D3UsedColumns, ones(numel(D3UsedColumns),1), u(D3UsedSpecies1) .* xvxinv(D3UsedSpecies2), nx*nu,1);
        xkronuinv = sparse(D4UsedColumns, ones(numel(D4UsedColumns),1), x(D4UsedSpecies1) .* uvuinv(D4UsedSpecies2), nu*nx,1);
        ukronuinv = sparse(D5UsedColumns, ones(numel(D5UsedColumns),1), u(D5UsedSpecies1) .* uvuinv(D5UsedSpecies2), nu*nu,1);
        
        val = dD1dk_rk_x * sparse(x) + dD2dk_rk_xx * xkronxinv + dD3dk_rk_ux * ukronxinv + dD4dk_rk_xu * xkronuinv + dD5dk_rk_uu * ukronuinv + dD6dk_rk_u * sparse(u); % fk_
        val = reshape(val, nr,nk) + dddk; % r_k
    end
end

function handle = d2rdx2Hidden(D2, D3, D4, D5, v_fun, dvdx_fun, d2vdx2_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nr = size(D2, 1);
nx = numel(vxInd);
nu = numel(vuInd);

% A2
D2UsedColumnsReverse = linearindexpermute(D2UsedColumns, [nx,nx], [2,1]);

Ixkronxdvxinvdx_rows = repmat(D2UsedColumns,[nx,1]);
Ixkronxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx, (D2UsedSpecies1-1)*nx));
 
xkronIxdvxinvdx_rows = Ixkronxdvxinvdx_rows;
xkronIxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx:nx*nx, D2UsedSpecies2-1));

% A3
ukronIxdvxinvdx_rows = repmat(D3UsedColumns,[nx,1]);
ukronIxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx:nx*nx, D3UsedSpecies2-1));

% A4
Ixkronudvuinvdx_rows = repmat(D4UsedColumns,[nx,1]);
Ixkronudvuinvdx_cols = vec(bsxfun(@plus, 1:nx, (D4UsedSpecies1-1)*nx));

% Return handle
handle = @d2rdx2;

    function val = d2rdx2(t, x, u)
%   d2rdx2 = 2*D2*(Ix kron diag(1/vx))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdx = dvdx_fun(t,x,u); % v_x
        d2vdx2 = d2vdx2_fun(t,x,u); % vx_x
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1)); % v_x
        d2vinvdx2 = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdx .* -2)), nv,1,nx), full(dvdx) .* -1), nv,nx*nx) + bsxfun(@times, (v .^ -2), reshape(d2vdx2, nv,nx*nx) .* -1)); % v_xx
        vxinv = 1 ./ v(vxInd); % x_
        dvxinvdx = dvinvdx(vxInd,:); % x_x
        d2vxinvdx2 = d2vinvdx2(vxInd,:); % x_xx
        dvuinvdx = dvinvdx(vuInd,:); % u_x
        d2vuinvdx2 = d2vinvdx2(vuInd,:); % u_xx
        
        % A2
        %IxkronIxvxinv = sparse([D2UsedColumns; D2UsedColumns], [D2UsedColumns; D2UsedColumnsReverse], [vxinv(D2UsedSpecies2); vxinv(D2UsedSpecies2)], nx*nx,nx*nx);
        IxkronIxvxinv_1 = sparse(D2UsedColumns, D2UsedColumnsReverse, vxinv(D2UsedSpecies2), nx*nx,nx*nx);
        IxkronIxvxinv_2 = sparse(D2UsedColumns, D2UsedColumns, vxinv(D2UsedSpecies2), nx*nx,nx*nx);
        
        %Ixkronxdvxinvdx = kron(eye(nx), diag(x) * dvxinvdx);
        Ixkronxdvxinvdx_1 = sparse(Ixkronxdvxinvdx_rows, Ixkronxdvxinvdx_cols, bsxfun(@times, x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:)), nx*nx,nx*nx);
        Ixkronxdvxinvdx_2 = Ixkronxdvxinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        %xkronIxdvxinvdx = reshape(kron(x, reshape(bsxfun(@times, eye(nx), reshape(full(dvxinvdx),[nx,1,nx])), [nx,nx*nx])), [nx*nx,nx*nx]);
        xkronIxdvxinvdx_1 = sparse(xkronIxdvxinvdx_rows, xkronIxdvxinvdx_cols, bsxfun(@times, x(D2UsedSpecies1), dvxinvdx(D2UsedSpecies2,:)), nx*nx,nx*nx);
        xkronIxdvxinvdx_2 = xkronIxdvxinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        %xkronxd2vxinvdx2 = kron(x, diag(x) * d2vxinvdx2);
        xkronxd2vxinvdx2 = sparse(nx*nx,nx*nx);
        xkronxd2vxinvdx2(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdx2(D2UsedSpecies2,:));

        % A3
        ukronIxdvxinvdx_1 = sparse(ukronIxdvxinvdx_rows, ukronIxdvxinvdx_cols, bsxfun(@times, u(D3UsedSpecies1), dvxinvdx(D3UsedSpecies2,:)), nx*nu,nx*nx);
        ukronIxdvxinvdx_2 = ukronIxdvxinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        ukronxd2vxinvdx2 = sparse(nx*nu,nx*nx);
        ukronxd2vxinvdx2(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdx2(D3UsedSpecies2,:));

        % A4
        Ixkronudvuinvdx_1 = sparse(Ixkronudvuinvdx_rows, Ixkronudvuinvdx_cols, bsxfun(@times, u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:)), nu*nx,nx*nx);
        Ixkronudvuinvdx_2 = Ixkronudvuinvdx_1(:,reshape(1:(nx*nx), nx,nx)');
        
        xkronud2vuinvdx2 = sparse(nu*nx,nx*nx);
        xkronud2vuinvdx2(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdx2(D4UsedSpecies2,:));
        
        % A5
        ukronud2vuinvdx2 = sparse(nu*nu,nx*nx);
        ukronud2vuinvdx2(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdx2(D5UsedSpecies2,:));

        val = D2 * (IxkronIxvxinv_1 + Ixkronxdvxinvdx_1 + IxkronIxvxinv_2 + xkronIxdvxinvdx_1 + Ixkronxdvxinvdx_2 + xkronIxdvxinvdx_2 + xkronxd2vxinvdx2)...
              + D3 * (ukronIxdvxinvdx_1 + ukronIxdvxinvdx_2 + ukronxd2vxinvdx2)...
              + D4 * (Ixkronudvuinvdx_1 + Ixkronudvuinvdx_2 + xkronud2vuinvdx2)...
              + D5 * ukronud2vuinvdx2; % r_xx
        val = reshape(val, nr*nx, nx); % rx_x
    end
end

function handle = d2rdu2Hidden(D2, D3, D4, D5, v_fun, dvdu_fun, d2vdu2_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nr = size(D2, 1);
nx = numel(vxInd);
nu = numel(vuInd);

% A3
Iukronxdvxinvdu_rows = repmat(D3UsedColumns,[nu,1]);
Iukronxdvxinvdu_cols = vec(bsxfun(@plus, 1:nu, (D3UsedSpecies1-1)*nu));

% A4
xkronIudvuinvdu_rows = repmat(D4UsedColumns,[nu,1]);
xkronIudvuinvdu_cols = vec(bsxfun(@plus, 1:nu:nu*nu, D4UsedSpecies2-1));

% A5
[u1ind,u2ind] = ind2sub([nu,nu], D5UsedColumns);
D5UsedColumnsReverse = sub2ind([nu, nu], u2ind, u1ind);

Iukronudvuinvdu_rows = repmat(D5UsedColumns,[nu,1]);
Iukronudvuinvdu_cols = vec(bsxfun(@plus, 1:nu, (D5UsedSpecies1-1)*nu));
 
ukronIudvuinvdu_rows = Iukronudvuinvdu_rows;
ukronIudvuinvdu_cols = vec(bsxfun(@plus, 1:nu:nu*nu, D5UsedSpecies2-1));

% Return handle
handle = @d2rdu2;

    function val = d2rdu2(t, x, u)
%   d2fdu2 = 2*A5*(Iu kron diag(1/vu))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdu = dvdu_fun(t,x,u); % v_u
        d2vdu2 = d2vdu2_fun(t,x,u); % vu_u
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1)); % v_u
        d2vinvdu2 = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdu .* -2)), nv,1,nu), full(dvdu) .* -1), nv,nu*nu) + bsxfun(@times, (v .^ -2), reshape(d2vdu2, nv,nu*nu) .* -1)); % v_uu
        dvxinvdu = dvinvdu(vxInd,:); % x_u
        d2vxinvdu2 = d2vinvdu2(vxInd,:); % x_uu
        vuinv = 1 ./ v(vuInd); % u_
        dvuinvdu = dvinvdu(vuInd,:); % u_u
        d2vuinvdu2 = d2vinvdu2(vuInd,:); % u_uu
        
        % A2
        xkronxd2vxinvdu2 = sparse(nx*nx,nu*nu);
        xkronxd2vxinvdu2(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdu2(D2UsedSpecies2,:));

        % A3
        Iukronxdvxinvdu_1 = sparse(Iukronxdvxinvdu_rows, Iukronxdvxinvdu_cols, bsxfun(@times, x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:)), nx*nu,nu*nu);
        Iukronxdvxinvdu_2 = Iukronxdvxinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        ukronxd2vxinvdu2 = sparse(nx*nu,nu*nu);
        ukronxd2vxinvdu2(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdu2(D3UsedSpecies2,:));
        
        % A4
        xkronIudvuinvdu_1 = sparse(xkronIudvuinvdu_rows, xkronIudvuinvdu_cols, bsxfun(@times, x(D4UsedSpecies1), dvuinvdu(D4UsedSpecies2,:)), nu*nx,nu*nu);
        xkronIudvuinvdu_2 = xkronIudvuinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        xkronud2vuinvdu2 = sparse(nu*nx,nu*nu);
        xkronud2vuinvdu2(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdu2(D4UsedSpecies2,:));

        % A5
        IukronIuvuinv_1 = sparse(D5UsedColumns, D5UsedColumnsReverse, vuinv(D5UsedSpecies2), nu*nu,nu*nu);
        IukronIuvuinv_2 = sparse(D5UsedColumns, D5UsedColumns, vuinv(D5UsedSpecies2), nu*nu,nu*nu);
        
        %Iukronudvuinvdu = kron(eye(nu), diag(u) * dvuinvdu);
        Iukronudvuinvdu_1 = sparse(Iukronudvuinvdu_rows, Iukronudvuinvdu_cols, bsxfun(@times, u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:)), nu*nu,nu*nu);
        Iukronudvuinvdu_2 = Iukronudvuinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        %ukronIudvuinvdu = reshape(kron(u, reshape(bsxfun(@times, eye(nu), reshape(full(dvuinvdu),[nu,1,nu])), [nu,nu*nu])), [nu*nu,nu*nu]);
        ukronIudvuinvdu_1 = sparse(ukronIudvuinvdu_rows, ukronIudvuinvdu_cols, bsxfun(@times, u(D5UsedSpecies1), dvuinvdu(D5UsedSpecies2,:)), nu*nu,nu*nu);
        ukronIudvuinvdu_2 = ukronIudvuinvdu_1(:,reshape(1:(nu*nu), nu,nu)');
        
        %ukronud2vuinvdu2 = kron(u, diag(u) * d2vuinvdu2);
        ukronud2vuinvdu2 = sparse(nu*nu,nu*nu);
        ukronud2vuinvdu2(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdu2(D5UsedSpecies2,:));

        val = D2 * xkronxd2vxinvdu2...
              + D3 * (Iukronxdvxinvdu_1 + Iukronxdvxinvdu_2 + ukronxd2vxinvdu2)...
              + D4 * (xkronIudvuinvdu_1 + xkronIudvuinvdu_2 + xkronud2vuinvdu2)...
              + D5 * (IukronIuvuinv_1 + Iukronudvuinvdu_1 + IukronIuvuinv_2 + ukronIudvuinvdu_1 + Iukronudvuinvdu_2 + ukronIudvuinvdu_2 + ukronud2vuinvdu2); % r_uu
        val = reshape(val, nr*nu, nu); % ru_u
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

function handle = d2rdudxHidden(D2, D3, D4, D5, v_fun, dvdx_fun, dvdu_fun, d2vdudx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nr = size(D2, 1);
nx = numel(vxInd);
nu = numel(vuInd);

% A2
Ixkronxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
Ixkronxdvxinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu, (D2UsedSpecies1-1)*nu)), [nu,nx], [2,1]);

xkronIxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
xkronIxdvxinvdu_cols = vec(bsxfun(@plus, 1:nx:nx*nu, D2UsedSpecies2-1));

% A3
Iukronxdvxinvdx_rows = repmat(D3UsedColumns,[nx,1]);
Iukronxdvxinvdx_cols = vec(bsxfun(@plus, 1:nx, (D3UsedSpecies1-1)*nx));

ukronIxdvxinvdu_rows = repmat(D3UsedColumns,[nu,1]);
ukronIxdvxinvdu_cols = vec(bsxfun(@plus, 1:nx:nx*nu, D3UsedSpecies2-1));

% A4
D4UsedColumnsReverse = linearindexpermute(D4UsedColumns, [nu, nx], [2,1]);

Ixkronudvuinvdu_rows = repmat(D4UsedColumns,[nu,1]);
Ixkronudvuinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu, (D4UsedSpecies1-1)*nu)), [nu,nx], [2,1]);

xkronIudvuinvdx_rows = repmat(D4UsedColumns,[nx,1]);
xkronIudvuinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu:nu*nx, D4UsedSpecies2-1)), [nu,nx], [2,1]);

% A5
Iukronudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
Iukronudvuinvdx_cols = vec(bsxfun(@plus, 1:nx, (D5UsedSpecies1-1)*nx));

ukronIudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
ukronIudvuinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nu:nu*nx, D5UsedSpecies2-1)), [nu,nx], [2,1]);

% Return handle
handle = @d2rdudx;

    function val = d2rdudx(t, x, u)
%   d2fdudx = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdx = dvdx_fun(t,x,u); % v_x
        dvdu = dvdu_fun(t,x,u); % v_u
        d2vdudx = d2vdudx_fun(t,x,u); % vx_x

        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1)); % v_x
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1)); % v_u
        d2vinvdudx = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdu .* -2)), nv,1,nu), full(dvdx) .* -1), nv,nx*nu) + bsxfun(@times, (v .^ -2), reshape(d2vdudx, nv,nx*nu) .* -1)); % v_xu

        vxinv = 1 ./ v(vxInd); % x_
        dvxinvdx = dvinvdx(vxInd,:); % x_x
        dvxinvdu = dvinvdu(vxInd,:); % x_u
        d2vxinvdudx = d2vinvdudx(vxInd,:); % x_xu
        
        vuinv = 1 ./ v(vuInd); % u_
        dvuinvdx = dvinvdx(vuInd,:); % u_x
        dvuinvdu = dvinvdu(vuInd,:); % u_u
        d2vuinvdudx = d2vinvdudx(vuInd,:); % u_uu
        
        % A2 derivatives
        Ixkronxdvxinvdu = sparse(Ixkronxdvxinvdu_rows, Ixkronxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nx*nu);
        xkronIxdvxinvdu = sparse(xkronIxdvxinvdu_rows, xkronIxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies1), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nx*nu);
        
        xkronxd2vxinvdudx = sparse(nx*nx,nx*nu);
        xkronxd2vxinvdudx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdudx(D2UsedSpecies2,:));
        
        % A3 derivatives
        IukronIxvxinv   = sparse(D3UsedColumns, D3UsedColumns, vxinv(D3UsedSpecies2), nx*nu,nx*nu);
        ukronIxdvxinvdu = sparse(ukronIxdvxinvdu_rows, ukronIxdvxinvdu_cols, bsxfun(@times, u(D3UsedSpecies1), dvxinvdu(D3UsedSpecies2,:)), nx*nu,nx*nu);
        Iukronxdvxinvdu = sparse(Iukronxdvxinvdx_rows, Iukronxdvxinvdx_cols, bsxfun(@times, x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:)), nx*nu,nx*nu);

        ukronxd2vxinvdudx = sparse(nx*nu,nx*nu);
        ukronxd2vxinvdudx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdudx(D3UsedSpecies2,:));

        % A4 derivatives
        IxkronIuvuinv_1   = sparse(D4UsedColumns, D4UsedColumnsReverse, vuinv(D4UsedSpecies2), nu*nx,nx*nu);
        Ixkronudvuinvdu_1 = sparse(Ixkronudvuinvdu_rows, Ixkronudvuinvdu_cols, bsxfun(@times, u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:)), nu*nx,nx*nu);
        xkronIudvuinvdx_1 = sparse(xkronIudvuinvdx_rows, xkronIudvuinvdx_cols, bsxfun(@times, x(D4UsedSpecies1), dvuinvdx(D4UsedSpecies2,:)), nu*nx,nx*nu);

        xkronud2vuinvdudx = sparse(nu*nx,nx*nu);
        xkronud2vuinvdudx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdudx(D4UsedSpecies2,:));
        
        % A5 derivatives
        Iukronudvuinvdx = sparse(Iukronudvuinvdx_rows, Iukronudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nx*nu);
        ukronIudvuinvdx = sparse(ukronIudvuinvdx_rows, ukronIudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies1), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nx*nu);
        
        ukronud2vuinvdudx = sparse(nu*nu,nx*nu);
        ukronud2vuinvdudx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdudx(D5UsedSpecies2,:));

        val = D2 * (Ixkronxdvxinvdu + xkronIxdvxinvdu + xkronxd2vxinvdudx) ...
              + D3 * (IukronIxvxinv + ukronIxdvxinvdu + Iukronxdvxinvdu + ukronxd2vxinvdudx) ...
              + D4 * (IxkronIuvuinv_1 + Ixkronudvuinvdu_1 + xkronIudvuinvdx_1 + xkronud2vuinvdudx) ...
              + D5 * (Iukronudvuinvdx + ukronIudvuinvdx + ukronud2vuinvdudx);
        val = reshape(val, nr*nx,nu); % rx_u
    end
end

function handle = d2rdxduHidden(D2, D3, D4, D5, v_fun, dvdx_fun, dvdu_fun, d2vdxdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nr = size(D2, 1);
nx = numel(vxInd);
nu = numel(vuInd);

% A2
Ixkronxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
Ixkronxdvxinvdu_cols = vec(bsxfun(@plus, 1:nu, (D2UsedSpecies1-1)*nu));

xkronIxdvxinvdu_rows = repmat(D2UsedColumns,[nu,1]);
xkronIxdvxinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx:nx*nu, D2UsedSpecies2-1)), [nx,nu], [2,1]);

% A3
D3UsedColumnsReverse = linearindexpermute(D3UsedColumns, [nx,nu], [2,1]);

Iukronxdvxinvdx_rows = repmat(D3UsedColumns,[nx,1]);
Iukronxdvxinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx, (D3UsedSpecies1-1)*nx)), [nx,nu], [2,1]);

ukronIxdvxinvdu_rows = repmat(D3UsedColumns,[nu,1]);
ukronIxdvxinvdu_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx:nx*nu, D3UsedSpecies2-1)), [nx,nu], [2,1]);

% A4
Ixkronudvuinvdu_rows = repmat(D4UsedColumns,[nu,1]);
Ixkronudvuinvdu_cols = vec(bsxfun(@plus, 1:nu, (D4UsedSpecies1-1)*nu));

xkronIudvuinvdx_rows = repmat(D4UsedColumns,[nx,1]);
xkronIudvuinvdx_cols = vec(bsxfun(@plus, 1:nu:nu*nx, D4UsedSpecies2-1));

% A5
Iukronudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
Iukronudvuinvdx_cols = linearindexpermute(vec(bsxfun(@plus, 1:nx, (D5UsedSpecies1-1)*nx)), [nx,nu], [2,1]);

ukronIudvuinvdx_rows = repmat(D5UsedColumns,[nx,1]);
ukronIudvuinvdx_cols = vec(bsxfun(@plus, 1:nu:nu*nx, D5UsedSpecies2-1));

% Return handle
handle = @d2rdudx;

    function val = d2rdudx(t, x, u)
%   d2fdudx = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1))
        % Compartment column
        v = v_fun(t,x,u); % v_
        nv = numel(v);
        dvdx = dvdx_fun(t,x,u); % v_x
        dvdu = dvdu_fun(t,x,u); % v_u
        d2vdxdu = d2vdxdu_fun(t,x,u); % vx_x

        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1)); % v_x
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1)); % v_u
        d2vinvdxdu = sparse(reshape(bsxfun(@times, reshape(full(bsxfun(@times, v .^ -3, dvdx .* -2)), nv,1,nx), full(dvdu) .* -1), nv,nu*nx) + bsxfun(@times, (v .^ -2), reshape(d2vdxdu, nv,nu*nx) .* -1)); % v_ux

        vxinv = 1 ./ v(vxInd); % x_
        dvxinvdx = dvinvdx(vxInd,:); % x_x
        dvxinvdu = dvinvdu(vxInd,:); % x_u
        d2vxinvdxdu = d2vinvdxdu(vxInd,:); % x_xu
        
        vuinv = 1 ./ v(vuInd); % u_
        dvuinvdx = dvinvdx(vuInd,:); % u_x
        dvuinvdu = dvinvdu(vuInd,:); % u_u
        d2vuinvdxdu = d2vinvdxdu(vuInd,:); % u_uu
        
        % A2
        Ixkronxdvxinvdu = sparse(Ixkronxdvxinvdu_rows, Ixkronxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nu*nx);
        xkronIxdvxinvdu = sparse(xkronIxdvxinvdu_rows, xkronIxdvxinvdu_cols, bsxfun(@times, x(D2UsedSpecies1), dvxinvdu(D2UsedSpecies2,:)), nx*nx,nu*nx);
        
        xkronxd2vxinvdudx = sparse(nx*nx,nx*nu);
        xkronxd2vxinvdudx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), d2vxinvdxdu(D2UsedSpecies2,:));
        
        % A3
        IukronIxvxinv = sparse(D3UsedColumns, D3UsedColumnsReverse, vxinv(D3UsedSpecies2), nx*nu,nu*nx);
        ukronIxdvxinvdu = sparse(ukronIxdvxinvdu_rows, ukronIxdvxinvdu_cols, bsxfun(@times, u(D3UsedSpecies1), dvxinvdu(D3UsedSpecies2,:)), nx*nu,nu*nx);
        Iukronxdvxinvdu = sparse(Iukronxdvxinvdx_rows, Iukronxdvxinvdx_cols, bsxfun(@times, x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:)), nx*nu,nu*nx);

        ukronxd2vxinvdudx = sparse(nx*nu,nx*nu);
        ukronxd2vxinvdudx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), d2vxinvdxdu(D3UsedSpecies2,:));

        % A4
        IxkronIuvuinv_1   = sparse(D4UsedColumns, D4UsedColumns, vuinv(D4UsedSpecies2), nu*nx,nx*nu);
        Ixkronudvuinvdu_1 = sparse(Ixkronudvuinvdu_rows, Ixkronudvuinvdu_cols, bsxfun(@times, u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:)), nu*nx,nu*nx);
        xkronIudvuinvdx_1 = sparse(xkronIudvuinvdx_rows, xkronIudvuinvdx_cols, bsxfun(@times, x(D4UsedSpecies1), dvuinvdx(D4UsedSpecies2,:)), nu*nx,nu*nx);

        xkronud2vuinvdudx = sparse(nu*nx,nx*nu);
        xkronud2vuinvdudx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), d2vuinvdxdu(D4UsedSpecies2,:));
        
        % A5
        Iukronudvuinvdx = sparse(Iukronudvuinvdx_rows, Iukronudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nu*nx);
        ukronIudvuinvdx = sparse(ukronIudvuinvdx_rows, ukronIudvuinvdx_cols, bsxfun(@times, u(D5UsedSpecies1), dvuinvdx(D5UsedSpecies2,:)), nu*nu,nu*nx);
        
        ukronud2vuinvdudx = sparse(nu*nu,nx*nu);
        ukronud2vuinvdudx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), d2vuinvdxdu(D5UsedSpecies2,:));

        val = D2 * (Ixkronxdvxinvdu + xkronIxdvxinvdu + xkronxd2vxinvdudx) ...
              + D3 * (IukronIxvxinv + ukronIxdvxinvdu + Iukronxdvxinvdu + ukronxd2vxinvdudx) ...
              + D4 * (IxkronIuvuinv_1 + Ixkronudvuinvdu_1 + xkronIudvuinvdx_1 + xkronud2vuinvdudx) ...
              + D5 * (Iukronudvuinvdx + ukronIudvuinvdx + ukronud2vuinvdudx);
        val = reshape(val, nr*nu,nx); % ru_x
    end
end

function handle = d2rdkdxHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, dD5dk_rk_uu, v_fun, dvdx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd, nr)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dD1dk_rk_x, 1) / nr;

% Return handle
handle = @d2rdkdx;

    function val = d2rdkdx(t, x, u)
%   d2rdxdk = dD1dk + dD2dk*(Ix kron x/vx) + dD2dk*(x kron diag(1/vx)) + dD3dk*(u kron diag(1/vx)) + dD4dk*(Ix kron u/vu)
        % Compartment column
        v = v_fun(t,x,u);
        dvdx = dvdx_fun(t,x,u);
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1));
        vxinv = 1 ./ v(vxInd);
        xvxinv = x .* vxinv;
        uvuinv = u ./ v(vuInd);
        dvxinvdx = dvinvdx(vxInd,:);
        dvuinvdx = dvinvdx(vuInd,:);
        
        % A2
        Ixkronxvxinv = sparse(D2UsedColumns, D2UsedSpecies1, xvxinv(D2UsedSpecies2), nx*nx,nx);
        xkronIxvxinv  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        xkronxdvxinvdx = sparse(nx*nx,nx);
        xkronxdvxinvdx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:));

        % A3
        ukronIxvxinv  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        ukronxdvxinvdx = sparse(nu*nx,nx);
        ukronxdvxinvdx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:));

        % A4
        Ixkronuvuinv = sparse(D4UsedColumns, D4UsedSpecies1, uvuinv(D4UsedSpecies2), nu*nx,nx);
        xkronudvuinvdx = sparse(nx*nu,nx);
        xkronudvuinvdx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:));
        
        % A5
        ukronudvuinvdx = sparse(nu*nu,nx);
        ukronudvuinvdx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:));

        val = dD1dk_rk_x ...
              + dD2dk_rk_xx * (Ixkronxvxinv + xkronIxvxinv + xkronxdvxinvdx) ...
              + dD3dk_rk_ux * (ukronIxvxinv + ukronxdvxinvdx) ...
              + dD4dk_rk_xu * (Ixkronuvuinv + xkronudvuinvdx) ...
              + dD5dk_rk_uu * ukronudvuinvdx; % rk_x
        val = spermute132(val, [nr,nk,nx], [nr*nx,nk]); % rx_k
    end
end

function handle = d2rdxdkHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, dD5dk_rk_uu, v_fun, dvdx_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2rdxdk;

    function val = d2rdxdk(t, x, u)
%   d2rdxdk = dD1dk + dD2dk*(Ix kron x/vx) + dD2dk*(x kron diag(1/vx)) + dD3dk*(u kron diag(1/vx)) + dD4dk*(Ix kron u/vu)
        % Compartment column
        v = v_fun(t,x,u);
        dvdx = dvdx_fun(t,x,u);
        dvinvdx = bsxfun(@times, (v .^ -2), (dvdx .* -1));
        vxinv = 1 ./ v(vxInd);
        xvxinv = x .* vxinv;
        uvuinv = u ./ v(vuInd);
        dvxinvdx = dvinvdx(vxInd,:);
        dvuinvdx = dvinvdx(vuInd,:);
        
        % A2
        Ixkronxvxinv = sparse(D2UsedColumns, D2UsedSpecies1, xvxinv(D2UsedSpecies2), nx*nx,nx);
        xkronIxvxinv  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vxinv(D2UsedSpecies2), nx*nx,nx);
        xkronxdvxinvdx = sparse(nx*nx,nx);
        xkronxdvxinvdx(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdx(D2UsedSpecies2,:));

        % A3
        ukronIxvxinv  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vxinv(D3UsedSpecies2), nx*nu,nx);
        ukronxdvxinvdx = sparse(nu*nx,nx);
        ukronxdvxinvdx(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdx(D3UsedSpecies2,:));

        % A4
        Ixkronuvuinv = sparse(D4UsedColumns, D4UsedSpecies1, uvuinv(D4UsedSpecies2), nu*nx,nx);
        xkronudvuinvdx = sparse(nx*nu,nx);
        xkronudvuinvdx(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdx(D4UsedSpecies2,:));
        
        % A5
        ukronudvuinvdx = sparse(nu*nu,nx);
        ukronudvuinvdx(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdx(D5UsedSpecies2,:));

        val = dD1dk_rk_x ...
              + dD2dk_rk_xx * (Ixkronxvxinv + xkronIxvxinv + xkronxdvxinvdx) ...
              + dD3dk_rk_ux * (ukronIxvxinv + ukronxdvxinvdx) ...
              + dD4dk_rk_xu * (Ixkronuvuinv + xkronudvuinvdx) ...
              + dD5dk_rk_uu * ukronudvuinvdx; % rk_x
    end
end

function handle = d2rdkduHidden(dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, dD5dk_rk_uu, dD6dk_rk_u, v_fun, dvdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd, nr)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dD3dk_rk_ux, 1) / nr;

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        v = v_fun(t,x,u);
        dvdu = dvdu_fun(t,x,u);
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1));
        xvxinv = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvuinv = u .* vuinv;
        dvxinvdu = dvinvdu(vxInd,:);
        dvuinvdu = dvinvdu(vuInd,:);
        
        % A2
        xkronxdvxinvdu = sparse(nx*nx,nu);
        xkronxdvxinvdu(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:));
        
        % A3
        Iukronxvxinv = sparse(D3UsedColumns, D3UsedSpecies1, xvxinv(D3UsedSpecies2), nx*nu,nu);
        ukronxdvxinvdu = sparse(nu*nx,nu);
        ukronxdvxinvdu(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:));
        
        % A4
        xkronIuvuinv  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        xkronudvuinvdu = sparse(nx*nu,nu);
        xkronudvuinvdu(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:));
        
        % A5
        Iukronuvuinv = sparse(D5UsedColumns, D5UsedSpecies1, uvuinv(D5UsedSpecies2), nu*nu,nu);
        ukronIuvuinv  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        ukronudvuinvdu = sparse(nu*nu,nu);
        ukronudvuinvdu(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:));
        
        val = dD2dk_rk_xx * xkronxdvxinvdu ...
              + dD3dk_rk_ux * (Iukronxvxinv + ukronxdvxinvdu) ...
              + dD4dk_rk_xu * (xkronIuvuinv + xkronudvuinvdu) ...
              + dD5dk_rk_uu * (Iukronuvuinv + ukronIuvuinv + ukronudvuinvdu) ...
              + dD6dk_rk_u; % fk_u
        val = spermute132(val, [nr,nk,nu], [nr*nu,nk]); % ru_k
    end
end

function handle = d2rdudkHidden(dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, dD5dk_rk_uu, dD6dk_rk_u, v_fun, dvdu_fun, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        v = v_fun(t,x,u);
        dvdu = dvdu_fun(t,x,u);
        dvinvdu = bsxfun(@times, (v .^ -2), (dvdu .* -1));
        xvxinv = x ./ v(vxInd);
        vuinv = 1 ./ v(vuInd);
        uvuinv = u .* vuinv;
        dvxinvdu = dvinvdu(vxInd,:);
        dvuinvdu = dvinvdu(vuInd,:);
        
        % A2
        xkronxdvxinvdu = sparse(nx*nx,nu);
        xkronxdvxinvdu(D2UsedColumns,:) = bsxfun(@times, x(D2UsedSpecies1) .* x(D2UsedSpecies2), dvxinvdu(D2UsedSpecies2,:));
        
        % A3
        Iukronxvxinv = sparse(D3UsedColumns, D3UsedSpecies1, xvxinv(D3UsedSpecies2), nx*nu,nu);
        ukronxdvxinvdu = sparse(nu*nx,nu);
        ukronxdvxinvdu(D3UsedColumns,:) = bsxfun(@times, u(D3UsedSpecies1) .* x(D3UsedSpecies2), dvxinvdu(D3UsedSpecies2,:));
        
        % A4
        xkronIuvuinv  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vuinv(D4UsedSpecies2), nu*nx,nu);
        xkronudvuinvdu = sparse(nx*nu,nu);
        xkronudvuinvdu(D4UsedColumns,:) = bsxfun(@times, x(D4UsedSpecies1) .* u(D4UsedSpecies2), dvuinvdu(D4UsedSpecies2,:));
        
        % A5
        Iukronuvuinv = sparse(D5UsedColumns, D5UsedSpecies1, uvuinv(D5UsedSpecies2), nu*nu,nu);
        ukronIuvuinv  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vuinv(D5UsedSpecies2), nu*nu,nu);
        ukronudvuinvdu = sparse(nu*nu,nu);
        ukronudvuinvdu(D5UsedColumns,:) = bsxfun(@times, u(D5UsedSpecies1) .* u(D5UsedSpecies2), dvuinvdu(D5UsedSpecies2,:));
        
        val = dD2dk_rk_xx * xkronxdvxinvdu ...
              + dD3dk_rk_ux * (Iukronxvxinv + ukronxdvxinvdu) ...
              + dD4dk_rk_xu * (xkronIuvuinv + xkronudvuinvdu) ...
              + dD5dk_rk_uu * (Iukronuvuinv + ukronIuvuinv + ukronudvuinvdu) ...
              + dD6dk_rk_u; % rk_u
    end
end

function array = make_room(array, needed)
% Make room in an array by doubling its rows if needed is larger than its
% current length. Fill the empty spots with zeros.
current = size(array, 1);
if needed > current
    addlength = max(current, needed-current);
    array = [array; zeros(addlength, size(array,2))];
end
end

function S = sparse(varargin)
% Redefined sparse() function that constructs full matrices using
% the sparse syntax if input arguments are symbolic.
S = sym_sparse(varargin{:});
end

function c = bsxfun(varargin)
% Redefined bsxfun() that works with symbolic arrays.
c = sym_bsxfun(varargin{:});
end
