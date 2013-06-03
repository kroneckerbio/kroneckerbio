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

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
assert(nargin >= 1, 'KroneckerBio:FinalizeModel:TooFewInputs', 'FinalizeModel requires at least 1 input argument')
assert(isscalar(m), 'KroneckerBio:FinalizeModel:MoreThanOneModel', 'The model structure must be scalar')

%% Place compartments
nvNew = m.add.nv;

% Check if item by this name already exists
handled = false(nvNew,1);
for iv = 1:nvNew
    % Check existing compartments
    matchPosition = find(strcmp(m.add.Compartments(iv).Name, {m.Compartments.Name}));
    if ~isempty(matchPosition)
        handled(iv) = true;
        m.Compartments(matchPosition) = m.add.Compartments(iv);
    end
    
    % Check newly added compartments
    matchPosition = find(strcmp(m.add.Compartments(iv).Name, {m.add.Compartments(1:iv-1).Name}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
    end
end

% Append new item
m.Compartments = [m.Compartments; m.add.Compartments(~handled)];

% Update count
m.nv = numel(m.Compartments);

%% Place species
nxuNew = m.add.nxu;

% Check if item by this name already exists
handled = false(nxuNew,1);
xuNewUsed = true(nxuNew,1);
for ixu = 1:nxuNew
    % Check existing states
    matchPosition = find(strcmp(m.add.Species(ixu).Name, {m.Species.Name}) & strcmp(m.add.Species(ixu).Compartment, {m.Species.Compartment}));
    if ~isempty(matchPosition)
        handled(ixu) = true;
        m.Species(matchPosition) = rmfield(m.add.Species(ixu), 'Units');
    end
    
    % Check newly added states
    matchPosition = find(strcmp(m.add.Species(ixu).Name, {m.add.Species(1:ixu-1).Name}) & strcmp(m.add.Species(ixu).Compartment, {m.add.Species(1:ixu-1).Compartment}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
        xuNewUsed(matchPosition) = false;
    end
end

% Append new item
m.Species = [m.Species; rmfield(m.add.Species(~handled), 'Units')];

% Update count
isu = cat(1, m.Species.IsInput, false(0,1));
m.isu = isu;
m.xInd = find(~isu);
m.uInd = find(isu);
m.nu = nnz(isu);
m.nx = numel(m.Species) - m.nu;
xuNamesFull = strcat({m.Species.Compartment}, '.', {m.Species.Name});
xNamesFull = strcat({m.Species(~isu).Compartment}, '.', {m.Species(~isu).Name});
uNamesFull = strcat({m.Species(isu).Compartment}, '.', {m.Species(isu).Name});

%% Place outputs
nyNew = m.add.ny;

% Check if item by this name already exists
handled = false(nyNew,1);
for iy = 1:nyNew
    % Check existing outputs
    matchPosition = find(strcmp(m.add.Outputs(iy).Name, {m.Outputs.Name}));
    if ~isempty(matchPosition)
        handled(iy) = true;
        m.Outputs(matchPosition) = m.add.Outputs(iy);
    end
    
    % Check newly added outputs
    matchPosition = find(strcmp(m.add.Outputs(iy).Name, {m.add.Outputs(1:iy-1).Name}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
    end
end

% Append new items
m.Outputs = [m.Outputs; m.add.Outputs(~handled)];

% Update count
m.ny = numel(m.Outputs);

%% Place parameters
nkNew = m.add.nk;

% Check if item by this name already exists
handled = false(nkNew,1);
for ik = 1:nkNew
    % Check existing parameters
    matchPosition = find(strcmp(m.add.Parameters(ik).Name, {m.Parameters.Name}));
    if ~isempty(matchPosition)
        handled(ik) = true;
        m.Parameters(matchPosition) = rmfield(m.add.Parameters(ik), 'Units');
    end
    
    % Check newly added parameters
    matchPosition = find(strcmp(m.add.Parameters(ik).Name, {m.add.Parameters(1:ik-1).Name}));
    if ~isempty(matchPosition)
        handled(matchPosition) = true;
    end
end

% Append new items
m.Parameters = [m.Parameters; rmfield(m.add.Parameters(~handled), 'Units')];

% Update count
m.nk = numel(m.Parameters);

%% Expand reactions
nrNew = m.add.nr; % Number of reaction specifications

% Loop over added reaction specifications
rNew = cell(nrNew,1); % Containers for all the spawned reactions
rCount = zeros(nrNew,1); % Number of reactions that have spawned
for ir = 1:nrNew
    % Determine possible compartments as a cell array of strings
    if isempty(m.add.Reactions(ir).Compartment)
        possibleComp = {m.Compartments.Name};
    elseif ischar(m.add.Reactions(ir).Compartment)
        possibleComp = {m.add.Reactions(ir).Compartment};
    elseif iscell(m.add.Reactions(ir).Compartment)
        possibleComp = m.add.Reactions(ir).Compartment;
    end
    
    % Container for all the reactions that come from this specification
    nvPossible = numel(possibleComp);
    rNew{ir} = emptyReactions(nvPossible);
    
    % Match reactants and products
    if any(m.add.Reactions(ir).Reactants{1} == '.')
        reactant1Complete = true;
        possibleReactant1 = strcmp(m.add.Reactions(ir).Reactants{1}, xuNamesFull);
    else
        reactant1Complete = false;
        possibleReactant1 = strcmp(m.add.Reactions(ir).Reactants{1}, {m.Species.Name});
    end
    if any(m.add.Reactions(ir).Reactants{2} == '.')
        reactant2Complete = true;
        possibleReactant2 = strcmp(m.add.Reactions(ir).Reactants{2}, xuNamesFull);
    else
        reactant2Complete = false;
        possibleReactant2 = strcmp(m.add.Reactions(ir).Reactants{2}, {m.Species.Name});
    end
    if any(m.add.Reactions(ir).Products{1} == '.')
        product1Complete = true;
        possibleProduct1 = strcmp(m.add.Reactions(ir).Products{1}, xuNamesFull);
    else
        product1Complete = false;
        possibleProduct1 = strcmp(m.add.Reactions(ir).Products{1}, {m.Species.Name});
    end
    if any(m.add.Reactions(ir).Products{2} == '.')
        product2Complete = true;
        possibleProduct2 = strcmp(m.add.Reactions(ir).Products{2}, xuNamesFull);
    else
        product2Complete = false;
        possibleProduct2 = strcmp(m.add.Reactions(ir).Products{2}, {m.Species.Name});
    end
    
    % Loop over possible compartments
    aReactionWasFound = false; % Used to throw a warning if this remains false
    for iv = 1:nvPossible
        % Match compartments
        possibleCompartment = strcmp(possibleComp{iv}, {m.Species.Compartment});
        
        % Find reactants and products for this compartment
        reactant1Exists = ~isempty(m.add.Reactions(ir).Reactants{1});
        reactant2Exists = ~isempty(m.add.Reactions(ir).Reactants{2});
        product1Exists  = ~isempty(m.add.Reactions(ir).Products{1});
        product2Exists  = ~isempty(m.add.Reactions(ir).Products{2});
        reactant1Found  = ~reactant1Exists || (any(possibleReactant1) && reactant1Complete) || any(possibleReactant1 & possibleCompartment);
        reactant2Found  = ~reactant2Exists || (any(possibleReactant2) && reactant2Complete) || any(possibleReactant2 & possibleCompartment);
        product1Found   = ~product1Exists  || (any(possibleProduct1)  && product1Complete)  || any(possibleProduct1 & possibleCompartment);
        product2Found   = ~product2Exists  || (any(possibleProduct2)  && product2Complete)  || any(possibleProduct2 & possibleCompartment);
        parameter1Exists = ~isempty(m.add.Reactions(ir).Parameters{1});
        parameter2Exists = ~isempty(m.add.Reactions(ir).Parameters{2});
        
        % Check for errors
        if ((reactant1Exists || reactant2Exists) && (reactant1Found && reactant2Found) && (~product1Found || ~product2Found) && parameter1Exists)
            error('KroneckerBio:FinalizeModel:MissingProduct', 'Reaction %s (#%i) has the necessary reactants in compartment %s but not the necessary products', m.add.Reactions(ir).Names{1}, ir, possibleComp{iv})
        end
        if ((product1Exists || product2Exists) && (product1Found && product2Found) && (~reactant1Found || ~reactant2Found) && parameter2Exists)
            error('KroneckerBio:FinalizeModel:MissingProduct', 'Reverse reaction %s (#%i) has the necessary reactants in compartment %s but not the necessary products', m.add.Reactions(ir).Names{2}, ir, possibleComp{iv})
        end
        
        % Check if this reaction takes place
        if ~(reactant1Found && reactant2Found && product1Found && product2Found)
            continue
        end
        
        % Forward reaction
        if parameter1Exists
            % Forward reaction exists in this compartment
            rCount(ir) = rCount(ir) + 1;
            
            % Add reaction to list
            rNew{ir}(rCount(ir)).Name         = m.add.Reactions(ir).Names{1};
            rNew{ir}(rCount(ir)).Reactants{1} = fixReactionSpecies(m.add.Reactions(ir).Reactants{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Reactants{2} = fixReactionSpecies(m.add.Reactions(ir).Reactants{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{1}  = fixReactionSpecies(m.add.Reactions(ir).Products{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{2}  = fixReactionSpecies(m.add.Reactions(ir).Products{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Parameter    = m.add.Reactions(ir).Parameters{1};
            
            aReactionWasFound = true;
        end
        
        % Reverse reaction
        if parameter2Exists
            % Forward reaction exists in this compartment
            rCount(ir) = rCount(ir) + 1;
            
            % Add reaction to list
            rNew{ir}(rCount(ir)).Name         = m.add.Reactions(ir).Names{2};
            rNew{ir}(rCount(ir)).Reactants{1} = fixReactionSpecies(m.add.Reactions(ir).Products{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Reactants{2} = fixReactionSpecies(m.add.Reactions(ir).Products{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{1}  = fixReactionSpecies(m.add.Reactions(ir).Reactants{1}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Products{2}  = fixReactionSpecies(m.add.Reactions(ir).Reactants{2}, possibleComp{iv});
            rNew{ir}(rCount(ir)).Parameter    = m.add.Reactions(ir).Parameters{2};

            aReactionWasFound = true;
        end
        
        % Exit if all compartments are specified
        if (reactant1Complete || ~reactant1Exists) && (reactant2Complete || ~reactant2Exists) && (product1Complete || ~product1Exists) && (product2Complete || ~product2Exists)
            break
        end
    end
    
    % Warn if no compartment had the species for this reaction specification
    if ~aReactionWasFound
        warning('KroneckerBio:FinalizeModel:UnusedReaction', 'Reaction %s (#%i) was not found to take place in any compartment', m.add.Reactions(ir).Names{1}, ir)
    end
end

%% Place reactions
% Total number of new elementary reactions
nrTotal = sum(rCount);

% Add room for new reactions
m.Reactions = [m.Reactions; emptyReactions(nrTotal)];

% Insert items
rEndIndex = m.nr;
for ir = 1:nrNew
    rStartIndex = rEndIndex + 1;
    rEndIndex   = rEndIndex + rCount(ir);
    m.Reactions(rStartIndex:rEndIndex) = rNew{ir}(1:rCount(ir));
end

% Update count
m.nr = numel(m.Reactions);

%% Useful information
% Constants
nv = m.nv;
nx = m.nx;
nu = m.nu;
nxu = nx + nu;
ny = m.ny;
nk = m.nk;
nr = m.nr;

% State compartments
m.vxInd = zeros(nx,1);
xInd = find(~isu);
for ix = 1:nx
    m.vxInd(ix) = find(strcmp(m.Species(xInd(ix)).Compartment, {m.Compartments.Name}));
end

% Input compartments
m.vuInd = zeros(nu,1);
uInd = find(isu);
for iu = 1:nu
    m.vuInd(iu) = find(strcmp(m.Species(uInd(iu)).Compartment, {m.Compartments.Name}));
end

% Map from species to x0 and u
x0uInd = zeros(nxu,1);
x0uInd(~isu) = 1:nx;
x0uInd(isu) = 1:nu;

%% Process compartments
% Dimensions
m.dv = [m.Compartments.Dimension].';

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
    nExpr = numel(m.Compartments(iv).Expressions);
    for iExpr = 1:nExpr
        % Find states that match the expression
        match = find(~cellfun(@isempty, regexp(xNamesFull, m.Compartments(iv).Expressions{iExpr}, 'once')));
        nAdd = numel(match);
        nB1Entries = nB1Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(B1Entries,1);
        if nB1Entries > currentLength
            addlength = max(currentLength, nAdd);
            B1Entries = [B1Entries; zeros(addlength,2)];
        end
        
        % Add entries
        B1Entries(nB1Entries-nAdd+1:nB1Entries,1) = iv;
        B1Entries(nB1Entries-nAdd+1:nB1Entries,2) = match;
        B1Values(nB1Entries-nAdd+1:nB1Entries) = m.Compartments(iv).Values(iExpr);

        % Find inputs that match the expression
        match = find(~cellfun(@isempty, regexp(uNamesFull, m.Compartments(iv).Expressions{iExpr}, 'once')));
        nAdd = numel(match);
        nB2Entries = nB2Entries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(B2Entries,1);
        if nB2Entries > currentLength
            addlength = max(currentLength, nAdd);
            B2Entries = [B2Entries; zeros(addlength,2)];
        end
        
        % Add entries
        B2Entries(nB2Entries-nAdd+1:nB2Entries,1) = iv;
        B2Entries(nB2Entries-nAdd+1:nB2Entries,2) = match;
        B2Values(nB2Entries-nAdd+1:nB2Entries) = m.Compartments(iv).Values(iExpr);

        % Find empty expressions, which are constants
        if isempty(m.Compartments(iv).Expressions{iExpr})
            nbEntries = nbEntries + 1;
            
            % Add more room in vector if necessary
            currentLength = size(bEntries,1);
            if nbEntries > currentLength
                addlength = max(currentLength, 1);
                bEntries = [bEntries; zeros(addlength,2)];
            end
            
            % Add entries
            bEntries(nbEntries,1) = iv;
            bEntries(nbEntries,2) = 1;
            bValues(nbEntries) = m.Compartments(iv).Values(iExpr);
        end
    end
end

% Remove duplicate entries
[B1Entries ind] = unique(B1Entries(1:nB1Entries,:), 'rows');
B1Values = B1Values(ind);

[B2Entries ind] = unique(B2Entries(1:nB2Entries,:), 'rows');
B2Values = B2Values(ind);

[bEntries ind] = unique(bEntries(1:nbEntries,:), 'rows');
bValues = bValues(ind);

% Construct matrices
m.B1 = sparse(B1Entries(:,1), B1Entries(:,2), B1Values, m.nv, m.nx);
m.B2 = sparse(B2Entries(:,1), B2Entries(:,2), B2Values, m.nv, m.nu);
m.b  = sparse(bEntries(:,1),  bEntries(:,2),  bValues,  m.nv, 1);

%% Compute non-concentration species values
x0Handled = true(nx,1);
uHandled = true(nu,1);
xuNewInd = zeros(nxuNew,1);
for ixuNew = find(xuNewUsed)'
    % Index in Species.Names that this new Species applies
    xuIndi = find(strcmp(m.add.Species(ixuNew).Name, {m.Species.Name}) & strcmp(m.add.Species(ixuNew).Compartment, {m.Species.Compartment}));
    xuNewInd(ixuNew) = xuIndi;
    
    if m.Species(xuIndi).IsInput
        [m.Species(xuIndi).Value.Function, uHandled(x0uInd(xuIndi))] = fixInputValueFunction(m.Species(xuIndi).Value.Function, m.add.Species(ixuNew).Units);
    else
        [m.Species(xuIndi).Value, x0Handled(x0uInd(xuIndi))] = fixStateInitialValue(m.Species(xuIndi).Value, m.add.Species(ixuNew).Units);
    end
end

%% Compute starting compartment volumes
% Assert that no unhandled species play a part in the compartment volume
if nnz(m.B1(:,~x0Handled)) > 0
    errorInd = find(m.B1(:,~x0Handled));
    errorInd = ind2sub(errorInd, m.nv, nnz(~x0Handled));
    error('KroneckerBio:FinalizeModel:ConcentrationAffectsVolume', 'Initial condition for species %s was given in concentration, but its value affects compartment %s; the amount in the concentration cannot be resolved', m.States(errorInd(2)).Name, m.Compartments(errorInd(1)).Name)
end
if nnz(m.B2(:,~uHandled)) > 0
    errorInd = find(m.B1(:,~x0Handled));
    errorInd = ind2sub(errorInd, m.nv, nnz(~x0Handled));
    error('KroneckerBio:FinalizeModel:ConcentrationAffectsVolume', 'Input %s was given in concentration, but its value affects a compartment; the amount in the concentration cannot be resolved', m.Inputs(errorInd(2)).Name, m.Compartments(errorInd(1)).Name)
end

x0 = cat(1, m.Species(~isu).Value, zeros(0,1));
u = completeInputFunction(m.Species(isu));
v0 = m.B1 * x0 + m.B2 * u(0) + m.b;

%% Compute remaining species values
for ixuNew = find(xuNewUsed)'
    % Index in m.Species matching this new Species
    xuIndi = xuNewInd(ixuNew);
    
    % Only do something if it was not already handled
    if m.Species(xuIndi).IsInput && ~uHandled(x0uInd(xuIndi))
        % Index in m.Compartments
        vInd = strcmp(m.Species(xuIndi).Compartment, {m.Compartments.Name});
        
        m.Species(xuIndi).Value.Function = fixInputValueFunction(m.Species(xuIndi).Value.Function, m.add.Species(ixuNew).Units, v0(vInd));
    elseif ~m.Species(xuIndi).IsInput && ~x0Handled(x0uInd(xuIndi))
        % Index in m.Compartments
        vInd = strcmp(m.Species(xuIndi).Compartment, {m.Compartments.Name});
        
        m.Species(xuIndi).Value = fixStateInitialValue(m.Species(xuIndi).Value, m.add.Species(ixuNew).Units, v0(vInd));
    end
end

%% Process states
% Put initial conditions into x0 vector
m.x0 = cat(1, m.Species(~isu).Value, zeros(0,1));

%% Process inputs
% Condense time varying inputs into u(t) function
m.u = completeInputFunction(m.Species(isu));

% Put input parameters into q vector
m.nqu = zeros(nu,1);
uInd = find(isu);
for iu = 1:nu
    m.nqu(iu) = numel(m.Species(uInd(iu)).Value.Parameters);
end
uValue = cat(1, m.Species(isu).Value, struct('Function', cell(0,1), 'Parameters', cell(0,1)));
m.q = cat(1, uValue.Parameters, zeros(0,1));
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
    nExpr = numel(m.Outputs(iy).Expressions);
    for iExpr = 1:nExpr
        % Find states that match the expression
        match = find(~cellfun(@isempty, regexp(xNamesFull, m.Outputs(iy).Expressions{iExpr}, 'once')));
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
        C1Values(nC1Entries-nAdd+1:nC1Entries) = m.Outputs(iy).Values(iExpr);
        
        % Find inputs that match the expression
        match = find(~cellfun(@isempty, regexp(uNamesFull, m.Outputs(iy).Expressions{iExpr}, 'once')));
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
        C2Values(nC2Entries-nAdd+1:nC2Entries) = m.Outputs(iy).Values(iExpr);

        % Find empty expressions, which are constants
        if isempty(m.Outputs(iy).Expressions{iExpr})
            ncEntries = ncEntries + 1;
            
            % Add more room in vector if necessary
            currentLength = size(cEntries,1);
            if nbEntries > currentLength
                addlength = max(currentLength, 1);
                cEntries = [cEntries; zeros(addlength,2)];
                cValues = [cValues; zeros(addlength,1)];
            end
            
            % Add entries
            cEntries(ncEntries,1) = iy;
            cEntries(ncEntries,2) = 1;
            cValues(ncEntries) = m.Outputs(iy).Values(iExpr);
        end
    end
end

% Remove duplicate entries
[C1Entries ind] = unique(C1Entries(1:nC1Entries,:), 'rows');
C1Values = C1Values(ind);

[C2Entries ind] = unique(C2Entries(1:nC2Entries,:), 'rows');
C2Values = C2Values(ind);

[cEntries ind] = unique(cEntries(1:ncEntries,:), 'rows');
cValues = cValues(ind);

% Construct matrices
m.C1 = sparse(C1Entries(:,1), C1Entries(:,2), C1Values, m.ny, m.nx);
m.C2 = sparse(C2Entries(:,1), C2Entries(:,2), C2Values, m.ny, m.nu);
m.c  = sparse(cEntries(:,1),  cEntries(:,2),  cValues,  m.ny, 1);

%% Process parameters
% Put rate parameter into k vector
m.k = cat(1, m.Parameters.Value, zeros(0,1));

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
    if isempty(m.Reactions(ir).Reactants{1})
        reactant1Exists = false;
    else
        reactant1Exists = true;
        reactant1IsState = true;
        reactant1 = find(strcmp(m.Reactions(ir).Reactants{1}, xNamesFull));
        if isempty(reactant1)
            reactant1IsState = false;
            reactant1 = find(strcmp(m.Reactions(ir).Reactants{1}, uNamesFull));
        end
    end
    
    if isempty(m.Reactions(ir).Reactants{2})
        reactant2Exists = false;
    else
        reactant2Exists = true;
        reactant2IsState = true;
        reactant2 = find(strcmp(m.Reactions(ir).Reactants{2}, xNamesFull));
        if isempty(reactant2)
            reactant2IsState = false;
            reactant2 = find(strcmp(m.Reactions(ir).Reactants{2}, uNamesFull));
        end
    end
    
    m.rOrder(ir) = reactant1Exists + reactant2Exists; % Store order
    
    product1 = find(strcmp(m.Reactions(ir).Products{1}, xNamesFull));
    product1IsState = ~isempty(product1);
    
    product2 = find(strcmp(m.Reactions(ir).Products{2}, xNamesFull));
    product2IsState = ~isempty(product2);
    
    parameter = find(strcmp(m.Reactions(ir).Parameter, {m.Parameters.Name}));
    assert(~isempty(parameter), 'KroneckerBio:FinalizeModel:MissingReactionParameter', 'Reaction %s (#%i) requires parameter %s, but no parameter by that name was found', m.Reactions(ir).Name, ir, m.Reactions(ir).Parameter)
    m.krInd(ir) = parameter; % Store index
    
    % Switch on reactant state
    if reactant1Exists && ~reactant2Exists && reactant1IsState
        % D1/A1 reaction
        % Add S entries
        nAdd = 1 + product1IsState + product2IsState;
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Subtract reactant
        SEntries(nSEntries-nAdd+1,1) = reactant1;
        SEntries(nSEntries-nAdd+1,2) = ir;
        SValues(nSEntries-nAdd+1)    = -1;
        
        % Add product 1
        if product1IsState
            SEntries(nSEntries-nAdd+2,1) = product1;
            SEntries(nSEntries-nAdd+2,2) = ir;
            SValues(nSEntries-nAdd+2)    = 1;
        end
        
        % Add product 2
        if product2IsState
            SEntries(nSEntries-nAdd+product1IsState+2,1) = product2;
            SEntries(nSEntries-nAdd+product1IsState+2,2) = ir;
            SValues(nSEntries-nAdd+product1IsState+2)    = 1;
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
        
        dD1dkEntries(nD1Entries,1) = sub2ind([nr,nx], ir, reactant1);
        dD1dkEntries(nD1Entries,2) = parameter;
        dD1dkValues(nD1Entries)    = 1;

    elseif reactant1Exists && reactant2Exists && reactant1IsState && reactant2IsState
        % D2/A2 reaction
        % Order reactants so that freest species is second
        if m.dv(m.vxInd(reactant1)) > m.dv(m.vxInd(reactant2)) || (m.dv(m.vxInd(reactant1)) == m.dv(m.vxInd(reactant2)) && reactant1 > reactant2)
            [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
        end
        
        % Add S entries
        nAdd = 2 + product1IsState + product2IsState;
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Subtract reactant 1
        SEntries(nSEntries-nAdd+1,1) = reactant1;
        SEntries(nSEntries-nAdd+1,2) = ir;
        SValues(nSEntries-nAdd+1)    = -1;
        
        % Subtract reactant 2
        SEntries(nSEntries-nAdd+2,1) = reactant2;
        SEntries(nSEntries-nAdd+2,2) = ir;
        SValues(nSEntries-nAdd+2)    = -1;

        % Add product 1
        if product1IsState
            SEntries(nSEntries-nAdd+3,1) = product1;
            SEntries(nSEntries-nAdd+3,2) = ir;
            SValues(nSEntries-nAdd+3)    = 1;
        end
        
        % Add product 2
        if product2IsState
            SEntries(nSEntries-nAdd+product1IsState+3,1) = product2;
            SEntries(nSEntries-nAdd+product1IsState+3,2) = ir;
            SValues(nSEntries-nAdd+product1IsState+3)    = 1;
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
        
        dD2dkEntries(nD2Entries,1) = sub2ind([nr,nx,nx], ir, reactant2, reactant1);
        dD2dkEntries(nD2Entries,2) = parameter;
        dD2dkValues(nD2Entries)    = 1;
        
    elseif reactant1Exists && reactant2Exists && xor(reactant2IsState, reactant1IsState)
        % Order reactants so that freest species is second
        if reactant1IsState
            if m.dv(m.vxInd(reactant1)) > m.dv(m.vuInd(reactant2))
                [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
                [reactant1IsState reactant2IsState] = deal(reactant2IsState, reactant1IsState);
            end
        else
            if m.dv(m.vuInd(reactant1)) >= m.dv(m.vxInd(reactant2))
                [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
                [reactant1IsState reactant2IsState] = deal(reactant2IsState, reactant1IsState);
            end
        end
        
        if reactant2IsState
            % D3/A3 reaction
            % Add S entries
            nAdd = 1 + product1IsState + product2IsState;
            nSEntries = nSEntries + nAdd;
            
            % Add more room in vector if necessary
            currentLength = size(nSEntries,1);
            if nSEntries > currentLength
                addlength = max(currentLength, nAdd);
                SEntries = [SEntries; zeros(addlength,2)];
                SValues  = [SValues;  zeros(addlength,1)];
            end
            
            % Subtract reactant 2
            SEntries(nSEntries-nAdd+1,1) = reactant2;
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = -1;
            
            % Add product 1
            if product1IsState
                SEntries(nSEntries-nAdd+2,1) = product1;
                SEntries(nSEntries-nAdd+2,2) = ir;
                SValues(nSEntries-nAdd+2)    = 1;
            end
            
            % Add product 2
            if product2IsState
                SEntries(nSEntries-nAdd+product1IsState+2,1) = product2;
                SEntries(nSEntries-nAdd+product1IsState+2,2) = ir;
                SValues(nSEntries-nAdd+product1IsState+2)    = 1;
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
            
            dD3dkEntries(nD3Entries,1) = sub2ind([nr,nx,nu], ir, reactant2, reactant1);
            dD3dkEntries(nD3Entries,2) = parameter;
            dD3dkValues(nD3Entries)    = 1;
            
        else%reactant1IsState
            % D4/A4 reaction
            % Add S entries
            nAdd = 1 + product1IsState + product2IsState;
            nSEntries = nSEntries + nAdd;
            
            % Add more room in vector if necessary
            currentLength = size(nSEntries,1);
            if nSEntries > currentLength
                addlength = max(currentLength, nAdd);
                SEntries = [SEntries; zeros(addlength,2)];
                SValues  = [SValues;  zeros(addlength,1)];
            end
            
            % Subtract reactant 1
            SEntries(nSEntries-nAdd+1,1) = reactant1;
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = -1;
            
            % Add product 1
            if product1IsState
                SEntries(nSEntries-nAdd+2,1) = product1;
                SEntries(nSEntries-nAdd+2,2) = ir;
                SValues(nSEntries-nAdd+2)    = 1;
            end
            
            % Add product 2
            if product2IsState
                SEntries(nSEntries-nAdd+product1IsState+2,1) = product2;
                SEntries(nSEntries-nAdd+product1IsState+2,2) = ir;
                SValues(nSEntries-nAdd+product1IsState+2)    = 1;
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
        
            dD4dkEntries(nD4Entries,1) = sub2ind([nr,nu,nx], ir, reactant2, reactant1);
            dD4dkEntries(nD4Entries,2) = parameter;
            dD4dkValues(nD4Entries)    = 1;
        
        end
        
    elseif reactant1Exists && reactant2Exists && ~reactant1IsState && ~reactant2IsState
        % D5/A5 reaction
        % Order reactants so that freest species is second
        if m.dv(m.vuInd(reactant1)) > m.dv(m.vuInd(reactant2)) || (m.dv(m.vuInd(reactant1)) == m.dv(m.vuInd(reactant2)) && reactant1 > reactant2)
            [reactant1 reactant2] = deal(reactant2, reactant1); % Swap
        end
        
        % Add S entries
        nAdd = product1IsState + product2IsState;
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Add product 1
        if product1IsState
            SEntries(nSEntries-nAdd+1,1) = product1;
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = 1;
        end
        
        % Add product 2
        if product2IsState
            SEntries(nSEntries-nAdd+product1IsState+1,1) = product2;
            SEntries(nSEntries-nAdd+product1IsState+1,2) = ir;
            SValues(nSEntries-nAdd+product1IsState+1)    = 1;
        end
        
        % Add D5 entry
        nD5Entries = nD5Entries + 1;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, 1);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        dD5dkEntries(nD5Entries,1) = sub2ind([nr,nu,nu], ir, reactant2, reactant1);
        dD5dkEntries(nD5Entries,2) = parameter;
        dD5dkValues(nD5Entries)    = 1;
        
    elseif reactant1Exists && ~reactant2Exists && ~reactant1IsState
        % D6/A6 reaction
        % Add S entries
        nAdd = product1IsState + product2IsState;
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Add product 1
        if product1IsState
            SEntries(nSEntries-nAdd+1,1) = product1;
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = 1;
        end
        
        % Add product 2
        if product2IsState
            SEntries(nSEntries-nAdd+product1IsState+1,1) = product2;
            SEntries(nSEntries-nAdd+product1IsState+1,2) = ir;
            SValues(nSEntries-nAdd+product1IsState+1)    = 1;
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
        
        dD6dkEntries(nD6Entries,1) = sub2ind([nr,nu], ir, reactant1);
        dD6dkEntries(nD6Entries,2) = parameter;
        dD6dkValues(nD6Entries)    = 1;
        
    elseif ~reactant1Exists && ~reactant2Exists
        % d/a reaction
        % Add S entries
        nAdd = product1IsState + product2IsState;
        nSEntries = nSEntries + nAdd;
        
        % Add more room in vector if necessary
        currentLength = size(nSEntries,1);
        if nSEntries > currentLength
            addlength = max(currentLength, nAdd);
            SEntries = [SEntries; zeros(addlength,2)];
            SValues  = [SValues;  zeros(addlength,1)];
        end
        
        % Add product 1
        if product1IsState
            SEntries(nSEntries-nAdd+1,1) = product1;
            SEntries(nSEntries-nAdd+1,2) = ir;
            SValues(nSEntries-nAdd+1)    = 1;
        end
        
        % Add product 2
        if product2IsState
            SEntries(nSEntries-nAdd+product1IsState+1,1) = product2;
            SEntries(nSEntries-nAdd+product1IsState+1,2) = ir;
            SValues(nSEntries-nAdd+product1IsState+1)    = 1;
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
        dddkValues(ndEntries)    = 1;
        
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
%m.A2 = reshape(mtimestall(m.dA2dk, kRand), nx,nx*nx); % Matlab fail
m.A3 = reshape(m.dA3dk * kRand, nx,nu*nx);
m.A4 = reshape(m.dA4dk * kRand, nx,nx*nu);
m.A5 = reshape(m.dA5dk * kRand, nx,nu*nu);
m.A6 = reshape(m.dA6dk * kRand, nx,nu);
m.a  = m.dadk * kRand;

m.D1 = reshape(m.dD1dk * kRand, nr,nx);
m.D2 = reshape(m.dD2dk * kRand, nr,nx*nx);
%m.D2 = reshape(mtimestall(m.dD2dk, kRand), nr,nx*nx);
m.D3 = reshape(m.dD3dk * kRand, nr,nu*nx);
m.D4 = reshape(m.dD4dk * kRand, nr,nx*nu);
m.D5 = reshape(m.dD5dk * kRand, nr,nu*nu);
m.D6 = reshape(m.dD6dk * kRand, nr,nu);
m.d  = m.dddk * kRand;

% Determine used columns of bimolecular matrices
[unused D2UsedColumns] = find(m.D2);
[unused D3UsedColumns] = find(m.D3);
[unused D4UsedColumns] = find(m.D4);
[unused D5UsedColumns] = find(m.D5);

D2UsedColumns = unique(D2UsedColumns);
D3UsedColumns = unique(D3UsedColumns);
D4UsedColumns = unique(D4UsedColumns);
D5UsedColumns = unique(D5UsedColumns);

[D2UsedSpecies2 D2UsedSpecies1] = ind2sub([nx,nx], D2UsedColumns);
[D3UsedSpecies2 D3UsedSpecies1] = ind2sub([nx,nu], D3UsedColumns);
[D4UsedSpecies2 D4UsedSpecies1] = ind2sub([nu,nx], D4UsedColumns);
[D5UsedSpecies2 D5UsedSpecies1] = ind2sub([nu,nu], D5UsedColumns);

%% Empty out added items
m.add.nv  = 0;
m.add.nxu = 0;
m.add.ny  = 0;
m.add.nk  = 0;
m.add.nr  = 0;
m.add.Compartments = growCompartments([],0);
m.add.Species      = growSpecies([],0);
m.add.Outputs      = growOutputs([],0);
m.add.Parameters   = growParameters([],0);
m.add.Reactions    = growReactions([],0);

%% Final build of model
m = final(m, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2);

end

function m = final(m, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2)

% Constants
nx = m.nx;
nu = m.nu;
nk = m.nk;
nr = m.nr;
isu = cat(1, m.Species.IsInput, false(0,1));

% Build kronecker matrices
sparsek = sparse(m.k);
m.A1 = reshape(m.dA1dk * sparsek, nx,nx);
%m.A2 = reshape(m.dA2dk * sparsek, nx,nx*nx);
m.A2 = reshape(mtimestall(m.dA2dk, sparsek), nx,nx*nx); % Matlab fail
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
m.f = fHidden(m.A1, m.A2, m.A3, m.A4, m.A5, m.A6, m.a, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.dfdx = dfdxHidden(m.A1, m.A2, m.A3, m.A4, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.dfdu = dfduHidden(m.A3, m.A4, m.A5, m.A6, m.B1, m.B2, m.b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.dfdk = dfdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.dadk, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.d2fdx2  = d2fdx2Hidden(m.A2, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies2, m.vxInd);
m.d2fdu2  = d2fdu2Hidden(m.A5, m.B1, m.B2, m.b, D5UsedColumns, D5UsedSpecies2, m.vuInd);
m.d2fdk2  = d2fdk2Hidden(m.nx, m.nk);
m.d2fdudx = d2fdudxHidden(m.A3, m.A4, m.B1, m.B2, m.b, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdxdu = d2fdxduHidden(m.A3, m.A4, m.B1, m.B2, m.b, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdkdx = d2fdkdxHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdxdk = d2fdxdkHidden(m.dA1dk_fk_x, m.dA2dk_fk_xx, m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdkdu = d2fdkduHidden(m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.B1, m.B2, m.b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.d2fdudk = d2fdudkHidden(m.dA3dk_fk_ux, m.dA4dk_fk_xu, m.dA5dk_fk_uu, m.dA6dk_fk_u, m.B1, m.B2, m.b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.r = rHidden(m.D1, m.D2, m.D3, m.D4, m.D5, m.D6, m.d, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.drdx = drdxHidden(m.D1, m.D2, m.D3, m.D4, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.drdu = drduHidden(m.D3, m.D4, m.D5, m.D6, m.B1, m.B2, m.b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);
m.drdk = drdkHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.dD5dk_rk_uu, m.dD6dk_rk_u, m.dddk, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, m.vxInd, m.vuInd);

m.d2rdx2  = d2rdx2Hidden(m.D2, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies2, m.vxInd);
m.d2rdk2  = d2rdk2Hidden(m.nr, m.nk);
m.d2rdxdk = d2rdxdkHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd);
m.d2rdkdx = d2rdkdxHidden(m.dD1dk_rk_x, m.dD2dk_rk_xx, m.dD3dk_rk_ux, m.dD4dk_rk_xu, m.B1, m.B2, m.b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, m.vxInd, m.vuInd, m.nr);

m.Ready = true;
m.Update = @Update;

    function mout = Update(k, x0, q)
        % Copy existing model
        mout = m;
        
        % Apply changes
        mout.k = k;
        mout.x0 = x0;
        mout.q = q;
        
        % Distribute values
        if m.nk >= 1
            k = num2cell(k);
            [mout.Parameters.Value] = k{:};
        end
        
        if m.nx >= 1
            x0 = num2cell(x0);
            [mout.Species(~isu).Value] = x0{:};
        end
        
        qIndex = 0;
        uInd = find(isu);
        for iu = 1:nu
            mout.Species(uInd(iu)).Value.Parameters = q(qIndex+1:qIndex+m.nqu(iu));
            qIndex = qIndex + m.nqu(iu);
        end
        
        % Rebuild model
        mout = final(mout, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2);
    end
end

function handle = fHidden(A1, A2, A3, A4, A5, A6, a, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @f;

    function val = f(t, x, u)
%   f = A1*x + A2*(x kron x/vx) + A3*(u kron x/vx) + A4*(x kron u/vu) + A5*(u kron u/vu) + A6*u + a
        % Compartment column
        v = B1 * x + B2 * u + b;
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

function handle = dfdxHidden(A1, A2, A3, A4, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdx;

    function val = dfdx(t, x, u)
%   dfdx = A1 + A2*(Ix kron x/vx) + A2*(x kron diag(1/vx)) + A3*(u kron diag(1/vx))
        % Compartment column
        v = B1 * x + B2 * u + b;
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

function handle = dfduHidden(A3, A4, A5, A6, B1, B2, b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @dfdu;

    function val = dfdu(t, x, u)
%   dfdu = A3*(Iu kron x/vx) + A4*(x kron diag(1/vu)) + A5*(Iu kron u/vu) + A5*(u kron diag(1/vu)) + A6
        % Compartment column
        v = B1 * x + B2 * u + b;
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

function handle = dfdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, dadk, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @dfdk;

    function val = dfdk(t, x, u)
%   dfdk = dA1dk*x + dA2dk*(x kron x/vx) + dA3dk*(u kron x/vx) + dA4dk*(x kron u/vu) + dA5dk*(u kron u/vu) + dA6dk*u + dadk
        % Compartment column
        v = B1 * x + B2 * u + b;
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

function handle = d2fdx2Hidden(A2, B1, B2, b, D2UsedColumns, D2UsedSpecies2, vxInd)
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
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        
        % Second-order derivative of (x kron x/vx)
        Ixkron1vx = sparse([D2UsedColumns; D2UsedColumns], [D2UsedColumns; D2UsedColumnsReverse], [vx(D2UsedSpecies2); vx(D2UsedSpecies2)], nx*nx,nx*nx);
        
        val = A2 * Ixkron1vx; % f_xx
        val = reshape(val, nx*nx, nx); % fx_x
    end
end

function handle = d2fdu2Hidden(A5, B1, B2, b, D5UsedColumns, D5UsedSpecies2, vuInd)
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
        v = B1 * x + B2 * u + b;
        vu = 1 ./ v(vuInd);
        
        % Second-order derivative of (u kron u/vu)
        Iukron1vu = sparse([D5UsedColumns; D5UsedColumns], [D5UsedColumns; D5UsedColumnsReverse], [vu(D5UsedSpecies2); vu(D5UsedSpecies2)], nu*nu,nu*nu);
        
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

function handle = d2fdudxHidden(A3, A4, B1, B2, b, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
Ix = speye(nx);
Iu = speye(nu);

% Return handle
handle = @d2fdudx;

    function val = d2fdudx(t, x, u)
%   d2fdudx = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1 - x *{xxx} dvxdx ^ -2)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1 - u *{uxu} dvudu ^ -2))
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = v(vxInd);
        vu = v(vuInd);
        vxinv = 1 ./ vx;
        vuinv = 1 ./ vu;
        dvxdx = B1(vxInd,:);
        dvudu = B2(vuInd,:);
        
        % Sparse and non-sparse kronecker multiplication
        Iukron1vx = sparse(D3UsedColumns, D3UsedColumns, vxinv(D3UsedSpecies2), nx*nu,nx*nu); % xu_xu
        if any(vec(dvxdx))
            Iukronxdvxdx = kron(Iu, bsxfun(@times, x .* vx .^ -2, dvxdx));
        else
            % Matlab's kron is expensive, don't do it if it is just zeros
            Iukronxdvxdx = sparse(nx*nu,nx*nu);
        end
        
        Ixkron1vu = sparse(D4UsedColumns, D4UsedColumns, vuinv(D4UsedSpecies2), nu*nx,nu*nx); % ux_ux
        if any(vec(dvudu))
            Ixkronudvudu = kron(Ix, bsxfun(@times, u .* vu .^ -2, dvudu));
        else
            % Matlab's kron is expensive, don't do it if it is just zeros
            Ixkronudvudu = sparse(nu*nx,nu*nx);
        end
        
        val = A3 * (Iukron1vx - Iukronxdvxdx) + spermute132(A4 * (Ixkron1vu - Ixkronudvudu), [nx,nu,nx],[nx,nx*nu]); % f_xu
        val = reshape(val, nx*nx,nu); % fx_u
    end
end

function handle = d2fdxduHidden(A3, A4, B1, B2, b, D3UsedColumns, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
Ix = speye(nx);
Iu = speye(nu);

% Return handle
handle = @d2fdxdu;

    function val = d2fdxdu(t, x, u)
%   d2fdxdu = A3 *{x.x;u.u} (Iu * (Ix *{xxx} vx ^ -1 - x *{xxx} dvxdx ^ -2)) + A4 *{x.x;u.u} (Ix * (Iu *{uxu} vu ^ -1 - u *{uxu} dvudu ^ -2))
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = v(vxInd);
        vu = v(vuInd);
        vxinv = 1 ./ vx;
        vuinv = 1 ./ vu;
        dvxdx = B1(vxInd,:);
        dvudu = B2(vuInd,:);
        
        % Sparse and non-sparse kronecker multiplication
        Iukron1vx = sparse(D3UsedColumns, D3UsedColumns, vxinv(D3UsedSpecies2), nx*nu,nx*nu); % xu_xu
        if any(vec(dvxdx))
            Iukronxdvxdx = kron(Iu, bsxfun(@times, x .* vx .^ -2, dvxdx));
        else
            % Matlab's kron is expensive, don't do it if it is just zeros
            Iukronxdvxdx = sparse(nx*nu,nx*nu);
        end
        
        Ixkron1vu = sparse(D4UsedColumns, D4UsedColumns, vuinv(D4UsedSpecies2), nu*nx,nu*nx); % ux_ux
        if any(vec(dvudu))
            Ixkronudvudu = kron(Ix, bsxfun(@times, u .* vu .^ -2, dvudu));
        else
            % Matlab's kron is expensive, don't do it if it is just zeros
            Ixkronudvudu = sparse(nu*nx,nu*nx);
        end
        
        val = spermute132(A3 * (Iukron1vx - Iukronxdvxdx), [nx,nx,nu],[nx,nu*nx]) + A4 * (Ixkron1vu - Ixkronudvudu); % f_ux
        val = reshape(val, nx*nu,nx); % fu_x
    end
end

function handle = d2fdkdxHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA1dk_fk_x, 1) / nx;

% Return handle
handle = @d2fdkdx;

    function val = d2fdkdx(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vx(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vx(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dA1dk_fk_x + dA2dk_fk_xx * (Ixkronxvx + xkron1vx) + dA3dk_fk_ux * ukron1vx + dA4dk_fk_xu * Ixkronuvu; % fx_k
        val = spermute132(val, [nx,nk,nx], [nx*nx,nk]); % fk_x 
    end
end

function handle = d2fdxdkHidden(dA1dk_fk_x, dA2dk_fk_xx, dA3dk_fk_ux, dA4dk_fk_xu, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdxdk;

    function val = d2fdxdk(t, x, u)
%   d2fdxdk = dA1dk + dA2dk*(Ix kron x/vx) + dA2dk*(x kron diag(1/vx)) + dA3dk*(u kron diag(1/vx)) + dA4dk*(Ix kron u/vu)
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vx(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vx(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dA1dk_fk_x + dA2dk_fk_xx * (Ixkronxvx + xkron1vx) + dA3dk_fk_ux * ukron1vx + dA4dk_fk_xu * Ixkronuvu; % fx_k
    end
end

function handle = d2fdkduHidden(dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, B1, B2, b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dA3dk_fk_ux, 1) / nx;

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        v = B1 * x + B2 * u + b;
        xvx = x ./ v(vxInd);
        vu = 1 ./ v(vuInd);
        uvu = u .* vu;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(D3UsedColumns, D3UsedSpecies1, xvx(D3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vu(D4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(D5UsedColumns, D5UsedSpecies1, uvu(D5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vu(D5UsedSpecies2), nu*nu,nu);
        
        val = dA3dk_fk_ux * Iukronxvx + dA4dk_fk_xu * xkron1vu + dA5dk_fk_uu * (Iukronuvu + ukron1vu) + dA6dk_fk_u; % fk_u
        val = spermute132(val, [nx,nk,nu], [nx*nu,nk]); % fu_k
    end
end

function handle = d2fdudkHidden(dA3dk_fk_ux, dA4dk_fk_xu, dA5dk_fk_uu, dA6dk_fk_u, B1, B2, b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2fdudk;

    function val = d2fdudk(t, x, u)
%   d2fdudk = 
        % Compartment column
        v = B1 * x + B2 * u + b;
        xvx = x ./ v(vxInd);
        vu = 1 ./ v(vuInd);
        uvu = u .* vu;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(D3UsedColumns, D3UsedSpecies1, xvx(D3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vu(D4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(D5UsedColumns, D5UsedSpecies1, uvu(D5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vu(D5UsedSpecies2), nu*nu,nu);
        
        val = dA3dk_fk_ux * Iukronxvx + dA4dk_fk_xu * xkron1vu + dA5dk_fk_uu * (Iukronuvu + ukron1vu) + dA6dk_fk_u; % fk_u
    end
end

function handle = rHidden(D1, D2, D3, D4, D5, D6, d, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @r;

    function val = r(t, x, u)
%   r = D1*x + D2*(x kron x/vx) + D3*(u kron x/vx) + D4*(x kron u/vu) + D5*(u kron u/vu) + D6*u + d
        % Compartment column
        v = B1 * x + B2 * u + b;
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

function handle = drdxHidden(D1, D2, D3, D4, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @drdx;

    function val = drdx(t, x, u)
%   drdx = D1 + D2*(Ix kron x/vx) + D2*(x kron diag(1/vx)) + D3*(u kron diag(1/vx))
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vx(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vx(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nx*nu,nx);
        
        val = D1 + D2 * (Ixkronxvx + xkron1vx) + D3 * ukron1vx + D4 * Ixkronuvu; % r_x
    end
end

function handle = drduHidden(D3, D4, D5, D6, B1, B2, b, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @drdu;

    function val = drdu(t, x, u)
%   drdu = D3*(Iu kron x/vx) + D4*(x kron diag(1/vu)) + D5*(Iu kron u/vu) + D5*(u kron diag(1/vu)) + D6
        % Compartment column
        v = B1 * x + B2 * u + b;
        xvx = x ./ v(vxInd);
        vu = 1 ./ v(vuInd);
        uvu = u .* vu;
        
        % Sparse kronecker multiplication
        Iukronxvx = sparse(D3UsedColumns, D3UsedSpecies1, xvx(D3UsedSpecies2), nx*nu,nu);
        xkron1vu  = sparse(D4UsedColumns, D4UsedSpecies2, x(D4UsedSpecies1) .* vu(D4UsedSpecies2), nu*nx,nu);
        Iukronuvu = sparse(D5UsedColumns, D5UsedSpecies1, uvu(D5UsedSpecies2), nu*nu,nu);
        ukron1vu  = sparse(D5UsedColumns, D5UsedSpecies2, u(D5UsedSpecies1) .* vu(D5UsedSpecies2), nu*nu,nu);
        
        val = D3 * Iukronxvx + D4 * xkron1vu + D5 * (Iukronuvu + ukron1vu) + D6; % r_u
    end
end

function handle = drdkHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, dD5dk_rk_uu, dD6dk_rk_u, dddk, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, D5UsedColumns, D5UsedSpecies1, D5UsedSpecies2, vxInd, vuInd)
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
        v = B1 * x + B2 * u + b;
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

function handle = d2rdx2Hidden(D2, B1, B2, b, D2UsedColumns, D2UsedSpecies2, vxInd)
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
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        
        % Second-order derivative of (x kron x/vx)
        Ixkron1vx = sparse([D2UsedColumns; D2UsedColumns], [D2UsedColumns; D2UsedColumnsReverse], [vx(D2UsedSpecies2); vx(D2UsedSpecies2)], nx*nx,nx*nx);
        
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

function handle = d2rdxdkHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);

% Return handle
handle = @d2rdxdk;

    function val = d2rdxdk(t, x, u)
%   d2rdxdk = dD1dk + dD2dk*(Ix kron x/vx) + dD2dk*(x kron diag(1/vx)) + dD3dk*(u kron diag(1/vx)) + dD4dk*(Ix kron u/vu)
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vx(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vx(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dD1dk_rk_x + dD2dk_rk_xx * (Ixkronxvx + xkron1vx) + dD3dk_rk_ux * ukron1vx + dD4dk_rk_xu * Ixkronuvu; % rx_k
    end
end

function handle = d2rdkdxHidden(dD1dk_rk_x, dD2dk_rk_xx, dD3dk_rk_ux, dD4dk_rk_xu, B1, B2, b, D2UsedColumns, D2UsedSpecies1, D2UsedSpecies2, D3UsedColumns, D3UsedSpecies1, D3UsedSpecies2, D4UsedColumns, D4UsedSpecies1, D4UsedSpecies2, vxInd, vuInd, nr)
% Constants
nx = numel(vxInd);
nu = numel(vuInd);
nk = size(dD1dk_rk_x, 1) / nr;

% Return handle
handle = @d2rdkdx;

    function val = d2rdkdx(t, x, u)
%   d2rdxdk = dD1dk + dD2dk*(Ix kron x/vx) + dD2dk*(x kron diag(1/vx)) + dD3dk*(u kron diag(1/vx)) + dD4dk*(Ix kron u/vu)
        % Compartment column
        v = B1 * x + B2 * u + b;
        vx = 1 ./ v(vxInd);
        xvx = x .* vx;
        uvu = u ./ v(vuInd);
        
        % Sparse kronecker multiplication
        Ixkronxvx = sparse(D2UsedColumns, D2UsedSpecies1, xvx(D2UsedSpecies2), nx*nx,nx);
        xkron1vx  = sparse(D2UsedColumns, D2UsedSpecies2, x(D2UsedSpecies1) .* vx(D2UsedSpecies2), nx*nx,nx);
        ukron1vx  = sparse(D3UsedColumns, D3UsedSpecies2, u(D3UsedSpecies1) .* vx(D3UsedSpecies2), nx*nu,nx);
        Ixkronuvu = sparse(D4UsedColumns, D4UsedSpecies1, uvu(D4UsedSpecies2), nu*nx,nx);
        
        val = dD1dk_rk_x + dD2dk_rk_xx * (Ixkronxvx + xkron1vx) + dD3dk_rk_ux * ukron1vx + dD4dk_rk_xu * Ixkronuvu; % rx_k
        val = spermute132(val, [nr,nk,nx], [nr*nx,nk]); % rk_x 
    end
end