function m = FinalizeModel(m, varargin)
%FinalizeModel Update the mathematical components of the model to reflect
%   changes made to the model.
%
%   m = FinalizeModel(m, ...)
%
%   This function contains generic model processing steps. Select the 
%   appropriate additional help file below depending on the type of the model.
%
%   Model.MassActionAmount
%       help finalizeModelMassActionAmount
%
%   Model.Analytic
%       help finalizeModelAnalytic

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
assert(nargin >= 1, 'KroneckerBio:FinalizeModel:TooFewInputs', 'FinalizeModel requires at least 1 input argument')
assert(isscalar(m), 'KroneckerBio:FinalizeModel:MoreThanOneModel', 'The model structure must be scalar')

%% Common finalization
% Model-specific common finalization modifications
if is(m, 'Model.MassActionAmount')
    rateName = 'Parameter';
elseif is(m, 'Model.Analytic')
    rateName = 'Rate';
else
    error('KroneckerBio:AddOutput:m', 'm must be a model')
end

%% Trim m.* components to only those actually added
m.Compartments = m.Compartments(1:m.nv);
m.Seeds        = m.Seeds(1:m.ns);
m.Parameters   = m.Parameters(1:m.nk);
m.Inputs       = m.Inputs(1:m.nu);
m.States       = m.States(1:m.nx);
m.Reactions    = m.Reactions(1:m.nr);
m.Rules        = m.Rules(1:m.nz);
m.Outputs      = m.Outputs(1:m.ny);

nv = numel(m.Compartments);
nk = numel(m.Parameters);
ns = numel(m.Seeds);
nu = numel(m.Inputs);
nx = numel(m.States);
nr = numel(m.Reactions);
nz = numel(m.Rules);
ny = numel(m.Outputs);
nxu = nx + nu;

%% Extract names
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

% Make species - compartment index mapping
vu_names = vec({m.Inputs.Compartment});
vuInd = lookupmember(vu_names, v_names);
vx_names = vec({m.States.Compartment});
vxInd = lookupmember(vx_names, v_names);
m.vxInd = vxInd;
m.vuInd = vuInd;

%% Error for repeated components
[~, ia, ~] = unique(v_names);
v_repeated = v_names(setdiff(1:nv, ia));
assert(isempty(v_repeated), 'KroneckerBio:FinalizeModel:RepeatCompartment', 'Compartment %s not unique', cellstr2str(v_repeated))

[~, ia, ~] = unique(k_names);
k_repeated = k_names(setdiff(1:nk, ia));
assert(isempty(k_repeated), 'KroneckerBio:FinalizeModel:RepeatParameter', 'Parameter %s not unique', cellstr2str(k_repeated))

[~, ia, ~] = unique(s_names);
s_repeated = s_names(setdiff(1:ns, ia));
assert(isempty(s_repeated), 'KroneckerBio:FinalizeModel:RepeatSeed', 'Seed %s not unique', cellstr2str(s_repeated))

[~, ia, ~] = unique(xu_full_names);
xu_repeated = xu_full_names(setdiff(1:nxu, ia));
assert(isempty(xu_repeated), 'KroneckerBio:FinalizeModel:RepeatSpecies', 'Species %s not unique', cellstr2str(xu_repeated))

[~, ia, ~] = unique(y_names);
y_repeated = y_names(setdiff(1:ny, ia));
assert(isempty(y_repeated), 'KroneckerBio:FinalizeModel:RepeatOutput', 'Output %s not unique', cellstr2str(y_repeated))

% Does it make sense to exclude repeated reactions and rules?

%% Resolve species compartments and standardize names

    function matches = getSpeciesFromExpr(expr)
        % Needed to get species that appear in rate expression but aren't
        %   reactants or products
        % Looks for unqualified species and qualified compartment.species.
        % Names in expr must be double-quoted if they contain invalid
        %   characters. Note: this implies qualified compartment.species must
        %   always be quoted because the dot is an invalid character.
        
        % Make species lookup list, quoting things that contain invalid characters
        namesLookup = [xu_names; xu_full_names];
        for iNames = 1:length(namesLookup)
            names_i = namesLookup{iNames};
            if regexp(names_i, '[^\w.]')
                namesLookup{iNames} = ['"' names_i '"'];
            end
        end
        
        % Tokenize expression into potentially substitutable parts and see if
        %   they're valid species names
        matches = {};
        parts = regexp(expr, '[\w\.]+|"[^"]+"', 'match');
        for iPart = 1:length(parts)
            part = parts{iPart};
            if any(ismember(namesLookup, part))
                part = strrep(part, '"', ''); % strip double-quotes because parts will go into species lists, which don't require them
                matches = [matches, part];
            end
        end
    end

    function [unambiguousSpecies, unqualified] = qualifyCompartment(species_all, compartment)
        if nargin < 2
            compartment = [];
        end
        
        % Incorporate reaction compartment
        unqualified = false(size(species_all));
        unambiguousSpecies = species_all;
        for j = 1:numel(unambiguousSpecies)
            
            species = species_all{j};
            
            if ismember('.', species) % qualified - make sure this species exists
                
                assert(ismember(species, xu_full_names), 'KroneckerBio:FinalizeModel:MissingQualifiedSpeciesName', 'The qualified name %s does not exist in the model', species)
                
            else % unqualified - apply default compartment rules
                
                assert(ismember(species, xu_names), 'KroneckerBio:FinalizeModel:MissingUnqualifiedSpeciesName', 'The unqualified name %s does not exist in the model', species)
                
                unqualified(j) = true;
                
                if isempty(compartment) % no reaction compartment
                    % Make sure species is unique - if not, it's ambiguous and throws an error
                    speciesPos = ismember(xu_names, species);
                    assert(sum(speciesPos) == 1, 'KroneckerBio:FinalizeModel:AmbiguousSpeciesName', 'The species name %s appears in multiple compartments in %s and no default compartment is specified', species, name)
                    unambiguousSpecies{j} = xu_full_names{speciesPos};
                else % reaction compartment present
                    unambiguousSpecies{j} = [compartment '.' species];
                    assert(ismember(unambiguousSpecies{j}, xu_full_names), 'KroneckerBio:FinalizeModel:MissingSpeciesInReactionCompartment', 'The species %s not found in compartment %s', species, compartment)
                end
                
            end
            
        end
    end

for i = 1:nr
    
    % Extract reaction
    reaction = m.Reactions(i);
    name        = reaction.Name;
    reactants   = reaction.Reactants;
    products    = reaction.Products;
    rate        = reaction.(rateName);
    compartment = reaction.Compartment;
    
    % Make sure compartment exists in model if specified
    if ~isempty(compartment)
        assert(ismember(compartment, v_names), 'KroneckerBio:FinalizeModel:MissingReactionCompartment', 'Compartment %s not found in reaction %s', compartment, name)
    end
    
    % Get species that appear in reaction
    %   In analytic models, this includes species that aren't reactants or products 
    reactionSpecies = [reactants, products];
    if is(m, 'Model.Analytic')
        reactionSpecies = [reactionSpecies, getSpeciesFromExpr(rate)];
    end
    
    % Incorporate reaction compartment
    [unambiguousSpecies, unqualified] = qualifyCompartment(reactionSpecies, compartment);
    
    % Update reactants and products
    reactants = unambiguousSpecies(1:numel(reactants));
    products = unambiguousSpecies(numel(reactants)+(1:numel(products)));
    
    % Rename in analytic model rate expression
    if is(m, 'Model.Analytic')
        rate = substituteQuotedExpressions(rate, reactionSpecies(unqualified), unambiguousSpecies(unqualified), true);
    end
    
    % Update reaction
    reaction.Reactants  = reactants;
    reaction.Products   = products;
    reaction.(rateName) = rate;
    m.Reactions(i) = reaction;
    
end

%% Resolve compartment sizes and standardize names
for i = 1:nv
    % Extract compartment
    compartment = m.Compartments(i);
    name = compartment.Name;
    expr = compartment.Size;
    
    % Check/qualify species in output expression
    if is(m, 'Model.MassActionAmount')
        contributor_names = expr(:,1);
        non_empty_contributors = ~cellfun(@isempty, contributor_names);
        unambiguous_names = qualifyCompartment(contributor_names(non_empty_contributors));
        expr(non_empty_contributors,1) = unambiguous_names;
    elseif is(m, 'Model.Analytic')
        contributor_names = getSpeciesFromExpr(expr);
        [unambiguousSpecies, unqualified] = qualifyCompartment(contributor_names);
        expr = substituteQuotedExpressions(expr, contributor_names(unqualified), unambiguousSpecies(unqualified), true);
    end
    
    % Update reaction
    compartment.Size = expr;
    m.Compartments(i) = compartment;
end

%% Resolve output compartments and standardize names
for i = 1:ny
    % Extract output
    output = m.Outputs(i);
    name = output.Name;
    expr = output.Expression;
    
    % Check/qualify species in output expression
    if is(m, 'Model.MassActionAmount')
        contributor_names = expr(:,1);
        non_empty_contributors = ~cellfun(@isempty, contributor_names);
        unambiguous_names = qualifyCompartment(contributor_names(non_empty_contributors));
        expr(non_empty_contributors,1) = unambiguous_names;
    elseif is(m, 'Model.Analytic')
        contributor_names = getSpeciesFromExpr(expr);
        [unambiguousSpecies, unqualified] = qualifyCompartment(contributor_names);
        expr = substituteQuotedExpressions(expr, contributor_names(unqualified), unambiguousSpecies(unqualified), true);
    end
    
    % Update reaction
    output.Expression = expr;
    m.Outputs(i) = output;
end

%% Type-specific finalization
if is(m, 'Model.MassActionAmount')
    m = finalizeModelMassActionAmount(m, varargin{:});
elseif is(m, 'Model.Analytic')
    m = finalizeModelAnalytic(m, varargin{:});
end
end
