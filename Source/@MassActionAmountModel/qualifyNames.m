function qualifyNames(this)
% Qualify species with compartments in all fields where unqualified species
% appear. Later steps expect this.

import Parser.*

nr = this.nr;
ny = this.ny;

% Make list of all compartment.species in model
x_full_names = vec(strcat({this.States.Compartment}, '.', {this.States.Name}));
u_full_names = vec(strcat({this.Inputs.Compartment}, '.', {this.Inputs.Name}));
xu_full_names = [x_full_names; u_full_names];

%% Reactions with reaction compartments
for i = 1:nr
    
    % Extract reaction
    reaction = this.Reactions(i);
    name        = reaction.Name;
    reactants   = reaction.Reactants;
    products    = reaction.Products;
    rate        = reaction.Parameter;
    compartment = reaction.Compartment;
    
    % Make sure compartment exists in model if specified
    if ~isempty(compartment)
        assert(ismember(compartment, m.v_names), 'KroneckerBio:FinalizeModel:MissingReactionCompartment', 'Compartment %s not found in reaction %s', compartment, name)
    end
    
    % Get species that appear in reaction
    species = [reactants, products];
    
    % Incorporate reaction compartment
    qualifiedSpecies = qualifyCompartments(species, xu_full_names, compartment);
    
    % Update reactants and products
    nReactants = numel(reactants);
    nProducts = numel(products);
    reactants = qualifiedSpecies(1:nReactants);
    products = qualifiedSpecies(nReactants+1:nReactants+nProducts);
    
    % Update reaction
    reaction.Reactants = reactants;
    reaction.Products  = products;
    reaction.Parameter = rate;
    this.Reactions(i)  = reaction;
    
end

%% Outputs
for i = 1:ny
    
    % Extract output
    output = this.Outputs(i);
    expr = output.Expression;
    
    % Check/qualify species in output expression
    species = expr(:,1)';
    [qualifiedSpecies, qualified] = qualifyCompartments(species, xu_full_names, compartment);
    expr(~qualified,1) = qualifiedSpecies(~qualified)';
    
    % Update reaction
    output.Expression = expr;
    this.Outputs(i) = output;
    
end
