function qualifyNames(this)
% Qualify species with compartments in all fields where unqualified species
% appear. Later steps expect this.

import Parser.*

nr = this.nr;
ny = this.ny;
nz = this.nz;

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
    rate        = reaction.Rate;
    compartment = reaction.Compartment;
    
    % Make sure compartment exists in model if specified
    if ~isempty(compartment)
        assert(ismember(compartment, m.v_names), 'KroneckerBio:FinalizeModel:MissingReactionCompartment', 'Compartment %s not found in reaction %s', compartment, name)
    end
    
    % Get species that appear in reaction, including species that aren't reactants or products 
    species = [reactants, products, getSpeciesFromExpr(rate, xu_full_names)];
    
    % Incorporate reaction compartment
    [qualifiedSpecies, qualified] = qualifyCompartments(species, xu_full_names, compartment);
    
    % Update reactants and products
    nReactants = numel(reactants);
    nProducts = numel(products);
    reactants = qualifiedSpecies(1:nReactants);
    products = qualifiedSpecies(nReactants+1:nReactants+nProducts);
    
    % Update species in rate expression
    rate = substituteQuotedExpressions(rate, species(~qualified), qualifiedSpecies(~qualified), true);
    
    % Update reaction
    reaction.Reactants = reactants;
    reaction.Products  = products;
    reaction.Rate      = rate;
    this.Reactions(i)  = reaction;
    
end

%% Outputs
for i = 1:ny
    
    % Extract output
    output = this.Outputs(i);
    expr = output.Expression;
    
    % Check/qualify species in output expression
    species = getSpeciesFromExpr(expr, xu_full_names);
    [qualifiedSpecies, qualified] = qualifyCompartments(species, xu_full_names);
    expr = substituteQuotedExpressions(expr, species(~qualified), qualifiedSpecies(~qualified), true);
    
    % Update output
    output.Expression = expr;
    this.Outputs(i) = output;
    
end

%% Rules
for i = 1:nz
    
    % Extract rule
    rule = this.Rules(i);
    name = rule.Name;
    expr = rule.Expression;
    
    % Check/qualify species in rule name
    species = getSpeciesFromExpr(name, xu_full_names);
    [qualifiedSpecies, qualified] = qualifyCompartments(species, xu_full_names);
    name = substituteQuotedExpressions(name, species(~qualified), qualifiedSpecies(~qualified)); % don't add quotes to output
    
    % Check/qualify species in rule expression
    species = getSpeciesFromExpr(expr, xu_full_names);
    [qualifiedSpecies, qualified] = qualifyCompartments(species, xu_full_names);
    expr = substituteQuotedExpressions(expr, species(~qualified), qualifiedSpecies(~qualified), true);
    
    % Update rule
    rule.Name       = name;
    rule.Expression = expr;
    this.Rules(i) = rule;
    
end