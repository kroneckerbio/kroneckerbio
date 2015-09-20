function m = addReactionAnalytic(m, name, reactants, products, forward, reverse, compartment)
%AddReaction Add a reaction to a Model.Analytic
%
%   m = AddReaction(m, name, reactants, products, kForward, kReverse, compartment)
%
%   A reaction is a conversion of zero or more reactant species into zero
%   or more product species associated with a rate expression. The reverse
%   reaction may or may not be specified at the same time. Each reactant
%   and product must be a single species in a single compartment (e.g.
%   cytoplasm.glucose). An optional reaction compartment may be supplied to
%   disambiguate species without a specified compartment.
%
%   The species and parameters must be added with the AddState, AddInput,
%   and AddParameter functions. They can be added in any order as long as
%   they are present before the next call to FinalizeModel.
%
%   Inputs
%   m: [ Model.Analytic struct ]
%       The model to which the reaction will be added
%   name: [ string | cell vector 2 of strings ]
%       A name for the reaction, or two different names for the forward and
%       reverse reactions. Reaction names do not need to be unique, so if
%       only one name is provided, it will be identical for the forward and
%       reverse reactions.
%   reactants: [ cell array of strings | string | empty ]
%       A list of species names consumed by this reaction. Must be a
%       species full name or a name that is unique among species in the
%       model, unless a reaction compartment is supplied.
%   products: [ cell array of strings | string | empty ]
%       A list of species names produced by this reaction. Like reactants.
%   kForward: [ string ]
%       Expression for forward reaction rate
%   kReverse: [ string | empty {''} ]
%       Expression for reverse reaction rate
%   compartment: [ string | empty {''} ]
%       A compartment name.
%       If empty, then each species referred to in the reactants, products,
%       forward, and reverse must unambiguously refer to a single species
%       in the model, either because it is unique or because it is a
%       species full name.
%       If supplied, then each species in the reactants and products that
%       does not have a compartment will be disambiguated with this
%       compartment. The reactants and products that are disambiguated will
%       also be disambiguated in forward and reverse expressions.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% Clean-up
if nargin < 7
    compartment = [];
    if nargin < 6
        reverse = [];
    end
end

% Set defaults
if isempty(compartment)
    compartment = '';
end

% Standardize reaction name
[name1, name2] = fixReactionName(name, forward, reverse);

% Standardize species names into compartment.species
reactants = fixReactionSpecies(reactants);
products =  fixReactionSpecies(products);

% Standardize reaction rate expressions
forward = fixRateExpressionAnalytic(forward);
reverse = fixRateExpressionAnalytic(reverse);

% Standardize compartment name
compartment = fixCompartmentName(compartment);

% Incorporate reaction compartment
if ~isempty(compartment)
    all_species = [reactants, products];
    
    unqualified = false(size(all_species));
    unambiguous_species = all_species;
    for i = 1:numel(unambiguous_species)
        if ~ismember('.', all_species{i})
            % It is unqualified
            unqualified(i) = true;
            unambiguous_species{i} = [compartment '.' all_species{i}];
        end
    end
    
    % Update reactants and products
    reactants = unambiguous_species(1:numel(reactants));
    products = unambiguous_species(numel(reactants)+(1:numel(products)));
    
    % Rename in expressions
    forward = substituteQuotedExpressions(forward, all_species(unqualified), unambiguous_species(unqualified), true);
    reverse = substituteQuotedExpressions(reverse, all_species(unqualified), unambiguous_species(unqualified), true);
end

% Add separate reactions for forward and reverse (if applicable)
if ~isempty(forward)
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactionsAnalytic(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name1;
    
    m.add.Reactions(nr).Reactants = reactants;
    
    m.add.Reactions(nr).Products = products;
    
    m.add.Reactions(nr).Rate = forward;
    
    m.Ready = false;
end

if ~isempty(reverse)
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactionsAnalytic(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name2;
    
    m.add.Reactions(nr).Reactants = products;
    
    m.add.Reactions(nr).Products = reactants;
    
    m.add.Reactions(nr).Rate = reverse;
    
    m.Ready = false;
end
