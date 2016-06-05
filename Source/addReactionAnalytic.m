function m = addReactionAnalytic(m, name, reactants, products, forward, reverse, compartment)
%AddReaction Add a reaction to a Model.Analytic
%
%   m = AddReaction(m, name, reactants, products, forward, reverse, compartment)
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
%   forward: [ string ]
%       Expression for forward reaction rate
%   reverse: [ string | empty {''} ]
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

if isempty(reverse)
    reverse = '';
end

if isempty(compartment)
    compartment = '';
end

% Standardize reaction name
[nameForward, nameReverse] = fixReactionName(name);

% Standardize reactions and products
reactants = fixReactionSpecies(reactants);
products  = fixReactionSpecies(products);

% Standardize reaction rate expressions
forward = fixRateExpressionAnalytic(forward);
reverse = fixRateExpressionAnalytic(reverse);

% Add separate reactions for forward and reverse (if applicable)
if ~isempty(forward)
    nr = m.nr + 1;
    m.nr = nr;
    m.Reactions = growReactionsAnalytic(m.Reactions, nr);
    
    m.Reactions(nr).Name        = nameForward;
    m.Reactions(nr).Reactants   = reactants;
    m.Reactions(nr).Products    = products;
    m.Reactions(nr).Rate        = forward;
    m.Reactions(nr).Compartment = compartment;
    
    m.Ready = false;
end

if ~isempty(reverse)
    nr = m.nr + 1;
    m.nr = nr;
    m.Reactions = growReactionsAnalytic(m.Reactions, nr);
    
    m.Reactions(nr).Name        = nameReverse;
    m.Reactions(nr).Reactants   = products;
    m.Reactions(nr).Products    = reactants;
    m.Reactions(nr).Rate        = reverse;
    m.Reactions(nr).Compartment = compartment;
    
    m.Ready = false;
end
