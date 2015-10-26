function AddReaction(this, name, reactants, products, forward, reverse, compartment)
%AddReaction Add a reaction to a AnalyticModel
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
%   m: [ scalar AnalyticModel ]
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

% Clean-up
if nargin < 7
    compartment = [];
    if nargin < 6
        reverse = [];
    end
end

if isempty(compartment)
    compartment = '';
end

% Standardize reaction name
[nameForward, nameReverse] = FieldValidator.ReactionName(name);

% Standardize reactions and products
reactants = FieldValidator.ReactionSpecies(reactants);
products  = FieldValidator.ReactionSpecies(products);

% Standardize reaction rate expressions - do nothing - just paste in

% Add separate reactions for forward and reverse (if applicable)
if ~isempty(forward)
    nr = this.nr + 1;
    this.nr = nr;
    this.growReactions;
    
    this.Reactions(nr).Name        = nameForward;
    this.Reactions(nr).Reactants   = reactants;
    this.Reactions(nr).Products    = products;
    this.Reactions(nr).Rate        = forward;
    this.Reactions(nr).Compartment = compartment;
    
    this.Ready = false;
end

if ~isempty(reverse)
    nr = this.nr + 1;
    this.nr = nr;
    this.growReactions;
    
    this.Reactions(nr).Name        = nameReverse;
    this.Reactions(nr).Reactants   = products;
    this.Reactions(nr).Products    = reactants;
    this.Reactions(nr).Rate        = reverse;
    this.Reactions(nr).Compartment = compartment;
    
    this.Ready = false;
end
