function AddReaction(this, name, reactants, products, kForward, kReverse, compartment)
%AddReaction Add a reaction to a MassActionAmountModel
%
%   A reaction is a conversion of zero, one, or two reactant species into
%   zero or more product species associated with a defined rate constant.
%   The reverse reaction may or may not be specified at the same time. If
%   there are more than two products then the reverse reaction must not be
%   specified.
%
%   The species and parameters must be added with the AddState, AddInput,
%   and AddParameter functions. They can be added in any order as long as
%   they are present before the next call to FinalizeModel.
%
%   Inputs
%   m: [ scalar MassActionAmountModel ]
%       The model to which the reaction will be added
%   name: [ string | cell vector 2 of strings ]
%       A name for the reaction, or two different names for the forward and
%       reverse reactions. Reaction names do not need to be unique, so if
%       only one name is provided, it will be identical for the forward and
%       reverse reactions.
%   reactants: [ cell array of strings | string | empty ]
%       A list of species in the model
%   products: [ string ]
%       Like reactants
%   kForward: [ string | cell vector ]
%       This is the kinetic parameter of this reaction along with a scaling
%       factor by which the kinetic parameter will be multiplied. It can be
%       provided as a cell vector of length 2, in which case the first
%       element is a string name of the parameter and the second element is
%       the scaling factor. If only a string is provided, it is the name of
%       the kinetic parameter and the scaling factor is assumed to be 1.
%   kReverse: [ string | cell vector | empty {''} ]
%       Like kForward.
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
        kReverse = [];
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

% Standardize reaction rate parameters
kForward  = FieldValidator.ReactionParameter(kForward);
kReverse  = FieldValidator.ReactionParameter(kReverse);

% Add separate reactions for forward and reverse (if applicable)
if ~isempty(kForward{1})
    nr = this.nr + 1;
    this.nr = nr;
    this.growReactions;
    
    this.Reactions(nr).Name        = nameForward;
    this.Reactions(nr).Reactants   = reactants;
    this.Reactions(nr).Products    = products;
    this.Reactions(nr).Parameter   = kForward;
    this.Reactions(nr).Compartment = compartment;
    
    this.Ready = false;
end

if ~isempty(kReverse{1})
    nr = this.nr + 1;
    this.nr = nr;
    this.growReactions;
    
    this.Reactions(nr).Name        = nameReverse;
    this.Reactions(nr).Reactants   = products;
    this.Reactions(nr).Products    = reactants;
    this.Reactions(nr).Parameter   = kReverse;
    this.Reactions(nr).Compartment = compartment;
    
    this.Ready = false;
end
