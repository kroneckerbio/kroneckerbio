function m = addReactionMassActionAmount(m, name, reactants, products, kForward, kReverse, compartment)
%AddReaction Add a reaction to a KroneckerBio model
%
%   m = AddReaction(m, name, reactants, products, kForward, kReverse)
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
%   m: [ model struct scalar ]
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
%       If empty, then each species referred to in the reactants and products
%       must unambiguously refer to a single species in the model, either
%       because it is unique or because it is a species full name. If supplied,
%       then each species in the reactants and products that does not have a
%       compartment will be disambiguated with this compartment.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

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
[nameForward, nameReverse] = fixReactionName(name);

% Standardize reactions and products
reactants = fixReactionSpecies(reactants);
products  = fixReactionSpecies(products);

% Standardize reaction rate parameters
kForward  = fixReactionParameter(kForward);
kReverse  = fixReactionParameter(kReverse);

% Add separate reactions for forward and reverse (if applicable)
if ~isempty(kForward{1})
    assert(numel(reactants) <= 2, 'KroneckerBio:AddReaction:ToomanyReactants', 'A reaction may have at most two reactants')
    nr = m.nr + 1;
    m.nr = nr;
    m.Reactions = growReactions(m.Reactions, nr);
    
    m.Reactions(nr).Name        = nameForward;
    m.Reactions(nr).Reactants   = reactants;
    m.Reactions(nr).Products    = products;
    m.Reactions(nr).Parameter   = kForward;
    m.Reactions(nr).Compartment = compartment;
    
    m.Ready = false;
end

if ~isempty(kReverse{1})
    assert(numel(products) <= 2, 'KroneckerBio:AddReaction:ToomanyReactants', 'A reaction may have at most two reactants')
    nr = m.nr + 1;
    m.nr = nr;
    m.Reactions = growReactions(m.Reactions, nr);
    
    m.Reactions(nr).Name        = nameReverse;
    m.Reactions(nr).Reactants   = products;
    m.Reactions(nr).Products    = reactants;
    m.Reactions(nr).Parameter   = kReverse;
    m.Reactions(nr).Compartment = compartment;
    
    m.Ready = false;
end
