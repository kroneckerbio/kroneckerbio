function m = addReactionMassActionAmount(m, name, reactants, products, kForward, kReverse)
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
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up
if nargin < 6
    kReverse = [];
end

% Standardize names
[name1, name2] = fixReactionName(name, kForward, kReverse);
reactants = fixReactionSpecies(reactants);
products = fixReactionSpecies(products);
kForward = fixReactionParameter(kForward);
kReverse = fixReactionParameter(kReverse);

% Increment counter
if ~isempty(kForward{1})
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactions(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name1;
    
    m.add.Reactions(nr).Reactants = reactants;
    
    m.add.Reactions(nr).Products = products;
    
    m.add.Reactions(nr).Parameter = kForward;
    
    m.Ready = false;
end

if ~isempty(kReverse{1})
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactions(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name2;
    
    m.add.Reactions(nr).Reactants = products;
    
    m.add.Reactions(nr).Products = reactants;
    
    m.add.Reactions(nr).Parameter = kReverse;
    
    m.Ready = false;
end
