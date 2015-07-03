function m = addReactionMassActionAmount(m, name, reactant1, reactant2, product1, product2, kForward, kReverse)
%AddReaction Add a reaction to a KroneckerBio model
%
%   m = AddReaction(m, name, reactant1, reactant2, product1,
%                   product2, kForward, kReverse)
%
%   A reaction is a conversion of one or two reactant species into one or
%   two product species associated with a defined rate constant. The
%   reverse reaction may or may not be specified at the same time. Each
%   reactant and product must be a single species in a single compartment
%   (e.g. cytoplasm.glucose).
%
%   Since most reactions happen in a single compartment, Kronecker does not
%   make you type out the compartment for every species; you can specify a
%   default compartment, which will be appended to each of the species that
%   do not have a compartment specified. Multiple compartments can be
%   specified and FinalizeModel will look in those compartments for the
%   species name. If no compartment is specified, it will look in every
%   compartment for that species name. If a species with that name is found
%   in several compartments that are searched, an error will be thrown.
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
%   reactant1: [ string ]
%       This is a species full name in the model or just a species name
%       that can accept a string in compartment as its compartment.
%   reactant2: [ string ]
%       Like reactant1.
%   product1: [ string ]
%       Like reactant1.
%   product2: [ string ]
%       Like reactant1.
%   kForward: [ string | cell vector {{'', 1}} ]
%       This is the kinetic parameter of this reaction along with a scaling
%       factor by which the kinetic parameter will be multiplied. It can be
%       provided as a cell vector of length 2, in which case the first
%       element is a string name of the parameter and the second element is
%       the scaling factor. If only a string is provided, it is the name of
%       the kinetic parameter and the scaling factor is assumed to be 1.
%   kReverse: [ string | cell vector {{'', 1}} ]
%       Like kForward.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up
if nargin < 8
    kReverse = [];
    if nargin < 7
        kForward = [];
    end
end

% Standardize names
[name1, name2] = fixReactionName(name, kForward, kReverse);
reactant1   = fixReactionSpecies(reactant1);
reactant2   = fixReactionSpecies(reactant2);
product1    = fixReactionSpecies(product1);
product2    = fixReactionSpecies(product2);
kForward    = fixReactionParameter(kForward);
kReverse    = fixReactionParameter(kReverse);

% Increment counter
if ~isempty(kForward{1})
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactions(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name1;
    
    m.add.Reactions(nr).Reactants = cell(1,0);
    if ~isempty(reactant1)
        m.add.Reactions(nr).Reactants = {reactant1};
    end
    if ~isempty(reactant2)
        m.add.Reactions(nr).Reactants = [m.add.Reactions(nr).Reactants, {reactant2}];
    end
    
    m.add.Reactions(nr).Products = cell(1,0);
    if ~isempty(product1)
        m.add.Reactions(nr).Products = {product1};
    end
    if ~isempty(product2)
        m.add.Reactions(nr).Products = [m.add.Reactions(nr).Products, {product2}];
    end
    
    m.add.Reactions(nr).Parameter = kForward;
    
    m.Ready = false;
end

if ~isempty(kReverse{1})
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactions(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name2;
    
    m.add.Reactions(nr).Reactants = cell(1,0);
    if ~isempty(product1)
        m.add.Reactions(nr).Reactants = {product1};
    end
    if ~isempty(product2)
        m.add.Reactions(nr).Reactants = [m.add.Reactions(nr).Reactants, {product2}];
    end
    
    m.add.Reactions(nr).Products = cell(1,0);
    if ~isempty(reactant1)
        m.add.Reactions(nr).Products = {reactant1};
    end
    if ~isempty(reactant2)
        m.add.Reactions(nr).Products = [m.add.Reactions(nr).Products, {reactant2}];
    end
    
    m.add.Reactions(nr).Parameter = kReverse;
    
    m.Ready = false;
end
