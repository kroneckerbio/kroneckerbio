function m = AddReaction(m, name, compartment, reactant1, reactant2, product1, product2, kForward, kReverse)
%AddReaction Add a reaction to a KroneckerBio model
%
%   m = AddReaction(m, name, compartment, reactant1, reactant2, product1,
%                   product2, kForward, kReverse)
%
%   A reaction is conversion of one or two reactant species into one or two
%   product species associated with a defined rate parameter. The reverse
%   reaction may or may not be specified at the same time. Each reactant
%   and product must be specified in a particular compartment by using the
%   species full name (e.g. cytoplasm.glucose). 
%
%   Since most reactions happen in a single compartment, Kronecker does not
%   make you type out the compartment for every species; you can specify a
%   default compartment, which will be appended to each of the species that
%   do not have a compartment specified. Multiple compartments can be
%   specified and the reaction will occur in all of those compartments (if
%   they have the reatants). If no compartment is specified, it will be
%   assumed that the reaction happens in every compartment in the model.
%
%   The species and parameters must be added with the AddSpecies and
%   AddParameter functions. They can be added in any order as long as they
%   are present before the next call to FinalizeModel.
%
%   In FinalizeModel, reactions that do not have the reactants required are
%   silently ignored. To say that glucose binds to glucokinase is true even
%   if there is no glucose in the compartment. Reactions that have the
%   reactants required but not the products are an error because the
%   reaction is said to occur but there is nowhere to put the results.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the reaction will be added
%   name: [ string | cell array 2 of strings {''} ]
%       A name for the reaction, or two different names for the forward and
%       reverse reactions. Reaction names do not need to be unique, so if
%       only one name is provided, it will be idnetical for the forward and
%       reverse reactions.
%   compartment: [ string | cell array of strings {''} ]
%       For all species that does not have a compartment specified, this
%       compartment will be used as their compartments. If a cell array of
%       strings is provided, a reaction will be added for each compartment
%       given. If the string is empty, a reaction will be added to every
%       compartment in the model. Note that if the compartment does not
%       have the reatants for the reaction, the reaction will be ignored.
%   reactant1: [ string ]
%       This is a species full name in the model or just a species name
%       that can accept a string in compartment as its compartment.
%   reactant2: [ string ]
%       Like reactant1.
%   product1: [ string ]
%       Like reactant1.
%   product2: [ string ]
%       Like reactant1.
%   kForward: [ string ]
%       This is a rate parameter name in the model.
%   kReverse: [ string ]
%       Like kForward.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up
if nargin < 9
    kReverse = [];
    if nargin < 8
        kForward = [];
    end
end

% Standardize names
[name1 name2] = fixReactionName(name, kForward, kReverse);
reactant1 = fixReactionSpecies(reactant1);
reactant2 = fixReactionSpecies(reactant2);
product1  = fixReactionSpecies(product1);
product2  = fixReactionSpecies(product2);
kForward  = fixReactionParameter(kForward);
kReverse  = fixReactionParameter(kReverse);
compartment = fixReactionCompartment(compartment);

% Increment counter
nr = m.add.nr + 1;
m.add.nr = nr;
m.add.Reactions = growReactions(m.add.Reactions, nr);

% Add item
m.add.Reactions(nr).Names = {name1, name2};
m.add.Reactions(nr).Compartment  = compartment;

if ~isempty(reactant1)
    m.add.Reactions(nr).Reactants = {reactant1, reactant2};
else
    % Prevent an order of {0, reactant}
    m.add.Reactions(nr).Reactants = {reactant2, reactant1};
end

if ~isempty(product1)
    m.add.Reactions(nr).Products = {product1, product2};
else
    % Prevent an order of {0, product}
    m.add.Reactions(nr).Products = {product2, product1};
end

m.add.Reactions(nr).Parameters = {kForward, kReverse};

m.Ready = false;