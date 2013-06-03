function m = AddLargeSizeReaction(m, name, compartment, reactants, products, kForward)
%AddReaction Add a reaction with an unlimited number of products to a
%   KroneckerBio model
%
%   m = AddLargeSizeReaction(m, name, compartment, reactants, products,
%                   kForward)
%
%   A reaction is conversion of one or two reactant species into an
%   arbitrary number of product species associated with a defined rate
%   parameter. The reverse reaction may or may not be specified at the same
%   time. Each reactant and product must be specified in a particular
%   compartment by using the species full name (e.g. cytoplasm.glucose).
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
%   reactants: [ string | cell array of strings {''} ]
%       These are species full names in the model or just species names
%       that can accept a string in compartment as its compartment.
%   products: [ string | cell array of strings {''}  ]
%       Like reactants.
%   kForward: [ string ]
%       This is a rate parameter name in the model.
%   kReverse: [ string ]
%       Like kForward.
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up
if nargin < 6
    kForward = [];
end

if isempty(reactants)
    reactants = cell(1,0);
end
if isempty(products)
    products = cell(1,0);
end

% Standardize names
[name1, unused] = fixReactionName(name, kForward, '');
for iReac = 1:numel(reactants)
    reactants{iReac} = fixReactionSpecies(reactants{iReac});
end
for iProd = 1:numel(products)
    products{iProd} = fixReactionSpecies(products{iProd});
end
kForward = fixReactionParameter(kForward);
compartment = fixReactionCompartment(compartment);

% Increment counter
nr = m.add.nr + 1;
m.add.nr = nr;
m.add.Reactions = growReactions(m.add.Reactions, nr);

% Add item
m.add.Reactions(nr).Name = name1;
m.add.Reactions(nr).Compartment  = compartment;
m.add.Reactions(nr).Reactants = row(reactants);
m.add.Reactions(nr).Products  = row(products);
m.add.Reactions(nr).Parameter = kForward;

m.Ready = false;
