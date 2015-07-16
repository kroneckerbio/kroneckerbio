function m = addReactionManySpeciesAnalytic(m, name, reactants, products, kForward, kReverse, compartment, id)
%AddReaction Add a reaction with reactants and products specified as a cell
%array to a Model.Analytic. Allows > 2 reactants and product stoichiometries.
%
%   m = AddReaction(m, name, reactants, products, kForward, kReverse)
%
%   A reaction is a conversion of >= 1 reactant species into >= 2
%   product species associated with a rate expression. The
%   reverse reaction may or may not be specified at the same time. Each
%   reactant and product must be a single species in a single compartment
%   (e.g. cytoplasm.glucose).
%
%   Specify stoichiometry > 1 by repeating the species in the reactants or
%   products list.
%
%   Since most reactions happen in a single compartment, Kronecker does not
%   make you type out the compartment for every species. You can specify a
%   compartment in which to search for species without prepended
%   compartments and which appear in multiple compartments. A species that
%   appears in only 1 compartment will ignore the default compartment and be
%   attached to the containing compartment. If the species is not found in the
%   specified compartment or the default compartment (1st compartment in the
%   model if not specified), an error is thrown.
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
%   reactants: [ string | cell array of strings | empty matrix [] ]
%       Species name or list of species names. Specified as "species" if 
%       species is unique or in default compartment or "compartment.species"
%       if species appears in multiple compartments
%   products: [ string ]
%       Like reactants
%   kForward: [ string ]
%       Expression for forward reaction rate
%   kReverse: [ string {''} ]
%       Expression for reverse reaction rate
%   compartment: [ string {'1st compartment of reaction'} ]
%       Default compartment for species w/o a specified compartment if species
%       isn't unique to another compartment
%   id: [ string {random UUID} | 1 x 2 cell vector of strings ]
%       A unique, valid variable name for the forward and reverse reactions
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% Clean-up
if nargin < 8
    id = [];
    if nargin < 7
        compartment = [];
        if nargin < 6
            kReverse = [];
            if nargin < 5
                kForward = [];
            end
        end
    end
end

% Get compartment names in model
vNames = [{m.Compartments.Name}, {m.add.Compartments.Name}];
vNames = unique(vNames(~cellfun('isempty',vNames)));
if isempty(vNames)
    error('addInputAnalytic: model has no compartments for state to reside in')
end

% Set defaults
if isempty(compartment)
    compartment = vNames{1};
end

% Standardize reactants and products into cell array of strings
if isempty(reactants) % 0th-order reaction
    reactants = {[]};
elseif ischar(reactants) % single reactant as a string
    reactants = {reactants};
end
if isempty(products) % degradation reaction
    products = {[]};
elseif ischar(products) % single products as a string
    products = {products};
end

% Standardize IDs
if iscell(id)
    switch numel(id)
        case 0
            id1 = genUID;
            id2 = genUID;
        case 1
            id1 = id{1};
            id2 = genUID;
        case 2
            id1 = id{1};
            id2 = id{2};
    end
elseif ischar(id)
    id1 = id;
    id2 = genUID;
elseif isempty(id)
    id1 = genUID;
    id2 = genUID;
end

% Standardize reaction name
[name1, name2] = fixReactionName(name, kForward, kReverse);

% Standardize species names into compartment.species
for i = 1:length(reactants)
    reactants{i} = fixSpeciesFullName(reactants{i}, compartment, m);
end
for i = 1:length(products)
    products{i} = fixSpeciesFullName(products{i}, compartment, m);
end

% Standardize reaction rate expressions
kForward    = fixRateExpressionAnalytic(kForward);
kReverse    = fixRateExpressionAnalytic(kReverse);

% Add separate reactions for forward and reverse (if applicable)
% Don't add reaction if rate is 0
if ~isempty(kForward)
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactionsAnalytic(m.add.Reactions, nr);
    m.add.Reactions(nr).Name = name1;
    m.add.Reactions(nr).ID = id1;
    m.add.Reactions(nr).Reactants = reactants;
    m.add.Reactions(nr).Products = products;
    m.add.Reactions(nr).Rate = kForward;
end

if ~isempty(kReverse)
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactionsAnalytic(m.add.Reactions, nr);
    m.add.Reactions(nr).Name = name2;
    m.add.Reactions(nr).ID = id2;
    m.add.Reactions(nr).Reactants = products;
    m.add.Reactions(nr).Products = reactants;
    m.add.Reactions(nr).Rate = kReverse;
end

m.Ready = false;
end
