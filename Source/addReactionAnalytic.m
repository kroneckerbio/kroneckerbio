function m = addReactionAnalytic(m, name, reactants, products, kForward, kReverse, id)
%AddReaction Add a reaction to a Model.Analytic
%
%   m = AddReaction(m, name, reactants, products, kForward, kReverse, id)
%
%   A reaction is a conversion of one or two reactant species into one or
%   two product species associated with a rate expression. The
%   reverse reaction may or may not be specified at the same time. Each
%   reactant and product must be a single species in a single compartment
%   (e.g. cytoplasm.glucose).
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
%   reactants: [ cell array of strings | string | empty ]
%       This is a species full name in the model or just a species name
%       that can accept a string in compartment as its compartment.
%   products: [ string ]
%       Like reactants
%   kForward: [ string ]
%       Expression for forward reaction rate
%   kReverse: [ string {''} ]
%       Expression for reverse reaction rate
%   id: [ string {[]} | 1 x 2 cell vector of strings ]
%       A unique, valid variable name for the forward and reverse reactions
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new reaction added.

% Clean-up
if nargin < 7
    id = [];
    if nargin < 6
        kReverse = [];
    end
end

% Set defaults
if isempty(id)
    id = '';
end

% Standardize IDs
if isempty(id)
    id1 = '';
    id2 = '';
elseif iscell(id)
    if numel(id) == 1
        assert(isempty(kForward) || isempty(kReverse), 'KroneckerBio:AddReaction:id', 'One id was provided, but forward and reverse reactions both need an id.')
        id1 = id{1};
        id2 = id{1};
    elseif numel(id) >= 2
        id1 = id{1};
        id2 = id{2};
    else
        id1 = [];
        id2 = [];
    end
elseif ischar(id)
    assert(isempty(kForward) || isempty(kReverse), 'KroneckerBio:AddReaction:id', 'One id was provided, but forward and reverse reactions both need an id.')
    id1 = id;
    id2 = id;
end

% Standardize reaction name
[name1, name2] = fixReactionName(name, kForward, kReverse);

% Standardize species names into compartment.species
reactants = fixReactionSpecies(reactants);
products =  fixReactionSpecies(products);

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
    
    if isempty(id1)
        id1 = '';
    elseif issym(id1)
        id1 = char(id1);
    end
    
    m.add.Reactions(nr).ID = id1;
    
    m.add.Reactions(nr).Reactants = reactants;
    
    m.add.Reactions(nr).Products = products;
    
    m.add.Reactions(nr).Rate = kForward;
    
    m.Ready = false;
end

if ~isempty(kReverse)
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactionsAnalytic(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name2;
    
    if isempty(id2)
        id2 = '';
    elseif issym(id2)
        id2 = char(id2);
    end
    m.add.Reactions(nr).ID = id2;
    
    m.add.Reactions(nr).Reactants = products;
    
    m.add.Reactions(nr).Products = reactants;
    
    m.add.Reactions(nr).Rate = kReverse;
    
    m.Ready = false;
end
