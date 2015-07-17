function m = addReactionAnalytic(m, name, reactant1, reactant2, product1, product2, kForward, kReverse, compartment, id)
%AddReaction Add a reaction to a Model.Analytic
%
%   m = AddReaction(m, name, reactant1, reactant2, product1,
%                   product2, kForward, kReverse)
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
if nargin < 10
    id = [];
    if nargin < 9
        compartment = [];
        if nargin < 8
            kReverse = [];
            if nargin < 7
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

% Standardize IDs
if iscell(id)
    if numel(id) == 1
        id1 = id{1};
        id2 = [];
    elseif numel(id) >= 2
        id1 = id{1};
        id2 = id{2};
    else
        id1 = [];
        id2 = [];
    end
elseif ischar(id)
    id1 = id;
    id2 = [];
elseif isempty(id)
    id1 = [];
    id2 = [];
end

% Standardize reaction name
[name1, name2] = fixReactionName(name, kForward, kReverse);

% Standardize species names into compartment.species
reactant1 = fixSpeciesFullName(reactant1, compartment, m);
reactant2 = fixSpeciesFullName(reactant2, compartment, m);
product1  = fixSpeciesFullName(product1, compartment, m);
product2  = fixSpeciesFullName(product2, compartment, m);

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
        id1 = genUID;
    end
    m.add.Reactions(nr).ID = id1;
    
    m.add.Reactions(nr).Reactants = cell(1,2);
    if ~isempty(reactant1)
        m.add.Reactions(nr).Reactants{1} = reactant1;
    end
    if ~isempty(reactant2)
        m.add.Reactions(nr).Reactants{2} = reactant2;
    end
    
    m.add.Reactions(nr).Products = cell(1,2);
    if ~isempty(product1)
        m.add.Reactions(nr).Products{1} = product1;
    end
    if ~isempty(product2)
        m.add.Reactions(nr).Products{2} = product2;
    end
    
    m.add.Reactions(nr).Rate = kForward;
    
    m.Ready = false;
end

if ~isempty(kReverse)
    nr = m.add.nr + 1;
    m.add.nr = nr;
    m.add.Reactions = growReactionsAnalytic(m.add.Reactions, nr);
    
    m.add.Reactions(nr).Name = name2;
    
    if isempty(id2)
        id2 = genUID;
    end
    m.add.Reactions(nr).ID = id2;
    
    m.add.Reactions(nr).Reactants = cell(1,1);
    if ~isempty(product1)
        m.add.Reactions(nr).Reactants{1} = product1;
    end
    if ~isempty(product2)
        m.add.Reactions(nr).Reactants{2} = product2;
    end
    
    m.add.Reactions(nr).Products = cell(1,2);
    if ~isempty(reactant1)
        m.add.Reactions(nr).Products{1} = reactant1;
    end
    if ~isempty(reactant2)
        m.add.Reactions(nr).Products{2} = reactant2;
    end
    
    m.add.Reactions(nr).Rate = kReverse;
    
    m.Ready = false;
end
