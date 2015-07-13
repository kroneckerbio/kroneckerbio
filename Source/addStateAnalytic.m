function m = addStateAnalytic(m, name, compartment, seed, id)
%AddState Add a state species to a KroneckerBio model. If the state's initial
%condition depends on a seed, the seed must already be added to the model.
%
%   m = AddState(m, name, compartment, seed)
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the state will be added
%   name: [ string ]
%       A name for the state. This is the name by which reactions will
%       refer to it.
%   compartment: [ string {1st compartment in model} ]
%       The name of the compartment to which it will be added
%   seed: [ string | nonnegative scalar {0} ]
%       String expression for the initial condition in terms of seeds or
%       nonnegative scalar initial amount.
%   id: [ string {random UUID} ]
%       A unique, valid variable name
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new state added.
%

% Clean up inputs
if nargin < 5
    id = [];
    if nargin < 4
        seed = [];
        if nargin < 3
            compartment = [];
        end
    end
end

% Get compartment names in model
vNames = [{m.Compartments.Name}, {m.add.Compartments.Name}];
vNames = unique(vNames(~cellfun('isempty',vNames)));
if isempty(vNames)
    error('addStateAnalytic: model has no compartments for state to reside in')
end

% Set defaults
if isempty(compartment)
    compartment = vNames{1};
end
if isempty(seed)
    seed = 0;
end
if isempty(id)
    id = genUID;
end

% Clean up initial condition
if ischar(seed)
    % pass
elseif isnumeric(seed) && isscalar(seed)
    seed = num2str(seed);
else
    error('addStateAnalytic: seed type not recognized')
end

% Increment counter
nx = m.add.nx + 1;
m.add.nx = nx;
m.add.States = growStatesAnalytic(m.add.States, nx);

% Add item
m.add.States(nx).Name = fixSpeciesName(name);
m.add.States(nx).ID = id;
m.add.States(nx).Compartment = fixCompartmentName(compartment);
m.add.States(nx).InitialValue = seed;

m.Ready = false;
