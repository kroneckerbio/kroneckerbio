function m = addInputAnalytic(m, name, compartment, default, id)
%AddInput Add an input species to a Model.Analytic
%
%   m = AddInput(m, name, compartment, default)
%
%   Inputs
%   m: [ symbolic model struct scalar ]
%       The model to which the input species will be added
%   name: [ string ]
%       A name for the input. This is the name by which reactions will
%       refer to it.
%   compartment: [ string {1st compartment in model} ]
%       The name of the compartment to which it will be added.
%   default: [ nonnegative scalar {0} ]
%       The default value for this input.
%   id: [ string {random UUID} ]
%       A unique, valid variable name
%
%   Outputs
%   m: [ model struct scalar ]
%       The model with the new input added.
if nargin < 5
    id = [];
    if nargin < 4
        default = [];
        if nargin < 3
            compartment = [];
        end
    end
end

% Get compartment names in model
vNames = [{m.Compartments.Name}; {m.add.Compartments.Name}];
vNames = unique(vNames(~cellfun('isempty',vNames)));
if isempty(vNames)
    error('addInputAnalytic: model has no compartments for state to reside in')
end

% Set defaults
if isempty(compartment)
    compartment = vNames{1};
end
if isempty(default)
    default = 0;
end
if isempty(id)
    id = genUID;
end

% Increment counter
nu = m.add.nu + 1;
m.add.nu = nu;
m.add.Inputs = growInputsAnalytic(m.add.Inputs, nu);

% Add item
m.add.Inputs(nu).Name = fixSpeciesName(name);
m.add.Inputs(nu).ID = id;
m.add.Inputs(nu).Compartment = fixCompartmentName(compartment);
m.add.Inputs(nu).DefaultValue = fixInputDefaultValue(default);

m.Ready = false;