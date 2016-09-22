function simbio = analytic2simbio(m)
% Convert kroneckerbio analytic model to Matlab SimBiology model
% The kroneckerbio model m should be finalized. This cleans up all the
%   expressions and checks model validity.
%
% Inputs:
%   m [ Model.Analytic struct ]
%       Kroneckerbio analytic model struct
%
% Outputs:
%   simbio [ (SimBio) Model ]
%       SimBiology Model object
%
% Notes:
% - State initial conditions can be arbitrary math expressions of seeds and
%   rate constants. An initial assignment rule is created for each state
%   whose initial condition isn't a simple constant number.
% - Outputs are added as dummy species with name "out_<output name>" as a
%   boundary condition with repeated assignment rule as its expression. We
%   can add an options struct to this function's input args to control this
%   behavior if desired.
%
% Limitations (these should be removed as features are implemented):
% - Compartments must be constant size, with size being (converted into) a
%   single number.
% - The const/default value of inputs are used as their initial amounts
% - kroneckerbio Rules aren't implemented yet. Priority TODO: figure out
%   how they work nowadays.
% - Stuff from experiments aren't in here. So the actual expressions for
%   inputs need to be pulled in later.

assert(strcmp(m.Type, 'Model.Analytic'), 'analytic2simbio:InvalidModelType', 'Input model must be a Model.Analytic')

% It may be possible to relax this, but this function relies on proper
%   species qualification
assert(m.Ready, 'analytic2simbio:UnfinalizedModel', 'Input model must be finalized') 

%% Initialize SimBio model
simbio = sbiomodel([m.Name m.Type]);

%% Add compartments
vNames = {m.Compartments.Name};
for i = 1:m.nv
    compartment = m.Compartments(i);
    name = compartment.Name;
    size = compartment.Size;
    if ischar(size)
        size = str2double(size);
    end
    
    addcompartment(simbio, name, size); % ConstantCapacity = true by default
end

%% Add parameters
% All parameters are entire-reaction-scoped
% Rate parameters are regular parameters
for i = 1:m.nk
    parameter = m.Parameters(i);
    name = parameter.Name;
    value = parameter.Value;
    
    addparameter(simbio, name, value);
end

% Seeds are also treated as regular parameters
for i = 1:m.ns
    parameter = m.Seeds(i);
    name = parameter.Name;
    value = parameter.Value;
    
    addparameter(simbio, name, value);
end

%% Add species
% States are modified by reactions
for i = 1:m.nx
    species = m.States(i);
    name = species.Name;
    compartment = species.Compartment;
    initial = str2double(species.InitialValue);
    initialConst = true;
    if isnan(initial) % didn't parse to a single constant number
        initial = 0;
        initialLHS = sprintf(name, i);
        initialRHS = species.InitialValue;
        initialConst = false;
    end
    
    vInd = find(ismember(vNames, compartment));
    addspecies(simbio.Compartments(vInd), name, initial);
    
    if ~initialConst
        rule_ = addrule(simbio, [initialLHS ' = ' initialRHS]);
        set(rule_, 'RuleType', 'initialAssignment');
    end
end

% Inputs aren't modified by reactions
% The const/default value of inputs are used as their initial amounts
% SimBio BoundaryCondition aren't modified by reactions but are modified
%   by rules, events, and doses. SimBio ConstantAmount aren't modified by
%   anything.
% See https://www.mathworks.com/help/simbio/ref/boundarycondition.html
for i = 1:m.nu
    species = m.Inputs(i);
    name = species.Name;
    compartment = species.Compartment;
    initial = species.DefaultValue;
    
    vInd = find(ismember(vNames, compartment));
    s_ = addspecies(simbio.Compartments(vInd), name, initial);
    set(s_, 'BoundaryCondition', true);
end

%% Add reactions
for i = 1:m.nr
    reaction = m.Reactions(i);
    reactants = reaction.Reactants;
    products = reaction.Products;
    rate = reaction.Rate;
    
    % Assemble reactants -> products expression
    reactionStr = [strjoin(reactants, ' + ') ' -> ' strjoin(products, ' = ')];
    r_ = addreaction(simbio, reactionStr);
    
    % Assemble reaction rate
    set(r_, 'ReactionRate', rate);
end

%% Add outputs using repeated assignment rules to dummy states
for i = 1:m.ny
    output = m.Outputs(i);
    name = output.Name;
    expr = output.Expression;
    
    dummyName = ['out_' name];
    
    s_ = addspecies(simbio, dummyName, 0);
    set(s_, 'BoundaryCondition', true);
    
    rule_ = addrule(simbio, [dummyName ' = ' expr]);
    set(rule_, 'RuleType', 'repeatedAssignment');
end

end

