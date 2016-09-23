function simbio = analytic2simbio(m, opts)
% Convert kroneckerbio analytic model to Matlab SimBiology model
% The kroneckerbio model m should be finalized. This cleans up all the
%   expressions and checks model validity.
%
% Inputs:
%   m [ Model.Analytic struct ]
%       Kroneckerbio analytic model struct
%   opts [ options struct ]
%       Options struct with the following fields:
%       .Verbose [ nonnegative integer {0} ]
%           Control printed information. Higher means more is printed.
%       .OutputsAsStates [ true | {false} ]
%           Whether to export outputs as states. If true, outputs are added as
%           boundary condition states with name "kboutput_<output name>" and
%           repeated assignment rule. This prefixing is necessary because of
%           conflicts with other states.
%
% Outputs:
%   simbio [ (SimBio) Model ]
%       SimBiology Model object
%
% Notes:
% - State initial conditions can be arbitrary math expressions of seeds and
%   rate constants. An initial assignment rule is created for each state
%   whose initial condition isn't a simple constant number.
% - 'null' is a placeholder value for a 0-th order production or a degradation
%   to nothing reaction. 'null' is not allowed as a species name.
%
% Limitations (these should be removed as features are implemented):
% - Compartments must be constant size, with size being (converted into) a
%   single number.
% - The const/default value of inputs are used as their initial amounts
% - kroneckerbio Rules aren't implemented yet. Priority TODO: figure out
%   how they work nowadays.
% - Stuff from experiments aren't in here. So the actual expressions for
%   inputs need to be pulled in later.

%% Clean up inputs
if nargin < 2
    opts = [];
end

opts_ = [];
opts_.Verbose = 0;
opts_.OutputsAsStates = false;
opts = mergestruct(opts_, opts);
verbose = logical(opts.Verbose);

assert(strcmp(m.Type, 'Model.Analytic'), 'analytic2simbio:InvalidModelType', 'Input model must be a Model.Analytic')

% It may be possible to relax this, but this function relies on proper
%   species qualification
assert(m.Ready, 'analytic2simbio:UnfinalizedModel', 'Input model must be finalized') 

if verbose; fprintf('Converting analytic model to SimBio...'); end

%% Initialize SimBio model
if isempty(m.Name)
    name = m.Type;
else
    name = m.Name;
end
simbio = sbiomodel(name);

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
        initialLHS = kroneckerbioExpr2SimbioExpr(name);
        initialRHS = kroneckerbioExpr2SimbioExpr(species.InitialValue);
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
    rate = kroneckerbioExpr2SimbioExpr(reaction.Rate);
    
    nR = length(reactants);
    for iR = 1:nR
        if ~isValidIdentifier(reactants{iR});
            reactants{iR} = kroneckerbioExpr2SimbioExpr(['"' reactants{iR} '"']);
        end
    end
    if nR == 0 % placeholder for 0-th order production
        reactants = {'null'};
    end
    
    nP = length(products);
    for iP = 1:nP
        if ~isValidIdentifier(products{iP});
            products{iP} = kroneckerbioExpr2SimbioExpr(['"' products{iP} '"']);
        end
    end
    if nP == 0 % placeholder for degradation to nothing
        products = {'null'};
    end
    
    % Assemble reactants -> products expression
    reactionStr = [strjoin(reactants, ' + ') ' -> ' strjoin(products, ' = ')];
    r_ = addreaction(simbio, reactionStr);
    
    % Assemble reaction rate
    set(r_, 'ReactionRate', rate);
end

%% Add outputs using repeated assignment rules to dummy states
if opts.OutputsAsStates
    for i = 1:m.ny
        output = m.Outputs(i);
        name = output.Name;
        expr = kroneckerbioExpr2SimbioExpr(output.Expression);
        
        dummyName = ['kboutput_' name];
        dummyName = kroneckerbioExpr2SimbioExpr(dummyName);
        
        s_ = addspecies(simbio, dummyName, 0);
        set(s_, 'BoundaryCondition', true);
        
        rule_ = addrule(simbio, [dummyName ' = ' expr]);
        set(rule_, 'RuleType', 'repeatedAssignment');
    end
end

end

