function m = sbml2analytic(sbml)

m = InitializeModelAnalytic(sbml.name);

%% Extract components
nv = numel(sbml.compartment);
v_ids = vec({sbml.compartment.id});
v_names = vec({sbml.compartment.name});
v_dims = double(vec([sbml.compartment.spatialDimensions])); % Sometimes Matlab freaks out if this is not a double

nk = numel(sbml.parameter);
k_ids = vec({sbml.parameter.id});
k_names = vec({sbml.parameter.name});

nxu = numel(sbml.species);
xu_ids = vec({sbml.species.id});
xu_names = vec({sbml.species.name});
vxu_ids = vec({sbml.species.compartment});
vxu_indexes = lookupmember(vxu_ids, v_ids);
xu_is_inputs = vec([sbml.species.boundaryCondition]) | vec([sbml.species.constant]);

nz = numel(sbml.rule);
z_ids = vec({sbml.rule.variable});
z_names = vec({sbml.rule.name});
z_values = vec({sbml.rule.formula});
z_types = vec({sbml.rule.typecode});

%% Deal with missing component values
% Default to NaN, rules can override this
v_set_size = logical(vec([sbml.compartment.isSetSize]));
v_sizes = repmat({nan}, nv,1);
v_sizes(v_set_size) = vec({sbml.compartment(v_set_size).size});

k_set_size = logical(vec([sbml.parameter.isSetValue]));
k_values = nan(nk,1);
k_values(k_set_size) = vec([sbml.parameter(k_set_size).value]);

% Species are super weird because they can have amounts or concentrations
% It is not clear how to deal with this, but this converts concentrations
% to amounts using the pre-rule-applied compartment sizes
xu_values = repmat({nan}, nv,1);
xu_set_amount = logical(vec([sbml.species.isSetInitialAmount]));
xu_set_conc = logical(vec([sbml.species.isSetInitialConcentration]));

xu_raw_amounts = vec([sbml.species.initialAmount]);
xu_raw_concentrations = vec([sbml.species.initialConcentration]);
xu_raw_compartment_sizes = vec([v_sizes{vxu_indexes}]);

xu_values(xu_set_conc) = num2cell(xu_raw_concentrations(xu_set_conc) ./ xu_raw_compartment_sizes(xu_set_conc));
xu_values(xu_set_amount) = num2cell(xu_raw_amounts(xu_set_amount));

%% Extract reactions
nr = numel(sbml.reaction);
r_names = vec({sbml.reaction.name});
r_reactants = arrayfun(@(ri){[vec({ri.reactant.species}), vec({ri.reactant.stoichiometry})]}, sbml.reaction);
r_products = arrayfun(@(ri){[vec({ri.product.species}), vec({ri.product.stoichiometry})]}, sbml.reaction);
rates = [sbml.reaction.kineticLaw];
r_rates = vec({rates.formula});

% Seperate out reaction parameters to be combined with regular parameters
reaction_parameters = vec([rates.parameter]);
nrk = numel(reaction_parameters);
rk_ids = vec({reaction_parameters.id});
rk_names = vec({reaction_parameters.name});

% Do same handling of unset parameters
rk_set_size = logical(vec([reaction_parameters.isSetValue]));
rk_values = nan(nrk,1);
rk_values(rk_set_size) = vec([reaction_parameters(rk_set_size).value]);

%% Append reaction parameters to regular parameters
nk = nk + nrk;
k_ids = [k_ids; rk_ids];
k_names = [k_names; rk_names];
k_values = [k_values; rk_values];

%% Copy ids to names if names are missing
v_names_missing = cellfun(@isempty, v_names);
v_names(v_names_missing) = v_ids(v_names_missing);

k_names_missing = cellfun(@isempty, k_names);
k_names(k_names_missing) = k_ids(k_names_missing);

xu_names_missing = cellfun(@isempty, xu_names);
xu_names(xu_names_missing) = xu_ids(xu_names_missing);

%% Collect species compartments and full names
vxu_names = v_names(vxu_indexes);
xu_full_names = strcat(vxu_names, '.', xu_names);

%% Copy stand-alone initial conditions to states
if isfield(sbml, 'initialAssignment')
    initials_ids = vec({sbml.initialAssignment.symbol});
    initials_values = vec({sbml.initialAssignment.math});
    
    initials_indexes = lookupmember(initials_ids, xu_ids);
    assert(all(initials_indexes ~= 0), 'KroneckerBio:SBML:UnsupportedInitialAssignment', 'Only species can be targets of initial assignment rules')

    xu_values(initials_indexes) = initials_values;
end

%% Process rules
z_handled = false(nz,1);

% Copy all initial assignment rules to initial conditions of states
z_is_initial = strcmp('SBML_INITIAL_RULE', z_types);
xu_getting_replaced = lookupmember(z_ids(z_is_initial), xu_ids);
assert(all(xu_getting_replaced ~= 0), 'KroneckerBio:SBML:UnsupportedInitialAssignment', 'Only species can be targets of initial assignment rules')
xu_values(xu_getting_replaced) = z_values(z_is_initial);
z_handled(z_is_initial) = true;

% Append rate rules are reactions
z_is_rate = strcmp('SBML_RATE_RULE', z_types);
nr = nr + nnz(z_is_rate);
r_names = [r_names; z_names(z_is_rate)];
r_reactants = [r_reactants; repmat({cell(0,2)}, nnz(z_is_rate),1)];
r_products = [r_products; cellfun(@(id){{id,1}}, z_ids(z_is_rate))];
r_rates = [r_rates; z_values(z_is_rate)];
z_handled(z_is_rate) = true;

% Assume that all rules left are assignment rules

% Compartment assignments get copied into compartment size
[z_is_compartment, v_getting_replaced] = ismember(z_ids, v_ids);
v_getting_replaced(v_getting_replaced == 0) = [];
v_sizes(v_getting_replaced) = z_values(z_is_compartment);
z_handled(z_is_compartment) = true;

% Parameter assignments delete the parameter
[z_is_parameter, k_handled] = ismember(z_ids, k_ids);
k_handled(k_handled == 0) = [];
z_names(z_is_parameter) = k_names(k_handled);

% Species assignments delete the species
[z_is_species, xu_handled] = ismember(z_ids, xu_ids);
z_is_species = z_is_species & ~z_handled; % Don't do anything with species already handled as initial conditions
xu_handled(xu_handled == 0 | z_handled) = [];
z_names(z_is_species) = xu_names(xu_handled);

% Purge handled objects
nz = nz - nnz(z_handled);
z_ids(z_handled) = [];
z_names(z_handled) = [];
z_values(z_handled) = [];

nk = nk - nnz(k_handled);
k_ids(k_handled) = [];
k_names(k_handled) = [];
k_values(k_handled) = [];

nxu = nxu - nnz(xu_handled);
xu_ids(xu_handled) = [];
xu_names(xu_handled) = [];
vxu_indexes(xu_handled) = [];
vxu_names(xu_handled) = [];
xu_full_names(xu_handled) = [];
xu_is_inputs(xu_handled) = [];
xu_values(xu_handled) = [];

%% Replace all IDs with names
all_ids = vec([v_ids; xu_ids; k_ids; z_ids]);
all_names = vec([v_names; xu_full_names; k_names; z_names]);

v_is_expression = cellfun(@ischar, v_sizes);
xu_is_expression = cellfun(@ischar, xu_values);
r_is_expression =  cellfun(@ischar, r_rates);

v_sizes(v_is_expression) =  substituteQuotedExpressions(v_sizes(v_is_expression), all_ids, all_names, true);
xu_values(xu_is_expression) = substituteQuotedExpressions(xu_values(xu_is_expression), all_ids, all_names, true);
r_rates(r_is_expression) = substituteQuotedExpressions(r_rates(r_is_expression), all_ids, all_names, true);
z_values = substituteQuotedExpressions(z_values, all_ids, all_names, true);

r_reactants = cellfun(@(ri){[xu_full_names(lookupmember(ri(:,1),xu_ids)), ri(:,2)]}, r_reactants);
r_products = cellfun(@(ri){[xu_full_names(lookupmember(ri(:,1),xu_ids)), ri(:,2)]}, r_products);

%% Add components to analytic model
for iv = 1:nv
    m = AddCompartment(m, v_names{iv}, v_dims(iv), v_sizes{iv});
end

for ik = 1:nk
    m = AddParameter(m, k_names{ik}, k_values(ik));
end

for ixu = 1:nxu
    if xu_is_inputs(ixu)
        m = AddInput(m, xu_names{ixu}, vxu_names{ixu}, xu_values{ixu});
    else
        m = AddState(m, xu_names{ixu}, vxu_names{ixu}, xu_values{ixu});
    end
end

for ir = 1:nr
    reactants_i = expand_reactants(r_reactants{ir});
    products_i = expand_reactants(r_products{ir});
    
    m = AddReaction(m, r_names{ir}, reactants_i, products_i, r_rates{ir});
end

for iz = 1:nz
    m = AddRule(m, z_names{iz}, z_values{iz});
end

assert(isempty(sbml.functionDefinition), 'KroneckerBio:SBML:functions', 'Model contains a function definition that is not surrently supported.')

end

function reactant_names = expand_reactants(reactants)
n = size(reactants,1);
reactant_names = cell(1,0);
for i = 1:n
    name = reactants{i,1};
    stoich = reactants{i,2};

    reactant_names = [reactant_names, repmat({name}, 1,stoich)];
end
end
