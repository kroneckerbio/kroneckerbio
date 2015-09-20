function m = sbml2analytic(sbml)

m = InitializeModelAnalytic(sbml.name);

%% Clean up model
% TODO: ensure that compartment, parameters, etc have unique names
for iv = 1:numel(sbml.compartment)
    vi = sbml.compartment(iv);
    if isempty(vi.name)
        sbml.compartment(iv).name = vi.id;
    end
end
v_ids = vec({sbml.compartment.id});
v_names = vec({sbml.compartment.name});

for ik = 1:numel(sbml.parameter)
    ki = sbml.parameter(ik);
    if isempty(ki.name)
        sbml.parameter(ik).name = ki.id;
    end
end
k_ids = vec({sbml.parameter.id});
k_names = vec({sbml.parameter.name});

for ixu = 1:numel(sbml.species)
    xui = sbml.species(ixu);
    if isempty(xui.name)
        sbml.species(ixu).name = xui.id;
    end
    compartment_index = lookupmember(xui.compartment, v_ids);
    sbml.species(ixu).compartment_index = compartment_index;
    sbml.species(ixu).compartment_name = sbml.compartment(compartment_index).name;
end
nxu = numel(sbml.species);
xu_ids = vec({sbml.species.id});
xu_names = vec({sbml.species.name});
vxu_names = vec({sbml.species.compartment_name});
xu_full_names = strcat(vxu_names, '.', xu_names);

for ir = 1:numel(sbml.reaction)
    ri = sbml.reaction(ir);
    % Reactions are allowed to have empty names
    
    for ik = 1:numel(ri.kineticLaw.parameter)
        ki = ri.kineticLaw.parameter(ik);
        if isempty(ki.name)
            sbml.reaction(ir).kineticLaw.parameter(ik).name = ki.id;
        end
    end
end

% Copy stand-alone initial conditions to states
if isfield(sbml, 'initialAssignment')
    for ii = 1:numel(sbml.initialAssignment)
        i_i = sbml.initialAssignment(ii);
        i_state = lookupmember(i_i.symbol, xu_ids);
        
        assert(i_state ~= 0, 'KroneckerBio:SBML:UnsupportedInitialAssignment', 'Only species can be targets of initial assignment rules')
        
        sbml.species(i_state).isSetInitialAmount = true;
        sbml.species(i_state).initialAmount = i_i.math;
    end
end

z_handled = false(numel(sbml.rule),1); % Stores the rules need to be deleted
for iz = 1:numel(sbml.rule)
    zi = sbml.rule(iz);
    type = zi.typecode;
    name = zi.name;
    metaid = zi.metaid;
    target = zi.variable;
    formula = zi.formula;
    
    if strcmp(type, 'SBML_ASSIGNMENT_RULE')
        % Normal rule
    elseif strcmp(type, 'SBML_INITIAL_RULE')
        % Copy to state initial condition
        i_state = lookupmember(target, xu_ids);
        
        assert(i_state ~= 0, 'KroneckerBio:SBML:UnsupportedInitialAssignment', 'Only species can be targets of initial assignment rules')
        
        sbml.species(i_state).isSetInitialAmount = true;
        sbml.species(i_state).initialAmount = formula;
        
        z_handled(iz) = true;
    elseif strcmp(type, 'SBML_RATE_RULE')
        % Append as a reaction with a single product
        product = struct;
        product.species = target;
        product.stoichiometry = 1;
        
        kinetic_law = struct;
        kinetic_law.formula = formula;

        rate = struct;
        rate.metaid = metaid;
        rate.name = name;
        rate.product = product;
        rate.kineticLaw = kinetic_law;
        
        sbml.reaction = insertstruct(sbml.reaction, rate, numel(sbml.reaction)+1,1);

        z_handled(iz) = true;
    else
        error('KroneckerBio:SBML:UnsupportedRuleType', 'Unsupported rule type %s', type)
    end
end

sbml.rule(z_handled) = [];

%% Extract model
for iv = 1:numel(sbml.compartment)
    vi = sbml.compartment(iv);
    name = vi.name;
    id = vi.id;
    dim = vi.spatialDimensions;
    if vi.isSetSize
        size = vi.size;
    else
        % Had better be set by a rule
        size = NaN;
    end
    
    m = AddCompartment(m, name, dim, size, id);
end

for ik = 1:numel(sbml.parameter)
    ki = sbml.parameter(ik);
    name = ki.name;
    id = ki.id;
    value = ki.value;
    
    m = AddParameter(m, name, value, id);
end

for ixu = 1:nxu
    xi = sbml.species(ixu);
    name = xi.name;
    id = xi.id;
    compartment_index = lookupmember(xi.compartment, v_ids);
    compartment = xi.compartment_name;
    is_input = xi.boundaryCondition || xi.constant;
    
    if xi.isSetInitialAmount
        value = xi.initialAmount;
    elseif xi.isSetInitialConcentration
        % TODO: determine if this is supposed to happen before or after initial assignments of compartment size
        value = xi.initialConcentration / sbml.compartment(compartment_index).size;
    else
        % Had better be set by a rule
        value = NaN;
    end
    
    if is_input
        m = AddInput(m, name, compartment, value, id);
    else
        m = AddState(m, name, compartment, value, id);
    end
end

for ir = 1:numel(sbml.reaction)
    ri = sbml.reaction(ir);
    name = ri.name;
    
    reactants = expand_reactants(ri.reactant, xu_full_names, xu_ids);
    products = expand_reactants(ri.product, xu_full_names, xu_ids);
    
    formula = ri.kineticLaw.formula;
    
    % Add parameters confined to reaction
    for ik = 1:numel(ri.kineticLaw.parameter)
        ki = ri.kineticLaw.parameter(ik);
        name_ik = ki.name;
        id = ki.id;
        value = ki.value;
        
        m = AddParameter(m, name_ik, value, id);
    end
    
    m = AddReaction(m, name, reactants, products, formula);
end

for iz = 1:numel(sbml.rule)
    % Only repeated assignment rules are left
    
    zi = sbml.rule(iz);
    target = zi.variable;
    formula = zi.formula;
    
    % See if it is a compartment size
    compartment_index = lookupmember(target, v_ids);
    if compartment_index ~= 0
        sbml.add.Compartments(compartment_index).Size = formula;
        continue
    end
    
    % See if it is a parameter value
    k_names = vec({m.add.Parameters(1:m.add.nk).Name});
    k_ids = vec({m.add.Parameters(1:m.add.nk).ID});
    parameter_index = lookupmember(target, k_ids);
    if parameter_index ~= 0
        % Extract parameter
        name = k_names{parameter_index};
        id = k_ids{parameter_index};
        
        % Remove parameter
        m = RemoveParameter(m, name);
        
        m = AddRule(m, name, formula, id);
        continue
    end
    
    % See if it is a species value
    xu_full_names = vec({m.add.Inputs(1:m.add.nu).Name, m.add.States(1:m.add.nx).Name});
    xu_ids = vec({m.add.Inputs(1:m.add.nu).ID, m.add.States(1:m.add.nx).ID});
    species_index = lookupmember(target, xu_ids);
    if species_index ~= 0
        % Extract species
        name = xu_full_names{species_index};
        id = xu_ids{species_index};
        
        % Remove species
        if species_index <= m.add.nu
            m = RemoveInput(m, name);
        else
            m = RemoveState(m, name);
        end
        
        m = AddRule(m, name, formula, id);
        continue
    end
end

assert(isempty(sbml.functionDefinition), 'KroneckerBio:SBML:functions', 'Model contains a function definition that is not surrently supported.')

%% Test code to make formulas readable
all_ids = vec({m.add.Compartments(1:m.add.nv).ID, m.add.Inputs(1:m.add.nu).ID, m.add.States(1:m.add.nx).ID, m.add.Parameters(1:m.add.nk).ID, m.add.Rules(1:m.add.nz).ID});
all_names = vec([{m.add.Compartments(1:m.add.nv).Name}, strcat({m.add.Inputs(1:m.add.nu).Compartment}, '.', {m.add.Inputs(1:m.add.nu).Name}), strcat({m.add.States(1:m.add.nx).Compartment}, '.', {m.add.States(1:m.add.nx).Name}), {m.add.Parameters(1:m.add.nk).Name}, {m.add.Rules(1:m.add.nz).Name}]);

for ix = 1:m.add.nx
    m.add.States(ix).InitialValue = substituteQuotedExpressions(m.add.States(ix).InitialValue, all_ids, all_names, true);
end
for ir = 1:m.add.nr
    m.add.Reactions(ir).Rate = substituteQuotedExpressions(m.add.Reactions(ir).Rate, all_ids, all_names, true);
end
for iz = 1:m.add.nz
    m.add.Rules(iz).Expression = substituteQuotedExpressions(m.add.Rules(iz).Expression, all_ids, all_names, true);
end

end

function reactant_names = expand_reactants(reactants, xu_full_names, xu_ids)
n = numel(reactants);
reactant_names = cell(1,0);
for i = 1:n
    assert(isempty(reactants(i).stoichiometryMath), 'KroneckerBio:SBML:StoichiometryMath', 'Use of stiochiomerty math is not supported')

    id = reactants(i).species;
    name = xu_full_names{lookupmember(id, xu_ids)};
    stoich = reactants(i).stoichiometry;

    reactant_names = [reactant_names, repmat({name}, 1,stoich)];
end
end
