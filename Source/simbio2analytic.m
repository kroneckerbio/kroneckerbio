function m = simbio2analytic(simbio)
% Convert Matlab SimBiology model to kroneckerbio analytic model

m = InitializeModelAnalytic(simbio.name);

%% Extract components
nv = numel(simbio.Compartments);
v_names = cell(nv,1);
v_sizes = cell(nv,1);
for iv = 1:nv
    vi = simbio.Compartments(iv);
    
    v_names{iv} = vi.Name;
    v_sizes{iv} = vi.Capacity;
end

nk = numel(simbio.Parameters);
k_names = cell(nk,1);
k_values = zeros(nk,1);
for ik = 1:nk
    ki = simbio.Parameters(ik);
    
    k_names{ik} = ki.Name;
    k_values(ik) = ki.Value;
end

nxu = numel(simbio.Species);
xu_names = cell(nxu,1);
vxu_names = cell(nxu,1);
xu_is_inputs = false(nxu,1);
xu_values = cell(nxu,1);
for ixu = 1:nxu
    xui = simbio.Species(ixu);
    
    xu_names{ixu} = xui.Name;
    vxu_names{ixu} = xui.Parent.Name;
    xu_is_inputs(ixu) = xui.BoundaryCondition || xui.ConstantAmount;
    xu_values{ixu} = xui.InitialAmount;
end
xu_full_names = strcat(vxu_names, '.', xu_names);

nz = numel(simbio.rule);
z_names = cell(nz,1);
z_values = cell(nz,1);
z_types = cell(nz,1);
for iz = 1:nz
    zi = simbio.Rules(iz);
    
    splits = regexp(zi.Rule, '=', 'split'); % TODO: make split smarter
    assert(numel(splits) == 2, 'KroneckerBio:simbio2analytic:InvalidRule', 'Rule #%i (%s) had an unparsible expression', iz, zi.Name)
    z_names{iz} = strtrim(splits{1});
    z_values{iz} = clean_simbio_expression(strtrim(splits{2}));
    
    z_types{iz} = zi.RuleType;
end

%% Extract reactions
nr = numel(simbio.reaction);
r_names = cell(nr,1);
r_reactants = cell(nr,1);
r_products = cell(nr,1);
r_rates = cell(nr,1);
for ir = 1:nr
    ri = simbio.Reactions(ir);
    
    r_names{ir} = ri.Name;
    
    n_reac = numel(ri.Reactants);
    r_reactants{ir} = cell(n_reac,2);
    for i_reac = 1:n_reac
        r_reactants{ir}(i_reac,:) = {[ri.Reactants(i_reac).Parent.Name '.' ri.Reactants(i_reac).Name], -ri.Stoichiometry(i_reac)}; % Simbio stoichiometry is negative for reactants
    end
    
    n_prod = numel(ri.Products);
    r_products{ir} = cell(n_prod,2);
    for i_prod = 1:n_prod
        r_products{ir}(i_prod,:) = {[ri.Products(i_prod).Parent.Name '.' ri.Products(i_prod).Name], ri.Stoichiometry(n_reac+i_prod)};
    end
    
    r_rates{ir} = clean_simbio_expression(ri.ReactionRate);
    
    kineticLaw = ri.kineticLaw; % Will be empty if no kinetic law parameters exist
    if ~isempty(kineticLaw)
        ki = kineticLaw.Parameters; % Will only fetch parameters unique to this kinetic law

        nkl = length(ki);
        for j = 1:nkl
            parameter = ki(j);
            
            nk = nk + 1;
            k_names = [k_names; parameter.Name];
            k_values = [k_values; parameter.Value];
        end
    end
end

%% Process rules
z_handled = false(nz,1);

% Copy all initial assignment rules to initial conditions of states
z_is_initial = strcmp('initialAssignment', z_types);
xu_getting_replaced = lookupmember(z_names(z_is_initial), xu_full_names); % All simbio species appear as full names
assert(all(xu_getting_replaced ~= 0), 'KroneckerBio:SBML:UnsupportedInitialAssignment', 'Only species can be targets of initial assignment rules')
xu_values(xu_getting_replaced) = z_values(z_is_initial);
z_handled(z_is_initial) = true;

% Append rate rules are reactions
z_is_rate = strcmp('rate', z_types);
r_names = [r_names; z_names(z_is_rate)];
r_reactants = [r_reactants; repmat({cell(0,2)}, nnz(z_is_rate),1)];
r_products = [r_products; cellfun(@(id){{id,1}}, z_names(z_is_rate))];
r_rates = [r_rates; z_values(z_is_rate)];
z_handled(z_is_rate) = true;

% Assume that all rules left are assignment rules

% Compartment assignments get copied into compartment size
[z_is_compartment, v_getting_replaced] = ismember(z_names, v_names);
v_getting_replaced(v_getting_replaced == 0) = [];
v_sizes(v_getting_replaced) = z_values(z_is_compartment);
z_handled(z_is_compartment) = true;

% Parameter assignments delete the parameter
[z_is_parameter, k_handled] = ismember(z_names, k_names);
k_handled(k_handled == 0) = [];

% Species assignments delete the species
[z_is_species, xu_handled] = ismember(z_names, xu_full_names);
z_is_species = z_is_species & ~z_handled; % Don't do anything with species already handled as initial conditions
xu_handled(xu_handled == 0 | z_handled) = [];
z_names(z_is_species) = xu_names(xu_handled);

% Species assignments also require removing the compartment from every
% expression in which it appears
z_values = substituteQuotedExpressions(z_values, xu_full_names(xu_handled), xu_names(xu_handled), true);
r_rates = substituteQuotedExpressions(r_rates, xu_full_names(xu_handled), xu_names(xu_handled), true);

% Purge handled objects
nz = nz - nnz(z_handled);
z_names(z_handled) = [];
z_values(z_handled) = [];

nk = nk - nnz(k_handled);
k_names(k_handled) = [];
k_values(k_handled) = [];

nxu = nxu - nnz(xu_handled);
xu_names(xu_handled) = [];
vxu_names(xu_handled) = [];
xu_is_inputs(xu_handled) = [];
xu_values(xu_handled) = [];

%% Add components to analytic model
for iv = 1:nv
    % Simbiology does not expose the dimension of the compartment
    m = AddCompartment(m, v_names{iv}, 3, v_sizes{iv});
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

function expressions = clean_simbio_expression(expressions)
% Assume Simbiology models cannot contain brackets in their names. These
% are special characters. Find all identifiers that are exclosed between
% brackets and replace the brackets with quotes.

temp = @transmute; % Matlab bug can't find local functions
expressions = regexprep(expressions, '(\[[^\[\]]*\])', '${temp($1)}');

    function output = transmute(input)
        if input(1) == '['
            % Quoted branch
            output = ['"' input(2:end-1) '"'];
        else
            % Dot branch
            output = ['"' input '"'];
        end
    end
end
