function SymModel = simbio2Symbolic(SimbioModel, opts)
%SSIMBIO2SYMBOLIC converts a Simbiology model into a structure of symbolic
%   variables so that differentiation and other mathematical manipulation
%   can be done on the system.
% 
%   SymModel = simbio2Symbolic(SimbioModel, opts)
%
%   Inputs
%       SimbioModel - A Simbiology model object
%       opts        - Optional function options
%           .Verbose [ true | {false} ]
%               Print progress to command window
%           .Order [ 0 | 1 | {2} | 3 ]
%               Determines how deep the derivatives should be taken with
%               respect to x and p. Each level increases the cost
%               exponentially, but may be required for certain
%               applications.
%
%   Outputs
%       SymModel
%       	.xuNames   - nx by 1 cell vector of strings containing the
%                        natural names of the species
%           .kNames    - nk by 1 cell vector of strings containing the
%                        natural names of the rate constants
%           .xu0        - nx by 1 vector defining the initial conditions
%           .k         - nk by 1 vector defining the values of the rate
%                        constants
%
%           .xuSyms     - nx by 1 vector of symbolics containing new
%                        differentiable names of the species
%           .kSyms     - nk by 1 vector of symbolics containing new
%                        differentiable names of the rate constants
%           xuSyms and kSyms are used in the symbolics below
%
%           .rSyms     - nr by 1 vector of symbolic expressions
%                        representing the rates of each reaction
%           .S         - nx by nr stoichiometry matrix representing how
%                        each reaction changes each species
%           .f         - nx by 1 vector of symbolics containing expressions
%                        for the odes defining the system
%           .dfdx      - nx by nx matrix of symbolics containing
%                        expressions for the partial derivatives on each
%                        ode (rows) with respect to each species (columns)
%           .dfdk      - nx by nk matrix of symbolics containing
%                        expressions for the partial derivatives on each
%                        ode (rows) with respect to each rate constant
%                        (columns)
%           .df2dx2    - nx*nx by nx matrix of symbolics containing
%                        expressions for the double partial derivative on
%                        vectorized dfdx (rows) with respect to the species
%                        (columns)
%           .df2dk2    - nx*nk by nk matrix of symbolics containing
%                        expressions for the double partial derivative on
%                        vectorized dfdk (rows) with respect to the rate
%                        constants (columns)
%           .df2dkdx   - nx*nx by nk matrix of symbolics containing
%                        expressions for the double partial derivative on
%                        vectorized dfdx (rows) with respect to the species
%                        (columns)
%           .df2dxdk   - nx*nk by nx matrix of symbolics containing
%                        expressions for the double partial derivative on
%                        vectorized dfdx (rows) with respect to the rate
%                        constants (columns)
%
%   Limitations:
%   Not all Simbiology features are compatible with this kind of analysis.
%   This function ignores any events, rules, and functions of the model.
%   Reaction.Reversible is ignored
%   Reaction.Active is ignored
%   Unit conversion is completely ignored.
%   All species are assumed to be in dimensions given by DefaultSpeciesDimension.
%   All reactions are assumed to be in dimensions of amount/time.
%   Warning: there may be other features missing as well. Do careful
%   comparisons to simbiology before trusting this code.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Options
% Resolve missing inputs
if nargin < 2
    opts = [];
end

% Copy model object
%SimbioModel = copyobj(SimbioModel);

%Options for displaying progress
defaultOpts.Verbose = 0;
defaultOpts.Order   = 2;

opts = mergestruct(defaultOpts, opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 1: Extracting the Model Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.Verbose; fprintf('Extracting model components...'); end
%% Model name
Name = SimbioModel.Name;

%% Build up the table of compartments
compartments = SimbioModel.Compartments;

nv = length(compartments);
vNames  = cell(nv,1);
vValues = zeros(nv,1);

%Get compartment names
for iv = 1:nv
    vNames{iv}  = compartments(iv).Name;
    vValues(iv) = compartments(iv).Capacity;
end

%% Build up the table of species
species = SimbioModel.Species;

nx = length(species);
xuNames  = cell(nx,1);
xu0      = zeros(nx,1);
xCap     = zeros(nx,1);
vxuNames = cell(nx,1);
isu      = false(nx,1);

%Get species names and compartment volume
for ix = 1:nx %Foreach species
    xuNames{ix}  = species(ix).Name;
    xu0(ix)      = species(ix).InitialAmount;
    xCap(ix)     = species(ix).Parent.Capacity;
    vxuNames{ix} = species(ix).Parent.Name;
    isu(ix)      = species(ix).BoundaryCondition || species(ix).ConstantAmount;
end

%% Build up the table of constants
constants = SimbioModel.Parameters;

nk = length(constants); %number of total parameters
nkm = nk; %length of model parameters
kNames = cell(nk,1);
k      = zeros(nk,1);

%Get model parameters
for ik = 1:nk %Foreach parameter
    kNames{ik} = constants(ik).Name;
    k(ik)      = constants(ik).Value;
end

%Get kinetic law parameters
reactions = SimbioModel.Reactions;

nr = length(reactions);

for ir = 1:nr
    kineticlaw = reactions(ir).KineticLaw; %will be empty if no kinetic law parameters exist
    if ~isempty(kineticlaw)
        constantskl = kineticlaw.Parameters; %will only fetch parameters unique to this kinetic law
        nkkl = length(constantskl);
        for j = 1:nkkl
            nk = nk + 1; %add one more parameter
            kNames{nk,1} = constantskl(j).Name;
            k(nk,1) = constantskl(j).Value;
        end
    end
end

if opts.Verbose; fprintf('done.\n'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 2: Renaming Everything to Allow for Symbolic Handling %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.Verbose; fprintf('Renaming variables...'); end
%% Convert compartment names to differentiable names
%The names to be used will be "co1x", "co2x" etc...
vSyms = sym(zeros(nv,1));
vNicestrs = cell(nv,1);

%We must be careful that the names we assign the variables are not already
%names being used by the model. This code checks if a potential name (i.e.
%"co1x" already exists. If it does, the code simply skips that number and
%tries the next one until it finds an unused name. Note that the "x" on the
%end of the variable names prevents string replace from seeing "sp1"
%in "sp10".

CurrentAttempt = 1; %The numerical appendix that will be tried
CurrentCompartment = 1; %Index through the compartments
AttemptIsGood = true; %becomes true when a name already exists

while(CurrentCompartment <= nv)%Until every compartment
    CurrentName = sprintf('co%dx', CurrentAttempt);
    %Check to make sure a compartment does not have our systematic name
    for check = 1:nv
        %compartments except the one being renamed must not have 'compartment#' as name
        if strcmp(CurrentName, vNames{check}) && CurrentCompartment ~= check
            AttemptIsGood = false;
        end
    end
    %Also check species names (you never know)
    for check = 1:nx
        if strcmp(CurrentName, xuNames{check})
            AttemptIsGood = false;
        end
    end
    %Also check parameter names
    for check = 1:nk
        if strcmp(CurrentName, kNames{check})
            AttemptIsGood = false;
        end
    end
    
    if AttemptIsGood %change it and move on
        SimbioModel.Compartments(CurrentCompartment).rename(CurrentName);
        vSyms(CurrentCompartment) = sym(CurrentName);
        vNicestrs{CurrentCompartment} = CurrentName;
        CurrentCompartment = CurrentCompartment + 1;
        CurrentAttempt = CurrentAttempt + 1;
    else %it failed; try again with different number
        CurrentAttempt = CurrentAttempt + 1;
        AttemptIsGood = true;
    end
end

%% Convert species names to differentiable names
%The names will be "sp1x", "sp2x", etc...
xuSyms = sym(zeros(nx,1));
xNicestrs = cell(nx,1);

%The important part here is to give every species a unique name, regardless
%of what compartment it is in. This way the compartments can be deleted
%without consequence.

CurrentAttempt = 1; %The numerical appendix that will be tried
CurrentSpecies = 1; %Index through the species
AttemptIsGood = true; %becomes true when a name already exists

while(CurrentSpecies <= nx) %Until every species
    CurrentName = sprintf('sp%dx', CurrentAttempt);
    %Check to make sure a compartment does not have our systematic name
    for check = 1:nv
        if strcmp(CurrentName, vNames{check})
            AttemptIsGood = false;
        end
    end
    %Check species names
    for check = 1:nx
        if strcmp(CurrentName, xuNames{check}) && CurrentSpecies ~= check
            AttemptIsGood = false;
        end
    end
    %Also check parameter names
    for check = 1:nk
        if strcmp(CurrentName, kNames{check})
            AttemptIsGood = false;
        end
    end
    
    if AttemptIsGood %change it and move on
        SimbioModel.Species(CurrentSpecies).rename(CurrentName);
        xuSyms(CurrentSpecies) = sym(CurrentName);
        xNicestrs{CurrentSpecies} = CurrentName;
        CurrentSpecies = CurrentSpecies + 1;
        CurrentAttempt = CurrentAttempt + 1;
    else %it failed; try again with different number
        CurrentAttempt = CurrentAttempt + 1;
        AttemptIsGood = true;
    end
end

%% Convert parameter names to things that won't confuse the symbolic toolbox
%The names will be "pa1x", "pa2x", etc...
kSyms = sym(zeros(nk,1));
kNicestrs = cell(nk,1);

CurrentAttempt = 1; %The numerical appendix that will be tried
CurrentParameter = 1; %Index through the name
AttemptIsGood = true; %becomes true when a name already exists

while(CurrentParameter <= nk) %Until every parameter
    CurrentName = sprintf('pa%dx', CurrentAttempt);
    %Check to make sure a compartment does not have our systematic name
    for check = 1:nv
        if strcmp(CurrentName, vNames{check})
            AttemptIsGood = false;
        end
    end
    %Check species names
    for check = 1:nx
        if strcmp(CurrentName, xuNames{check})
            AttemptIsGood = false;
        end
    end
    %Check parameter names
    for check = 1:nk
        if strcmp(CurrentName, kNames{check}) && CurrentParameter ~= check
            AttemptIsGood = false;
        end
    end
    
    if AttemptIsGood %change it and move on
        if CurrentParameter <= nkm %if it is a model paramter
            SimbioModel.Parameters(CurrentParameter).rename(CurrentName);
            
        else %otherwise cycle through the kinetic law parameters
            %This loop cycles through the reactions, finding how many
            %kinetic law parameters are in each. ind starts off as the
            %number of non-model parameters and is decremented by the
            %number of kinetic law parameters in each reaction until ind is
            %less than the kinetic law parameters defined in a specific
            %reaction. This means that ind now refers to the parameter that
            %needs to be changed.
            ind = CurrentParameter - nkm;
            ParameterNameNotChanged = true;
            CurrentReaction = 1;
            while(ParameterNameNotChanged)
                kineticlaw = reactions(CurrentReaction).KineticLaw; % Will be empty if this reaction does not use a kinetic law
                if ~isempty(kineticlaw)
                    nkkl = length(kineticlaw.Parameters);
                else
                    nkkl = 0;
                end
                if ind <= nkkl
                    kineticlaw.Parameters(ind).rename(CurrentName);
                    ParameterNameNotChanged = false;
                else
                    ind = ind - nkkl;
                    CurrentReaction = CurrentReaction + 1;
                end
            end
        end
        
        kSyms(CurrentParameter) = sym(CurrentName);
        kNicestrs{CurrentParameter} = CurrentName;
        CurrentParameter = CurrentParameter + 1;
        CurrentAttempt = CurrentAttempt + 1;
    else %it failed; try again with different number
        CurrentAttempt = CurrentAttempt + 1;
        AttemptIsGood = true;
    end
end

if opts.Verbose; fprintf('done.\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 3: Building the Diff Eqs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.Verbose; fprintf('Building the differential equations...'); end
%% Build up rate string and stoichiometry
reactions = SimbioModel.Reactions;
speciesDimension = SimbioModel.getconfigset.CompileOptions.DefaultSpeciesDimension;

nr = length(reactions);
S = sparse(nx, nr);
rNames = cell(nr,1);
rStrs = cell(nr,1);

%Get each reaction and build stochiometry matrix
for i = 1:nr %Foreach reaction
    % Get reaction name
    rNames{i} = reactions(i).Name;
    
    %Get the reaction rate
    rStrs{i,1} = reactions(i).Reactionrate;
    
    %Build the stochiometry matrix
    reactants  = reactions(i).Reactants;
    products   = reactions(i).Products;
    stoichio   = reactions(i).Stoichiometry; % = [Reactants, Products] in order

    nreact     = length(reactants);
    nProd      = length(products);

    for j = 1:nreact %Foreach reactant
        name   = reactants(j).Name;
        ind    = find(strcmp(xNicestrs, name));
        if strcmp(speciesDimension, 'substance')
            % both stoichiometry and species are in amount
            S(ind, i) = S(ind, i) + stoichio(j);
        else % speciesDimension == 'concentration'
            % stoichiometry is in concentration, reactions are in amount
            S(ind, i) = S(ind, i) + stoichio(j) / xCap(ind);
        end
    end
    
    for j = 1:nProd %Foreach product
        name   = products(j).Name;
        ind    = find(strcmp(xNicestrs, name));
        if strcmp(speciesDimension, 'substance')
            % both stoichiometry and species are in amount
            S(ind, i) = S(ind, i) + stoichio(j+nreact);
        else % speciesDimension == 'concentration'
            % stoichiometry is in concentration, reactions are in amount
            S(ind, i) = S(ind, i) + stoichio(j+nreact) / xCap(ind);
        end
    end
end

%% Assemble rules expressions
rules = SimbioModel.Rules;
nRules = numel(rules);

% Initialize empty symbolic vector
targetStrs = cell(nRules,1);
valueStrs  = cell(nRules,1);

for iRule = 1:nRules
    rule = rules(iRule);
    
    if strcmp(rule.RuleType, 'repeatedAssignment')
        % Split string rule into target and value
        splits = regexp(rule.Rule, '=', 'split');
        assert(numel(splits) == 2, 'KroneckerBio:simbio2Symbolic:InvalidRepeatedAssignment', 'Rule %i had an unparsible repeated assignment rule', iRule)
        targetStrs{iRule} = splits{1};
        valueStrs{iRule}  = splits{2};
    else
        error('KroneckerBio:simbio2Symbolic:UnsupportedRuleType', 'Kronecker only supports repeatedAssignment rules, rule %i has type %s', iRule, rule.RuleType)
    end
end

%% Convert reactions and rules to symbolics
%Delete compartment prefixes "co1x."
for i = 1:nv
    %'.' is a special character in regexprep, use '\.' to really mean '.'
    rStrs = regexprep(rStrs, [vNicestrs{i} '\.'], '');
end

% Replace string variables with symbolics
for iv = 1:nv
    rStrs = regexprep(rStrs, vNicestrs{iv}, ['sym(''' vNicestrs{iv} ''')']);
    valueStrs = regexprep(valueStrs, vNicestrs{iv}, ['sym(''' vNicestrs{iv} ''')']);
end
for ix = 1:nx
    rStrs = regexprep(rStrs, xNicestrs{ix}, ['sym(''' xNicestrs{ix} ''')']);
    valueStrs = regexprep(valueStrs, xNicestrs{ix}, ['sym(''' xNicestrs{ix} ''')']);
end
for ik = 1:nk
    rStrs = regexprep(rStrs, kNicestrs{ik}, ['sym(''' kNicestrs{ik} ''')']);
    valueStrs = regexprep(valueStrs, kNicestrs{ik}, ['sym(''' kNicestrs{ik} ''')']);
end
for iRule = 1:nRules
    rStrs = regexprep(rStrs, targetStrs{iRule}, ['sym(''' targetStrs{iRule} ''')']);
    valueStrs = regexprep(valueStrs, targetStrs{iRule}, ['sym(''' targetStrs{iRule} ''')']);
end

%Create symbolic reaction rates
rSyms = repmat(sym('empty'), nr,1);
for ir = 1:nr
    rSyms(ir) = eval(rStrs{ir});
end

%Create symbolic rules
targetSyms = repmat(sym('empty'), nRules,1);
valueSyms  = repmat(sym('empty'), nRules,1);
for iRule = 1:nRules
    targetSyms(iRule)  = sym(targetStrs{iRule});
    valueSyms(iRule) = eval(valueStrs{iRule});
end

%% Replace reaction rates with symbolics
% This may require up to nRules iterations of substitution
for iRule = 1:nRules
    rSyms = subs(rSyms, targetSyms, valueSyms, 0);
end

% Delete rule parameters
found = lookup(targetSyms, kSyms);
kSyms(found(found ~= 0)) = [];
kNames(found(found ~= 0)) = [];
k(found(found ~= 0)) = [];
nk = numel(kSyms);

% Convert all rule species to inputs
found = lookup(targetSyms, xuSyms);
isu(found(found ~= 0)) = true;

%% Generate the ode of the model
%Odes of the model for all variables
f = S*rSyms;

if opts.Verbose; fprintf('done.\n'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 4: Build Symbolic Model %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SymModel.Type       = 'Model.Symbolic';
SymModel.Name       = Name;

SymModel.nv         = nv;
SymModel.nxu        = nx;
SymModel.nk         = nk;
SymModel.nr         = nr;

SymModel.vSyms      = vSyms;
SymModel.vNames     = vNames;
SymModel.dv         = zeros(nv,1) + 3; % TODO: from units
SymModel.vValues    = vValues;

SymModel.xuSyms     = xuSyms;    %sym nxx1
SymModel.xuNames    = xuNames;   %cell nxx1
SymModel.vxuNames   = vxuNames;
SymModel.isu        = isu;
SymModel.xu0        = xu0;       %num nxx1

SymModel.kSyms      = kSyms;    %sym nxx1
SymModel.kNames     = kNames;   %cell nkx1
SymModel.k          = k;        %num nxx1

SymModel.rNames     = rNames;
SymModel.r          = rSyms;    %sym nrx1
SymModel.S          = S;        %num nxxnr

SymModel.f          = f;        %sym nxx1
