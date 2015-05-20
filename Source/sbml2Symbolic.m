function symModel = sbml2Symbolic(sbmlModel, opts)
%SBML2SYMBOLIC Inport SBML model and covert to kroneckerbio symbolic model.

%   Inputs
%   SimbioModel: [ Simbiology model scalar ]
%       A Simbiology model object
%   opts: [ options struct scalar {} ]
%       .Verbose [ logical scalar {false} ]
%       	Print progress to command window
%       .UseNames [logical scalar {false} ]
%           Whether to convert SBML IDs to Names and autogenerate new IDs
%           Use this when the supplied SBML model uses "nice" names as IDs
%       .EvaluateExternalFunctions [ logical scalar {false} ]
%           Determines whether to evaluate calls to external functions in
%           the reaction rate expressions. The external functions are
%           evaluated with symbolic input arguments to obtain a symbolic
%           version of the output that can be differentiated. If set to
%           false, derivatives of external function calls cannot be
%           computed.
%
%   Outputs
%   SymModel
%       .Type: [ 'Model.SymbolicReactions' ]
%       .Name: [ string ]
%           Copied from SimbioModel.Name
%       .nv: [ nonegative integer scalar ]
%           Number of compartments
%       .nk: [ nonegative integer scalar ]
%           Number of kinetic parameters
%       .ns: [ nonegative integer scalar ]
%           Number of seed parameters
%       .nq: [ nonegative integer scalar ]
%           Number of input control parameters
%       .nu: [ nonegative integer scalar ]
%           Number of inputs
%       .nx: [ nonegative integer scalar ]
%           Number of states
%       .nr: [ nonegative integer scalar ]
%           Number of reactions
%       .vSyms: [ symbolic vector nv ]
%           Symbolic name of each compartment
%       .vNames: [ cell vector of strings nv ]
%           Natural names of the compartments
%       .v: [ positive vector nv ]
%           Sizes of the compartments
%       .kSyms: [ symbolic vector nk ]
%           Symbolic name of each kinetic parameter
%       .kNames: [ cell vector nk of strings ]
%           Natural names of the kinetic parameters
%       .k: [ nonegative vector nk ]
%           Kinetic parameter values
%       .sNames: [ cell vector of strings ns ]
%           Natural names of the seed parameters
%       .s: [ nonegative vector ns ]
%           Seed parameter values
%       .q: [ nonegative vector nq ]
%           Input control parameter values
%       .uSyms: [ symbolic vector nu ]
%           Symbolic name of each input species
%       .uNames: [ cell vector of strings nu ]
%           Natural names of the input species
%       .uInd [ positive integer vector nu ]
%           Index of the compartment to which the inputs belong
%       .u: [ symbolic vector nu ]
%           Symbolic representation of each input species
%       .xSyms: [ symbolic vector nx ]
%           Symbolic name of each state species
%       .xInd [ positive integer vector nx ]
%           Index of the compartment to which the states belong
%       .xNames: [ cell vector of strings nx ]
%           Natural names of the state species
%       .dx0ds: [ nonegative matrix nx by ns ]
%           The influence of each seed parameter on each state species's
%           initial amount
%       .x0 [ symbolic vector nx ]
%           Initial conditions of the state species
%       .r [ symbolic vector nr ]
%           Symbolic representation of each reaction rate using the
%           symbolic species and parameters
%       .S [ matrix nx by nr ]
%           The stoichiometry matrix of the reactions
%       .Su [ matrix nx by nr ]
%           A stoichiometry matrix of how the reactions are trying to alter
%           the inputs

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Options
% Resolve missing inputs
if nargin < 2
    opts = [];
end

% Options for displaying progress
defaultOpts.Verbose = 0;
defaultOpts.Validate = false;
defaultOpts.EvaluateExternalFunctions = false;
defaultOpts.UseNames = false;
defaultOpts.AddRulesAsOutputs = false;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Call libSBML to import SBML model
if verbose; fprintf('Convert SBML model using libSBML...\n'); end

sbml = TranslateSBML(sbmlModel, double(opts.Validate), opts.Verbose);

if verbose; fprintf('done.\n'); end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 1: Extracting the Model Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose; fprintf('Extracting model components...'); end

%% Model name
if opts.UseNames
    name = sbml.id;
else
    name = sbml.name;
end

%% Compartments
nv = length(sbml.compartment);
vIDs   = cell(nv,1);
vNames = cell(nv,1);
v      = zeros(nv,1);
dv     = zeros(nv,1);

for i = 1:nv
    vIDs{i}   = sbml.compartment(i).id;
    vNames{i} = sbml.compartment(i).name;
    
    if sbml.compartment(i).isSetSize
        v(i) = sbml.compartment(i).size;
    else
        % TODO: Consider variable compartment size by assignment or fitting
        warning('Warning:sbml2Symbolic:CompartmentSizeNotSet: Compartment size not set, setting default size = 1.')
        v(i) = 1;
    end
    
    dv(i) = sbml.compartment(i).spatialDimensions;
    
end

%% Species
nxu = length(sbml.species);
xuIDs   = cell(nxu,1);
xuNames = cell(nxu,1);
xu0     = zeros(nxu,1);
vxuInd  = zeros(nxu,1);
isu     = false(nxu,1);
xuSubstanceUnits = false(nxu,1);

for i = 1:nxu
    species = sbml.species(i);
    
    xuIDs{i}   = species.id;
    xuNames{i} = species.name;
    
    if species.isSetInitialAmount
        xu0(i) = species.initialAmount;
    elseif species.isSetInitialConcentration
        xu0(i) = species.initialConcentration;
    else
        warning('sbml2Symbolic:InitialConcentrationNotSet: Initial species conc. not set for %s, setting default conc. = 0.', xuIDs{i})
        xu0(i) = 0;
    end
    
    % Get species compartment by id
    vxuInd(i) = find(strcmp(species.compartment, vIDs));
    
    % Species is input/conc. doesn't change due to reactions, etc.
    isu(i) = species.boundaryCondition || species.constant;
    
    % Species substance units in amount/true or conc./false
    xuSubstanceUnits(i) = logical(species.hasOnlySubstanceUnits);
end

%% Parameters
nk  = length(sbml.parameter); % Total number of parameters (may change)
kIDs   = cell(nk,1);
kNames = cell(nk,1);
k      = zeros(nk,1);

% Get model parameters
for ik = 1:nk
    kIDs{ik}   = sbml.parameter(ik).id;
    kNames{ik} = sbml.parameter(ik).name;
    k(ik)      = sbml.parameter(ik).value;
end

%% Reactions
nr = length(sbml.reaction);

% Local parameters
for i = 1:nr
    kineticLaw = sbml.reaction(i).kineticLaw; % Will be empty if no kinetic law parameters exist
    if ~isempty(kineticLaw)
        constantskl = kineticLaw.parameter; % Will only fetch parameters unique to this kinetic law
        nkkl = length(constantskl);
        for j = 1:nkkl
            nk = nk + 1; % Add one more parameter
            kIDs{nk,1}   = constantskl(j).id;
            kNames{nk,1} = constantskl(j).name;
            k(nk,1)      = constantskl(j).value;
        end
    end
end

%% Done extracting model variables
if verbose; fprintf('done.\n'); end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 2: Convert to symbolics %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose; fprintf('Converting to symbolics...'); end

%% Compartments
vSyms = sym(vIDs);

%% Species
xuSyms = sym(xuIDs);

%% Parameters
kSyms = sym(kIDs);

%% Done converting to symbolics
if verbose; fprintf('done.\n'); end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 3: Building the Diff Eqs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose; fprintf('Building diff eqs...'); end

%% Assemble rules expressions
% Orders ruls by repeated, then initial assignments
% Note: libSBML breaks original ordering?
nRules = length(sbml.rule);

% Rule types:
%   0 = repeated assignment
%   1 = initial assignment
assignmentTypes = zeros(nRules,1);
if isfield(sbml, 'initialAssignment')
    nInitialAssignments = length(sbml.initialAssignment);
    nRules = nRules + nInitialAssignments;
    assignmentTypes = [assignmentTypes; ones(nInitialAssignments,1)];
end

% Store targets and values
targetStrs = {sbml.rule.variable}';
valueStrs  = {sbml.rule.formula}' ;

if isfield(sbml, 'initialAssignment')
    targetStrs = [targetStrs; {sbml.initialAssignment.symbol}'];
    valueStrs  = [valueStrs;  {sbml.initialAssignment.math}'];
end

% Turn into symbolics
targetSyms = sym(targetStrs);
valueSyms  = sym(valueStrs);

% Make lookup function indicating which values are constant
allIDs = [xuIDs; kIDs; vIDs];
allConstants = [isu; logical([sbml.parameter.constant]'); logical([sbml.compartment.constant]')];
isConstant = @(var) allConstants(strcmp(var, allIDs));

substitute = false(nRules,1);
makeoutput = false(nRules,1);
setseedvalue = false(nRules,1);
setparametervalue = false(nRules,1);
setcompartmentvalue = false(nRules,1);

for i = 1:nRules
    
    if assignmentTypes(i) == 0 % Repeated assignments
        
        substitute(i) = true;
        makeoutput(i) = true;
        
    elseif assignmentTypes(i) == 1 % Initial assignments
        
        % Get strings of all variables in valueSyms
        thisvalueStrs = arrayfun(@char, symvar(valueSyms(i)), 'UniformOutput', false);
        
        % Check whether initialAssignment is from constants to constants.
        % Initial assignment to constants from constants can be treated the
        % same as repeatedAssignment, since neither the assignees or
        % assigners change with time. Initial assignments to species from
        % constants can be treated as seed parameters if the assigner (1)
        % only assigns to species and (2) only appears with other seed
        % parameters. Any other case will not work quite the same as it
        % does in SimBiology, since there is currently no way in
        % KroneckerBio to enforce the rule's constraint after the model is
        % built.
        
        % Determine whether target and values are constants
        targetIsConstant = isConstant(targetStrs{i});
        valueIsConstant = arrayfun(isConstant, thisvalueStrs);
        
        % Determine variable types of target and values
        targetIsSpecies       = any(strcmp(targetStrs{i}, xuIDs));
        targetIsParameter     = any(strcmp(targetStrs{i}, kIDs));
        targetIsCompartment   = any(strcmp(targetStrs{i}, vIDs));
%         valuesAreParameters   = any(cellfun(@strcmp, repmat(thisvalueStrs(:), 1, nk), repmat(kIDs(:)', length(thisvalueStrs), 1)), 2);
%         valuesAreSpecies      = any(cellfun(@strcmp, repmat(thisvalueStrs(:), 1, nxu), repmat(xuIDs(:)', length(thisvalueStrs), 1)), 2);
%         valuesAreCompartments = any(cellfun(@strcmp, repmat(thisvalueStrs(:), 1, nv), repmat(vIDs(:)', length(thisvalueStrs), 1)), 2);
        
        % If all the associated values are constants...
        if targetIsConstant && all(valueIsConstant)
            
            % Perform substitution to enforce the rule
            substitute(i) = true;
            
            % If the target is a species, set up an output. Otherwise
            % don't.
            if targetIsSpecies
                makeoutput(i) = true;
            end
            
        else % If some values are not constant...
            
            warning([targetStrs{i} ' will be set to ' valueStrs{i} ' initially, but if ' strjoin(thisvalueStrs, ',') ' is/are changed following model initialization, ' targetStrs{i} ' must be updated manually to comply with the rule.'])
            
            if targetIsSpecies
                setseedvalue(i) = true;
            elseif targetIsParameter
                setparametervalue(i) = true;
            elseif targetIsCompartment
                setcompartmentvalue(i) = true;
            end
            
        end
    end
    
end

%% States and Outputs
nx = nnz(~isu);
xIDs   = xuIDs(~isu);
xNames = xuNames(~isu);
xSyms  = xuSyms(~isu);
vxInd  = vxuInd(~isu);

% Represent every state's initial condition with a seed
ns = nx;
if opts.UseNames
    sNames = xIDs;
else
    sNames = xNames;
end
sIDs = cell(ns,1);
for i = 1:ns
    sIDs{i} = genUID;
end
sSyms  = sym(sIDs);
s      = xu0(~isu);
x0     = sSyms; % initial conditions syms which can be updated in experiments

nu = nnz(isu);
uIDs   = xuIDs(isu);
uNames = xuNames(isu);
uSyms  = xuSyms(isu);
u      = xu0(isu);
vuInd  = vxuInd(isu);

% Input parameters don't have an analog in SBML
nq = 0;
qSyms = sym(zeros(0,1));
qIDs   = cell(0,1);
qNames = cell(0,1);
q = zeros(0,1);

%% Reactions
% Need to apply assignment rules to rate forms
nSEntries = 0;
SEntries  = zeros(0,3);
rIDs   = cell(nr,1);
rNames = cell(nr,1);
rStrs  = cell(nr,1);

% Get each reaction and build stochiometry matrix
for i = 1:nr
    reaction = sbml.reaction(i);
    
    % Get reaction name
    rIDs{i}   = reaction.id;
    rNames{i} = reaction.name;
    
    % Get reaction rate
    %   Formula and math are similar, but in 1 case, formula uses "power" while math uses "pow"
    rStrs{i,1} = reaction.kineticLaw.formula; % check this or formula
    
    % Tally new entries
    nReactants = length(reaction.reactant);
    nProducts = length(reaction.product);
    nAdd = nReactants + nProducts;
    
    % Add more room in vector if necessary
    currentLength = size(SEntries,1);
    if nSEntries + nAdd > currentLength
        addLength = max(currentLength, 1);
        SEntries = [SEntries; zeros(addLength,3)];
    end
    
    % Build stoichiometry matrix
    for j = 1:nReactants
        reactant = reaction.reactant(j).species;
        stoich = -reaction.reactant(j).stoichiometry;
        ind = find(strcmp(xuIDs, reactant));
        
        nSEntries = nSEntries + 1;
        SEntries(nSEntries,1) = ind;
        SEntries(nSEntries,2) = i;
        
        if xuSubstanceUnits(ind) % Both stoichiometry and species are in amount
            SEntries(nSEntries,3) = stoich;
        else % Stoichiometry is in concentration, reactions are in amount
            SEntries(nSEntries,3) = stoich / v(vxuInd(ind));
        end
        
    end
    
    for j = 1:nProducts
        product = reaction.product(j).species;
        stoich = reaction.product(j).stoichiometry;
        ind = find(strcmp(xuIDs, product));
        
        nSEntries = nSEntries + 1;
        SEntries(nSEntries,1) = ind;
        SEntries(nSEntries,2) = i;
        
        if xuSubstanceUnits(ind) % Both stoichiometry and species are in amount
            SEntries(nSEntries,3) = stoich;
        else % Stoichiometry is in concentration, reactions are in amount
            SEntries(nSEntries,3) = stoich / v(vxuInd(ind));
        end
    end
end

% Symbolically evaluate r
r = sym(rStrs);

% Assemble stoichiometry matrix
S = sparse(SEntries(1:nSEntries,1), SEntries(1:nSEntries,2), SEntries(1:nSEntries,3), nxu, nr);
Su = S(isu,:);
S = S(~isu,:);

%% Substitute assignment rules into reaction rates
% This may require up to nRules iterations of substitution
% nRules^2 time complexity overall
subsRules = find(substitute);
for i = subsRules(:)'
    r = subs(r, targetSyms(subsRules), valueSyms(subsRules));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 5: Cleanup symbolic names (if specified) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If UseNames = true, write IDs into Names and generate UUIDS for IDs
%   Useful if SBML file's IDs potentially lead to name collisions
if opts.UseNames
    if verbose; fprintf(' generating unique IDs...'); end
    
    %% Generate unique IDs
    % Compartments
    vNames = vIDs;
    for i = 1:nv
        vIDs{i} = genUID;
    end
    
    % States
    xNames = xIDs;
    for i = 1:nx
        xIDs{i} = genUID;
    end
    
    % Inputs
    uNames = uIDs;
    for i = 1:nu
        uIDs{i} = genUID;
    end
    
    % Parameters
    kNames = kIDs;
    for i = 1:nk
        kIDs{i} = genUID;
    end
    
    % Seed parameters are already initialized with unique names
    
    % Input parameters aren't implemented
    
    %% Replace IDs in sym vars everywhere with new unique IDs
    oldSyms = [vSyms; xSyms; uSyms; kSyms];
    newIDs = [vIDs; xIDs; uIDs; kIDs];
    newSyms = sym(newIDs);
    
    vSyms = subs(vSyms, oldSyms, newSyms);
    xSyms = subs(xSyms, oldSyms, newSyms);
    uSyms = subs(uSyms, oldSyms, newSyms);
    kSyms = subs(kSyms, oldSyms, newSyms);
    r = subs(r, oldSyms, newSyms);
    
end

% Reevaluate terms to change external functions to symbolic expressions
if opts.EvaluateExternalFunctions
    r = evaluateExternalFunctions(r, [vIDs; xIDs; uIDs; kIDs]);
end

% Delete rule parameters
[kSyms, kNames, k, nk] = deleteRuleParameters(kSyms, kNames, k, targetSyms, substitute);
[xSyms, xNames, s, nx, found] = deleteRuleParameters(xSyms, xNames, s, targetSyms, substitute);
sSyms(found(found ~= 0)) = [];
sNames(found(found ~= 0)) = [];
[uSyms, uNames, u, nu] = deleteRuleParameters(uSyms, uNames, u, targetSyms, substitute);

% Convert rule terms to outputs
%   This is optional but convenient
%   Make sure to add additional outputs as desired when building model
if opts.AddRulesAsOutputs
    y = valueSyms(makeoutput);
    yNames = arrayfun(@char, targetSyms, 'UniformOutput', false);
else
    y = sym([]);
    yNames = {};
end

% Substitute values for initial assignments
for i = 1:nRules
    % If this is an initial assignment rule...
    if setseedvalue(i) || setparametervalue(i) || setcompartmentvalue(i)
        target = targetSyms(i);
        valuestoassign = valueSyms(i);
        valuestoassign = subs(valuestoassign, [xSyms; uSyms; kSyms; vSyms], [s; u; k; v]);
        if setseedvalue(i)
            seedtargets_i = find(logical(target == xSyms));
            if isempty(seedtargets_i)
                inputtargets_i = find(logical(target == uSyms));
                u(inputtargets_i) = double(valuestoassign);
            else
                s(seedtargets_i) = double(valuestoassign);
            end
        elseif setparametervalue(i)
            paramtargets_i = find(logical(target == kSyms));
            k(paramtargets_i) = double(valuestoassign);
        elseif setcompartmentvalue(i)
            compartmenttargets_i = find(logical(target == vSyms));
            v(compartmenttargets_i) = double(valuestoassign);
        end
    end
end

%% Done building diff eqs
if verbose; fprintf('done.\n'); end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 6: Build Symbolic Model %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symModel.Type       = 'Model.SymbolicReactions';
symModel.Name       = name;

symModel.nv         = nv;
symModel.nk         = nk;
symModel.ns         = ns;
symModel.nq         = nq;
symModel.nu         = nu;
symModel.nx         = nx;
symModel.nr         = nr;

symModel.vSyms      = vSyms;
symModel.vNames     = vNames;
symModel.dv         = dv;
symModel.v          = v;

symModel.kSyms      = kSyms;
symModel.kNames     = kNames;
symModel.k          = k;

symModel.sSyms      = sSyms;
symModel.sNames     = sNames;
symModel.s          = s;

symModel.qSyms      = qSyms;
symModel.qNames     = qNames;
symModel.q          = q;

symModel.uSyms      = uSyms;
symModel.uNames     = uNames;
symModel.vuInd      = vuInd;
symModel.u          = u;

symModel.xSyms      = xSyms;
symModel.xNames     = xNames;
symModel.vxInd      = vxInd;
symModel.x0         = x0;

symModel.rNames     = rNames;
symModel.r          = r;
symModel.S          = S;
symModel.Su         = Su;

symModel.y          = y;
symModel.yNames     = yNames;

end

%% %%%%%%%%%%%%%%%%%
% Helper Functions %
%%%%%%%%%%%%%%%%%%%%
function [kSyms, kNames, k, nk, found] = deleteRuleParameters(kSyms, kNames, k, targetSyms, substitute)

if isempty(kSyms)
    nk = numel(kSyms);
    found = [];
    return
end

found = lookup(targetSyms(substitute), kSyms);
kSyms(found(found ~= 0)) = [];
kNames(found(found ~= 0)) = [];
k(found(found ~= 0)) = [];
nk = numel(kSyms);

end

function rOut = evaluateExternalFunctions(rIn, ids)
% Evaluate symbolic functions/pull in functions defined in path
%   Necessary for "power" and other MathML function translation
% Initialize symbolic variables
syms(ids{:});

% Evaluate the expressions to remove function calls
rOut = eval(rIn);

% Clear the symbolic variables
% clear(ids{:})
end