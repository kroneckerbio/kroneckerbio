function symModel = sbml2Symbolic(sbmlModel, opts)
%SBML2SYMBOLIC Inport SBML model and covert to kroneckerbio symbolic model.
%   Detailed explanation goes here
% Currently based heavily on simbio2Symbolic (which MathWorks clearly adpated the libSBML API)
% Note: Doesn't convert variable names in Part 2, assuming that SBML IDs
%   are all allowed symbolic variable names.
%   Uses IDs as the internal identifier in this function.

% TEST: convert SBML -> SimBio -> symbolic for comparison; remove when done
% When done, deprecate the sbmlimport function (it depends on the SimBiology toolbox)
% symbolicTest = simbio2Symbolic(sbmlimport(sbmlModel));

%% Options
% Resolve missing inputs
if nargin < 2
    opts = [];
end

% Options for displaying progress
defaultOpts.Verbose = 0;
defaultOpts.Validate = false;
defaultOpts.UseNames = false;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Call libSBML to import SBML model
if verbose; fprintf('Convert SBML model using libSBML...\n'); end

sbml = TranslateSBML(sbmlModel, double(opts.Validate), opts.Verbose);

if verbose; fprintf('done.\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 1: Extracting the Model Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: In SBML, `id` is a required globally scoped reference while `name`
% is an optional human readable identifier. This function uses ids throughout for
% parsing/converting and outputs ids by default. You can specify that ids
% be converted to names with the option ConvertIDs. However, this may
% introduce ambiguity and/or errors if names aren't valid Matlab strings or
% breaks scoping rules that I don't know about.
if verbose; fprintf('Extracting model components...'); end

%% Model name
if ~isempty(sbml.name) % In this case, name is more useful if present
    name = sbml.name;
elseif ~isempty(sbml.id)
    name = sbml.id;
else
    name = 'Inported SBML model';
end

%% Compartments
nv = length(sbml.compartment);
vIDs   = cell(nv,1);
vNames = cell(nv, 1);
v      = zeros(nv,1);

for iv = 1:nv
    vIDs{iv}   = sbml.compartment(iv).id;
    vNames{iv} = sbml.compartment(iv).name;
    
    if sbml.compartment(iv).isSetSize
        v(iv) = sbml.compartment(iv).size;
    else
        % TODO: Consider variable compartment size by assignment or fitting
        warning('Warning:sbml2Symbolic:CompartmetSizeNotSet: Compartment size not set, setting default size = 1.')
        v(iv) = 1;
    end
    
end

%% Species
nxu = length(sbml.species);
xuIDs   = cell(nxu,1);
xuNames = cell(nxu,1);
xu0     = zeros(nxu,1);
vxuInd  = zeros(nxu,1);
isu     = false(nxu,1);
xuSubstanceUnits = false(nxu,1);

for ixu = 1:nxu
    species = sbml.species(ixu);
    
    xuIDs{ixu}   = species.id;
    xuNames{ixu} = species.name;
    
    if species.isSetInitialAmount
        xu0(ixu) = species.initialAmount;
    elseif species.isSetInitialConcentration
        xu0(ixu) = species.initialConcentration;
    else
        warning('sbml2Symbolic:InitialConcentrationNotSet: Initial species conc. not set for %s, setting default conc. = 0.', xuIDs{ixu})
        xu0(ixu) = 0;
    end
    
    % Get species compartment by id
    vxuInd(ixu) = find(strcmp(species.compartment, vIDs));
    
    % Species is input/conc. doesn't change due to reactions, etc.
    isu(ixu) = species.boundaryCondition || species.constant;
    
    % Species substance units in amount/true or conc./false
    xuSubstanceUnits(ixu) = logical(species.hasOnlySubstanceUnits);
end

%% Parameters
nk  = length(sbml.parameter); % Total number of parameters (may change)
nkm = nk; % Number of model parameters
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
for ir = 1:nr
    kineticLaw = sbml.reaction(ir).kineticLaw; % Will be empty if no kinetic law parameters exist
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 2: Convert to symbolics %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose; fprintf('Converting to symbolics...'); end

%% Compartments
vSyms = sym(zeros(nv,1));
for i = 1:nv
    vSyms(i) = sym(vIDs{i});
end

%% Species
xuSyms = sym(zeros(nxu,1));
for i = 1:nxu
    xuSyms(i) = sym(xuIDs{i});
end

%% Parameters
kSyms = sym(zeros(nk,1));
for i = 1:nk
    kSyms(i) = sym(kIDs{i});
end
%% Done converting to symbolics
if verbose; fprintf('done.\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 3: Building the Diff Eqs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose; fprintf('Building diff eqs...'); end

%% Rules
% Assume s.rule field contains repeated assignments - why?
% s.rule field always exists, but can be empty
nRules = length(sbml.rule);

if nRules > 0
    warning('sbml2Symbolic: Only repeated assignment and initial assignment rules are supported at this time.')
end

% Initialize empty symbolic vector
targetStrs = cell(nRules,1);
valueStrs  = cell(nRules,1);

% Assign repeated ruels
for iRule = 1:nRules
    targetStrs{iRule} = sbml.rule(iRule).variable;
    valueStrs{iRule}  = sbml.rule(iRule).formula;
end

% s.initialAssignment field is separate if it exists - why?
if isfield(sbml, 'initialAssignment')
    warning('sbml2Symbolic:TODO: Initial assignment rules enforced only at model initialization, not after modifying model.')
    
    nInitialAssignments = length(sbml.initialAssignment);
    
    targetStrsInit = cell(nRules,1);
    valueStrsInit  = cell(nRules,1);
    
    for i = 1:nInitialAssignments
        targetStrsInit{i} = sbml.initialAssignment(i).symbol;
        valueStrsInit{i}  = sbml.initialAssignment(i).math;
    end
    
    targetStrs = [targetStrs; targetStrsInit];
    valueStrs  = [valueStrs;  valueStrsInit ];
    
    nRules = nRules + nInitialAssignments;
end

%% States and Outputs
nx = nnz(~isu);
xIDs   = xuIDs(~isu);
xNames = xuNames(~isu);
xSyms  = xuSyms(~isu);
x0     = sym(xu0(~isu)); % Default is no seed
vxInd  = vxuInd(~isu);

ns = 0;
sSyms = sym(zeros(0,1));
sIDs   = cell(0,1);
sNames = cell(0,1);
s = zeros(0,1);

nu = nnz(isu);
uIDs   = xuIDs(isu);
uNames = xuNames(isu);
uSyms  = xuSyms(isu);
u = sym(xu0(isu)); % Default is no time varying inputs
vuInd   = vxuInd(isu);

% Input parameters don't have an analog in SBML
nq = 0;
qSyms = sym(zeros(0,1));
qIDs   = cell(0,1);
qNames = cell(0,1);
q = zeros(0,1);

%% Reactions
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
    rStrs{i,1} = reaction.kineticLaw.math; % check this or formula
    
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

%% Convert reactions and rules to symbolics
% TODO: Currently doesn't work in simbio2Symbolic
%   Need to put contents of rules in scope
% Create symbolic rules
targetSyms = sym(zeros(nRules,1));
valueSyms  = sym(zeros(nRules,1));
for iRule = 1:nRules
    targetSyms(iRule) = sym(targetStrs{iRule});
    valueSyms(iRule)  = eval(valueStrs{iRule});
end

%% Replace reaction rates with symbolics
% This may require up to nRules iterations of substitution
for iRule = 1:nRules
    r = subs(r, targetSyms, valueSyms, 0);
end

% Delete rule parameters
found = lookup(targetSyms, kSyms);
kSyms(found(found ~= 0)) = [];
kIDs(found(found ~= 0)) = [];
k(found(found ~= 0)) = [];
nk = numel(kSyms);

% Convert all rule species to inputs
found = lookup(targetSyms, xuSyms);

%% Done building diff eqs
if verbose; fprintf('done.\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 5: Use human readable names instead of IDs if desired %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.UseNames
    if verbose; fprintf('Converting IDs to names...'); end
    vIDs = vNames;
    kIDs = kNames;
    sIDs = sNames;
    qIDs = qNames;
    uIDs = uNames;
    xIDs = xNames;
    rIDs = rNames;
    if verbose; fprintf('done.\n'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
symModel.vNames     = vIDs;
symModel.dv         = zeros(nv,1) + 3; % TODO: from units
symModel.v          = v;

symModel.kSyms      = kSyms;
symModel.kNames     = kIDs;
symModel.k          = k;

symModel.sSyms      = sSyms;
symModel.sNames     = sIDs;
symModel.s          = s;

symModel.qSyms      = qSyms;
symModel.qNames     = qIDs;
symModel.q          = q;

symModel.uSyms      = uSyms;
symModel.uNames     = uIDs;
symModel.vuInd      = vuInd;
symModel.u          = u;

symModel.xSyms      = xSyms;
symModel.xNames     = xIDs;
symModel.vxInd      = vxInd;
symModel.x0         = x0;

symModel.rNames     = rIDs;
symModel.r          = r;
symModel.S          = S;
symModel.Su         = Su;
