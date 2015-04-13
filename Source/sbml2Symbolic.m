function SymModel = sbml2Symbolic(sbmlModel, opts)
%SBML2SYMBOLIC Inport SBML model and covert to kroneckerbio symbolic model.
%   Detailed explanation goes here
% Currently based heavily on simbio2Symbolic (which MathWorks clearly adpated the libSBML API)
% Note: Doesn't convert variable names in Part 2, assuming that SBML IDs
% are all allowed symbolic variable names (check this).
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
defaultOpts.ConvertIDs = false;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

%% Call libSBML to import SBML model
if verbose; fprintf('Convert SBML model using libSBML...\n'); end
sbml = TranslateSBML(sbmlModel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 1: Extracting the Model Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: In SBML, `id` is a required globally scoped reference while `name`
% is an optional human readable identifier. This function uses ids throughout for
% parsing/converting and outputs ids by default. You can specify that ids
% be converted to names with the option ConvertIDs. However, this may
% introduce ambiguity and/or errors if names aren't valid Matlab strings or
% breaks scoping rules that I don't know about.
if verbose; fprintf('Extracting model components...\n'); end

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
vNames  = cell(nv,1);
v = zeros(nv,1);

for iv = 1:nv
    vNames{iv} = sbml.compartment(iv).id;
    
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
xuNames  = cell(nxu,1);
xu0      = zeros(nxu,1);
vxuInd   = zeros(nxu,1);
isu      = false(nxu,1);
xuSubstanceUnits = false(nxu,1);

for ixu = 1:nxu
    species = sbml.species(ixu);
    
    xuNames{ixu} = species.id;
    
    if species.isSetInitialAmount
        xu0(ixu) = species.initialAmount;
    elseif species.isSetInitialConcentration
        xu0(ixu) = species.initialConcentration;
    else
        warning('sbml2Symbolic:InitialConcentrationNotSet: Initial species conc. not set for %s, setting default conc. = 0.', xuNames{ixu})
        xu0(ixu) = 0;
    end
    
    % Get species compartment by id
    vxuInd(ixu) = find(strcmp(species.compartment, vNames));
    
    % Species is input/conc. doesn't change due to reactions, etc.
    isu(ixu) = species.boundaryCondition || species.constant;
    
    % Species substance units in amount/true or conc./false
    xuSubstanceUnits(ixu) = logical(species.hasOnlySubstanceUnits);
end

%% Parameters
nk     = length(sbml.parameter); % Total number of parameters (may change)
nkm    = nk; % Number of model parameters
kNames = cell(nk,1);
k      = zeros(nk,1);

% Get model parameters
for ik = 1:nk
    kNames{ik} = sbml.parameter(ik).id;
    k(ik) = sbml.parameter(ik).value;
end

%% Reactions
nr = length(sbml.reaction);

for ir = 1:nr
    kineticLaw = sbml.reaction(ir).kineticLaw; % Will be empty if no kinetic law parameters exist
    if ~isempty(kineticLaw)
        constantskl = kineticLaw.parameter; % Will only fetch parameters unique to this kinetic law
        nkkl = length(constantskl);
        for j = 1:nkkl
            nk = nk + 1; % Add one more parameter
            kNames{nk,1} = constantskl(j).id;
            k(nk,1) = constantskl(j).value;
        end
    end
end

%% Done extracting model variables
if verbose; fprintf('done.\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 2: Convert to symbolics %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shouldn't need to convert names - all ids should be valid symbolic variable names

%% Compartments
vSyms = sym(zeros(nv,1));
for i = 1:nv
    vSyms(i) = sym(vNames{i});
end

%% Species
xuSyms = sym(zeros(nxu,1));
for i = 1:nxu
    xuSyms(i) = sym(xuNames{i});
end

%% Parameters
kSyms = sym(zeros(nk,1));
for i = 1:nk
    kSyms(i) = sym(kNames{i});
end

%% Done converting to symbolics
if verbose; fprintf('done.\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 3: Building the Diff Eqs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
xNames = xuNames(~isu);
xSyms  = xuSyms(~isu);
x0      = sym(xu0(~isu)); % Default is no seed
vxInd   = vxuInd(~isu);

ns = 0;
sSyms = sym(zeros(0,1));
sNames = cell(0,1);
s = zeros(0,1);

nu = nnz(isu);
uNames = xuNames(isu);
uSyms  = xuSyms(isu);
u = sym(xu0(isu)); % Default is no time varying inputs
vuInd   = vxuInd(isu);

nq = 0;
qSyms = sym(zeros(0,1));
qNames = cell(0,1);
q = zeros(0,1);

%% Reactions
nSEntries = 0;
SEntries  = zeros(0,3);
rNames    = cell(nr,1);
rStrs     = cell(nr,1);

% Get each reaction and build stochiometry matrix
for i = 1:nr
    reaction = sbml.reaction(i);
    
    % Get reaction name
    rNames{i} = reaction.id;
    
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
        ind = find(strcmp(xuNames, reactant));
        
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
        ind = find(strcmp(xuNames, product));
        
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
kNames(found(found ~= 0)) = [];
k(found(found ~= 0)) = [];
nk = numel(kSyms);

% Convert all rule species to inputs
found = lookup(targetSyms, xuSyms);

%% Done building diff eqs
if verbose; fprintf('done.\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 5: Convert IDs to human readable names if desired %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert all IDs to names
if opts.ConvertIDs
    error('Error:NotImplemented: Converting IDs to human readable names not implemented yet.')
    fprintf('Converting IDs to names...\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Part 5: Build Symbolic Model %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SymModel.Type       = 'Model.SymbolicReactions';
SymModel.Name       = name;

SymModel.nv         = nv;
SymModel.nk         = nk;
SymModel.ns         = ns;
SymModel.nq         = nq;
SymModel.nu         = nu;
SymModel.nx         = nx;
SymModel.nr         = nr;

SymModel.vSyms      = vSyms;
SymModel.vNames     = vNames;
SymModel.dv         = zeros(nv,1) + 3; % TODO: from units
SymModel.v          = v;

SymModel.kSyms      = kSyms;
SymModel.kNames     = kNames;
SymModel.k          = k;

SymModel.sSyms      = sSyms;
SymModel.sNames     = sNames;
SymModel.s          = s;

SymModel.qSyms      = qSyms;
SymModel.qNames     = qNames;
SymModel.q          = q;

SymModel.uSyms      = uSyms;
SymModel.uNames     = uNames;
SymModel.vuInd      = vuInd;
SymModel.u          = u;

SymModel.xSyms      = xSyms;
SymModel.xNames     = xNames;
SymModel.vxInd      = vxInd;
SymModel.x0         = x0;

SymModel.rNames     = rNames;
SymModel.r          = r;
SymModel.S          = S;
SymModel.Su         = Su;
