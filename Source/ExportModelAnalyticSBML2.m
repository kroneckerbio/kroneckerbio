function ExportModelAnalyticSBML2(m, filename, opts)
%ExportModelAnalyticSBML2 Export kroneckerbio analytic model to SBML file
%   using libSBML's Matlab interface. Requires that the input model be
%   finalized.
%
% Inputs:
%   m [ Model.Analytic struct ]
%       Finalized kroneckerbio analytic model
%   filename [ string ]
%       Name of output SBML file
%   opts [ options struct ]
%       Options struct allowing the following fields:
%       .Verbose [ nonnegative integer {0} ]
%           Control printed information. Higher means more is printed.
%
% Notes:
% - Old versions of this code caused the libSBML mex functions to segfault
%   *A LOT*. This updated version doesn't seem to have those issues (so
%   far). A major change is that a dummy base model is loaded and its
%   struct fields overwritten. The old method built up the struct ab
%   initio, which may have caused these problems.
% - The input model must be finalized because that qualifies species with
%   compartments and other cleanup and validation.
% 
% TODO:
% - Break up the struct building part into a separate function for
%   easier unit testing
% - Do we use amounts or concs for species? Right now, it's all set to
%   amounts but I don't remember if this is right, or if we should have the
%   option for either. The old code supported the option.
% - Implement repeated assignment rules

%% Options
% Resolve missing inputs
if nargin < 3
    opts = [];
end

% Options for displaying progress
opts_.Verbose = 0;
opts = mergestruct(opts_, opts);
verbose = logical(opts.Verbose);

% Sanity check
assert(strcmp(m.Type, 'Model.Analytic'), 'ExportModelAnalyticSBML2:InvalidModelType', 'Input model must be a Model.Analytic')

assert(m.Ready, 'ExportModelAnalyticSBML2:UnfinalizedModel', 'Input model must be finalized')

if verbose; fprintf('Converting analytic model to SBML...'); end

%% Create name <-> ID mapping for all model components
nx = m.nx;
nu = m.nu;
nv = m.nv;
nk = m.nk;
ns = m.ns;

xNames = {m.States.Name};
uNames = {m.Inputs.Name};
vNames = {m.Compartments.Name};
kNames = {m.Parameters.Name};
sNames = {m.Seeds.Name};
names = [xNames uNames vNames kNames sNames];

xCompartments = {m.States.Compartment};
uCompartments = {m.Inputs.Compartment};
compartments = [xCompartments uCompartments];

xIds = strcat(repmat({'x'}, nx, 1), strtrim(cellstr(num2str((1:nx).'))))';
uIds = strcat(repmat({'u'}, nu, 1), strtrim(cellstr(num2str((1:nu).'))))';
vIds = strcat(repmat({'v'}, nv, 1), strtrim(cellstr(num2str((1:nv).'))))';
kIds = strcat(repmat({'k'}, nk, 1), strtrim(cellstr(num2str((1:nk).'))))';
sIds = strcat(repmat({'s'}, ns, 1), strtrim(cellstr(num2str((1:ns).'))))';
ids = [xIds uIds vIds kIds sIds];

xCompartmentIds = cell(1,nx);
for i = 1:nx
    xCompartmentIds(i) = vIds(ismember(vNames, xCompartments{i}));
end
uCompartmentIds = cell(1,nu);
for i = 1:nu
    uCompartmentIds(i) = vIds(ismember(vNames, uCompartments{i}));
end

% Qualify and quote (if necessary) species names
%   These appear in arbitrary expressions in reaction rates, state initial
%   conditions, outputs, and rules.
xuQualified = strcat([xCompartments, uCompartments]', repmat({'.'}, nx+nu, 1), [xNames, uNames]')';
for i = 1:nx+nu
xuIds = [xIds uIds];

allQualified = [xuQualified vNames kNames sNames];
allIds = [xuIds vIds kIds sIds];

%% Load dummy base struct
% Note/Debug: Can load a more complete SBML model to see what the structs
%   should look like during development
path = fileparts(mfilename('fullpath'));
sbml = TranslateSBML(fullfile(path, 'base_model.xml'));

%% Main SBML struct
% TODO: annotate with version of kroneckerbio version and commit - also it
%   doesn't look like sbml.notes is actually included in the output SBML
%   file
sbml.name  = m.Name;
sbml.id    = regexprep(m.Name, '\W', '');
sbml.notes = 'Converted from kroneckerbio model';

%% Compartments
vBase = sbml.compartment;
v = vBase;
v(1) = [];
for i = 1:nv
    compartment = m.Compartments(i);
    name = compartment.Name;
    dimension = compartment.Dimension;
    size = compartment.Size;
    if ischar(size)
        size = str2double(size);
    end
    
    v(i) = vBase;
    v(i).name = name;
    v(i).id = vIds{i};
    v(i).spatialDimensions = dimension;
    v(i).size = size;
end
sbml.compartment = v;

%% Parameters
% Rate parameters are regular parameters
kBase = sbml.parameter;
k = kBase;
k(1) = [];
for i = 1:nk
    parameter = m.Parameters(i);
    name = parameter.Name;
    value = parameter.Value;
    
    k(i) = kBase;
    k(i).name = name;
    k(i).id = kIds{i};
    k(i).value = value;
end

% Seeds are also treated as regular parameters
for i = 1:ns
    parameter = m.Seeds(i);
    name = parameter.Name;
    value = parameter.Value;
    
    j = m.nk + i;
    k(j) = kBase;
    k(j).name = name;
    k(j).id = sIds{i};
    k(j).value = value;
end
sbml.parameter = k;

%% Species
% States are modified by reactions
sBase = sbml.species;
s = sBase;
s(1) = [];
initialAssignmentsBase = sbml.initialAssignment;
initialAssignments = initialAssignmentsBase;
initialAssignments(1) = [];
initialAssignmentInd = 1;
for i = 1:nx
    species = m.States(i);
    name = species.Name;
    initial = str2double(species.InitialValue);
    initialConst = true;
    if isnan(initial) % didn't parse to a single constant number
        initial = 0;
        initialLHS = xIds{i};
        
        % RHS expression an only be a function of rate parameters and seeds
        %   so a simple sub is OK
        % TODO: May need to symbolically evaluate to take care of certain
        %   math functions
        initialRHS = substituteQuotedExpressions(species.InitialValue, names, ids);
        initialConst = false;
    end
    
    s(i) = sBase;
    s(i).name = name;
    s(i).id = xIds{i};
    s(i).compartment = xCompartmentIds{i};
    s(i).hasOnlySubstanceUnits = 1;
    s(i).isSetInitialConcentration = 0;
    s(i).isSetInitialAmount = 1;
    s(i).initialConcentration = NaN;
    if initialConst
        s(i).initialAmount = initial;
    else
        s(i).initialAmount = 0;
    end
    
    if ~initialConst
        initialAssignments(initialAssignmentInd) = initialAssignmentsBase;
        initialAssignments(initialAssignmentInd).symbol = initialLHS;
        initialAssignments(initialAssignmentInd).math = initialRHS;
        initialAssignmentInd = initialAssignmentInd + 1;
    end
end
sbml.initialAssignment = initialAssignments;

% Inputs aren't modified by reactions
% The const/default value of inputs are used as their initial amounts
% BoundaryCondition aren't modified by reactions but are modified
%   by rules. ConstantAmount aren't modified by anything.
for i = 1:nu
    species = m.Inputs(i);
    name = species.Name;
    initial = species.DefaultValue;
    
    j = m.nx + i;
    s(j) = sBase;
    s(j).name = name;
    s(j).id = uIds{i};
    s(j).compartment = uCompartmentIds{i};
    s(j).hasOnlySubstanceUnits = 1;
    s(j).isSetInitialConcentration = 0;
    s(j).isSetInitialAmount = 1;
    s(j).initialConcentration = NaN;
    s(j).initialAmount = initial;
    s(j).boundaryCondition = true;
end
sbml.species = s;

%% Reactions
rBase = sbml.reaction;
reactantBase = rBase.reactant;
kineticLawBase = rBase.kineticLaw;
r = rBase;
r(1) = [];
for i = 1:m.nr
    reaction = m.Reactions(i);
    name = reaction.Name;
    reactants_ = reaction.Reactants;
    products_ = reaction.Products;
    rate = substituteQuotedExpressions(reaction.Rate, allQualified, allIds);
    
    nR = length(reactants_);
    reactants = reactantBase;
    reactants(1) = [];
    for iR = 1:nR
        reactants(iR) = reactantBase;
        if ~isValidIdentifier(reactants_{iR});
            reactants_{iR} = ['"' reactants_{iR} '"'];
        end
        reactants(iR).species = substituteQuotedExpressions(reactants_{iR}, xuQualified, xuIds);
    end
    
    nP = length(products_);
    products = reactantBase;
    products(1) = [];
    for iP = 1:nP
        products(iP) = reactantBase;
        if ~isValidIdentifier(products_{iP});
            products_{iP} = ['"' products_{iP} '"'];
        end
        products(iP).species = substituteQuotedExpressions(products_{iP}, xuQualified, xuIds);
    end
    
    kineticLaw = kineticLawBase;
    kineticLaw.formula = rate;
    kineticLaw.math = rate;
    
    r(i) = rBase;
    r(i).name = name;
    r(i).id = sprintf('r%i', i);
    r(i).reactant = reactants;
    r(i).product = products;
    r(i).kineticLaw = kineticLaw;
end
sbml.reaction = r;

%% Repeated assignment rules
ruleBase = sbml.rule;
rules = ruleBase;
rules(1) = [];
for i = 1:m.nr
    % TODO: implement this
end
sbml.rule = rules;

%% Test for validity
[valid, mes] = isSBML_Model(sbml);
if ~valid
    error('ExportModelAnalyticSBML2:SbmlStructCreationInvalid', ['libSBML model validator failed with message: ' mes '. Please file a bug report.'])
end

%% Output SBML
OutputSBML(sbml, filename);

end

