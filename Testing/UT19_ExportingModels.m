function tests = UT19_ExportingModels()
% Warning: Finalizing the models is slow
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testKbName2SimBioName(a)
% This is sort of redundant with a test in UT15, but adds the 'a."b!"' case
inputs = {'abc', ...
    'comp.species', ...
    '"a@"', ...
    '"a@.b!"', ...
    '"a@".bc', ...
    'a."b!"'}; % this is also valid now

outputs = {'abc', ...
    'comp.species', ...
    '[a@]', ...
    '[a@].[b!]', ...
    '[a@].bc', ...
    'a.[b!]'};

for i = 1:length(inputs)
    out = kroneckerbioExpr2SimbioExpr(inputs{i});
    a.verifyEqual(out, outputs{i});
end
end

function testEquilibriumMassActionExport(a)
m = LoadModelMassAction('Equilibrium.txt'); % this function finalizes the model
validateExporter(a, m, 'massaction');
end

function testSimpleMassActionExport(a)
m = LoadModelMassAction('Simple.txt');
validateExporter(a, m, 'massaction');
end

function testTestModelExport(a)
m = LoadModelSbmlAnalytic('test.xml');
m = FinalizeModel(m);
validateExporter(a, m, 'analytic');
end

function testMMModelExport(a)
m = michaelis_menten_model();
validateExporter(a, m, 'analytic');
end

function testQuotedIdentifiersModelExport(a)
m = InitializeModelAnalytic('QuotedIdentifiersModel');
m = AddCompartment(m, 'v', 3, 1);
m = AddCompartment(m, 'v!', 3, 1);
m = AddInput(m, 'E#', 'v', 1);
m = AddState(m, 'S', 'v', 'S0^2');
m = AddState(m, 'P:P', 'v', 'P0 + 10');
m = AddState(m, 'D', 'v!', 1.1);
m = AddParameter(m, 'K_m', 10);
m = AddParameter(m, 'kc@', 2);
m = AddSeed(m, 'S0', 5);
m = AddSeed(m, 'P0', 1.5);
m = AddReaction(m, 'r1', 'S', 'P:P', '"kc@"*"E#"*S/(K_m+S)');
m = AddReaction(m, 'r2', 'D', [], '"kc@"*D');
m = FinalizeModel(m);

validateExporter(a, m, 'analytic');

end

%% Helper functions
function validateExporter(a, m, modelType)
% Exports m and tests it. SBML exporters are tested by importing, finalizing
%   model, and comparing. m must be finalized. modelType is either
%   'massaction' or 'analytic'.
% TODO: actually simulate the loaded models after export and see if they
%   give the same results
nv = m.nv;
nk = m.nk;
ns = m.ns;
nr = m.nr;
nx = m.nx;
nu = m.nu;

vNames = {m.Compartments.Name};
kNames = {m.Parameters.Name};
sNames = {m.Seeds.Name};
xNames = {m.States.Name};
uNames = {m.Inputs.Name};

% analytic -> simbio
switch modelType
    case 'massaction'
        simbio = ExportModelMassActionSimBio(m);
    case 'analytic'
        simbio = ExportModelAnalyticSimBio(m);
end

a.verifyEqual(length(simbio.Compartments), nv);
a.verifyEqual(length(simbio.Parameters), nk+ns);
a.verifyEqual(length(simbio.Reactions), nr);
a.verifyEqual(length(simbio.Species), nx+nu);

% analytic -> sbml using simbio
sbmlFilename = getRandomFilename('xml');
cleanup = onCleanup(@() cleanupFile(sbmlFilename)); % may want to disable this when writing the tests
switch modelType
    case 'massaction'
        ExportModelMassActionSBML(m, sbmlFilename);
    case 'analytic'
        ExportModelAnalyticSBML(m, sbmlFilename);
end
m2 = LoadModelSbmlAnalytic(sbmlFilename, struct('Validate', true)); % only the analytic loaders are fully working right now
m2 = FinalizeModel(m2);

validateLoadedSbmlModel(m2);

% analytic -> sbml using libsbml
sbmlFilename2 = getRandomFilename('xml');
cleanup2 = onCleanup(@() cleanupFile(sbmlFilename2));
switch modelType
    case 'massaction'
        return % not implemented yet
%         ExportModelMassActionSBML2(m, sbmlFilename2);
    case 'analytic'
        ExportModelAnalyticSBML2(m, sbmlFilename2); % don't know how to suppress the 'Document written' message
end
m3 = LoadModelSbmlAnalytic(sbmlFilename2, struct('Validate', true));
m3 = FinalizeModel(m3);

validateLoadedSbmlModel(m3)

    function validateLoadedSbmlModel(mLoaded)
        a.verifyEqual(mLoaded.nv, nv);
        a.verifyEqual(mLoaded.nk, nk+ns); % conversion turns seeds into rate parameters
        a.verifyEqual(mLoaded.nr, nr);
        a.verifyEqual(mLoaded.nx, nx);
        a.verifyEqual(mLoaded.nu, nu); % importer recognizes boundary condition species
        
        a.verifyEqual({mLoaded.Compartments.Name}, vNames);
        a.verifyEqual(sort({mLoaded.Parameters.Name}), sort([kNames, sNames])); % components may be reordered
        a.verifyEqual(sort({mLoaded.States.Name}), sort(xNames));
        a.verifyEqual(sort({mLoaded.Inputs.Name}), sort(uNames));
    end
end

function filename = getRandomFilename(ext)
symbols = ['a':'z' 'A':'Z' '0':'9'];
STR_LEN = 10;
nums = randi(numel(symbols), [1 STR_LEN]);
st = symbols (nums);
filename = [st '.' ext];
end

function cleanupFile(filename)
if exist(filename, 'file') == 2
    delete(filename);
end
end