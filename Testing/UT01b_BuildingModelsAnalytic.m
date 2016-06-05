function tests = UT01b_BuildingModelsAnalytic()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testInitializingModel(a)
model_name = 'TestAnalytic';
m = InitializeModelAnalytic('TestAnalytic');
a.verifyEqual(m.Name, model_name);
end

function testInitializingModelWithEmptyName(a)
m = InitializeModelAnalytic();
a.verifyEqual(m.Name, '');
end

function testAddCompartment0(a)
m = InitializeModelAnalytic();
m = AddCompartment(m, 'test', 0, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 0);
end

function testNonUnityZeroCompartment(a)
m = InitializeModelAnalytic();
a.verifyError(@()AddCompartment(m, 'test0', 0, 2), 'KroneckerBio:Compartment:NonUnityZeroCompartment');
end

function testAddCompartments(a)
m = InitializeModelAnalytic();
m = AddCompartment(m, 'test1', 1, 1);
m = AddCompartment(m, 'test2', 2, 1);
m = AddCompartment(m, 'test3', 3, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 1);
a.verifyEqual(m.Compartments(2).Dimension, 2);
a.verifyEqual(m.Compartments(3).Dimension, 3);
end

function testAddCompartment4(a)
m = InitializeModelAnalytic();
a.verifyError(@()AddCompartment(m, 'test4', 4, 1), 'KroneckerBio:Compartment:Dimension');
end

function testStoichiometriyMatrixInRightOrder(a)
m = InitializeModelAnalytic();
m = AddCompartment(m, 'v1', 3, 1);
m = AddState(m, 'x1', 'v1', 2);
m = AddInput(m, 'u1', 'v1', 3);
m = AddParameter(m, 'k1', 4);
m = AddReaction(m, '', 'x1', 'u1', 'x1*k1');
m = FinalizeModel(m);

a.verifyEqual(m.f(0, 5, 3), -5*4)
end

function [m, x0, u0] = model_with_some_species()
m = InitializeModelAnalytic();
m = AddCompartment(m, 'v1', 3, 1);
m = AddState(m, 'x1', 'v1', 2);
m = AddState(m, 'x2', 'v1', 3);
m = AddState(m, 'x3', 'v1', 5);
m = AddState(m, 'x4', 'v1', 7);
m = AddParameter(m, 'k1', 4);
m = AddParameter(m, 'k2', 3);
m = FinalizeModel(m);
x0 = m.x0(m.s);
u0 = vec([]);
end

function testAddReactionFwd(a)
[m, x0, u0] = model_with_some_species();
m = AddReaction(m, 'test', {'x1', 'x2'}, 'x3', 'k1*x1*x2'); % must be specified explicitly in analytic model
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [-24;-24;24;0])
end

function testAddReactionRev(a)
[m, x0, u0] = model_with_some_species();
m = AddReaction(m, 'test', {'x1', 'x2'}, 'x3', '', 'k2*x3');
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [15;15;-15;0])
end

function testAddRule(a)
[m, x0, u0] = model_with_some_species();
m = AddRule(m, 'z1', 'sqrt(x1)');
m = AddReaction(m, 'test', 'x1', 'x2', 'z1');
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [-sqrt(2);sqrt(2);0;0])
end

function testRepeatedComponents(a)
% Consider adding repeated reaction and rule tests
m = InitializeModelAnalytic();
m = AddCompartment(m, 'v1', 3, 1);
m = AddParameter(m, 'k', 4);
m = AddSeed(m, 's1', 1);
m = AddState(m, 'x1', 'v1', 1);
m = AddInput(m, 'u1', 'v1', 1);
m = AddOutput(m, 'o1', 'x1');

test = AddCompartment(m, 'v1', 3, 1);
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:RepeatCompartment');

test = AddParameter(m, 'k', 4);
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:RepeatParameter');

test = AddSeed(m, 's1', 1);
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:RepeatSeed');

test = AddState(m, 'x1', 'v1', 1);
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:RepeatSpecies');

test = AddInput(m, 'u1', 'v1', 1);
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:RepeatSpecies');

test = AddOutput(m, 'o1', 'x1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:RepeatOutput');
end

function testDuplicateSpeciesNames(a)
m = InitializeModelAnalytic();
m = AddCompartment(m, 'v1', 3, 1);
m = AddCompartment(m, 'v2', 3, 1);
m = AddParameter(m, 'k', 4);
m = AddState(m, 'x1', 'v1');
m = AddState(m, 'x1', 'v2');
m = AddState(m, 'x2', 'v2');

a.verifyError(@()AddReaction(m, '', {'', 'x1'}, {}, 'k*x1'), 'KroneckerBio:fixReactionSpecies:InvalidBlankName')

a.verifyError(@()AddReaction(m, '', {1, 'x1'}, {}, 'k*x1'), 'KroneckerBio:fixReactionSpecies:InvalidName')

a.verifyError(@()AddReaction(m, '', 'v1.x1.y1', {}, 'k*x1'), 'KroneckerBio:fixReactionSpecies:InvalidNameDots')

a.verifyError(@()AddReaction(m, '', 'v1.x1"y1', {}, 'k*x1'), 'KroneckerBio:fixReactionSpecies:InvalidNameDoubleQuotes')

test = AddReaction(m, '', 'v1.x1', {}, 'k*v1.x1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1) % no error on finalization

test = AddReaction(m, '', 'x1', {}, 'k*x1', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', {}, 'x1', 'k*x1', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'x1', {}, 'k*v1.x1', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'x1', {}, 'k*v2.x2', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'v1.x2', {}, 'k*v1.x2');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingQualifiedSpeciesName')

test = AddReaction(m, '', 'x3', {}, 'k*x3');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingUnqualifiedSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k*x1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:AmbiguousSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k*v1.x1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:AmbiguousSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k*x2', '', 'v1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingSpeciesInReactionCompartment')

test = AddReaction(m, '', 'x2', {}, 'k*x2', '', 'v1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingSpeciesInReactionCompartment')

test = AddReaction(m, '', 'x1', {}, 'k*x1', '', 'v3');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingReactionCompartment')
end

function testExtraQuotesWarnings(a)
m = InitializeModelAnalytic();
m = AddCompartment(m, 'v1-', 3, 1);
m = AddCompartment(m, 'v2', 3, 1);
m = AddState(m, 'x1*', 'v1');
m = AddState(m, 'x2', 'v1');
m = AddParameter(m, 'k', 4);

a.verifyWarning(@()AddReaction(m, '', 'x1*', '', 'k*"v1-"."x1*"'), 'KroneckerBio:UnnecessaryQuotedSpecies')
a.verifyWarningFree(@()AddReaction(m, '', 'x1*', '', 'k*"v1-.x1*"'))
a.verifyWarning(@()AddReaction(m, '', 'x2', '', 'k*"v2.x2"'), 'KroneckerBio:UnnecessaryQuotedSpecies')
a.verifyWarning(@()AddReaction(m, '', 'x2', '', 'k*"v2"."x2"'), 'KroneckerBio:UnnecessaryQuotedSpecies')
a.verifyWarningFree(@()AddReaction(m, '', 'x2', '', 'k*v2.x2'))
a.verifyWarningFree(@()AddReaction(m, '', 'x1*', '', 'k*"x1*"'))
a.verifyWarning(@()AddReaction(m, '', 'x2', '', 'k*"x2"'), 'KroneckerBio:UnnecessaryQuotedSpecies')
a.verifyWarningFree(@()AddReaction(m, '', 'x2', '', 'k*x2'))

% Not bothering to check all combinations because they all call the same code
a.verifyWarning(@()AddCompartment(m, 'v3', 3, '"v1-"."x1*"'), 'KroneckerBio:UnnecessaryQuotedSpecies')
a.verifyWarning(@()AddRule(m, 'z1', 'k*k*"v1-"."x1*"'), 'KroneckerBio:UnnecessaryQuotedSpecies')
a.verifyWarning(@()AddOutput(m, 'z1', 'k*k*"v1-"."x1*"'), 'KroneckerBio:UnnecessaryQuotedSpecies')
end

function testAddReactionsWithEmptyNames(a)
% Test various ways that reactions with empty reaction and product species lists
% are added. Does not test finalization.
m = InitializeModelAnalytic();
m = AddCompartment(m, 'v1', 3, 1);
m = AddParameter(m, 'k1', 4);
m = AddParameter(m, 'k2', 3);
m = AddState(m, 'x1', 'v1', 2);
m = AddState(m, 'x2', 'v1', 3);

%% Forms of empty species lists
emptySpeciesForms = {[], {}, '', cell(1,0)};
nForms = length(emptySpeciesForms);
nr = m.nr;

% Empty reactants
for i = 1:nForms
    m = AddReaction(m, 'r1', emptySpeciesForms{i}, 'x2', 'k1');
    a.verifyEqual(m.Reactions(nr+i).Reactants, cell(1,0))
end

% Empty products
for i = 1:nForms
    m = AddReaction(m, 'r1', 'x1', emptySpeciesForms{i}, 'k1');
    a.verifyEqual(m.Reactions(nr+nForms+i).Products, cell(1,0))
end

%% Forms of some empty species lists that should error
emptySpeciesForms = {{[]}, {''}};
nForms = length(emptySpeciesForms);

% Empty reactants
for i = 1:nForms
    a.verifyError(@()AddReaction(m, 'r1', emptySpeciesForms{i}, 'x2', 'k1'), 'KroneckerBio:fixReactionSpecies:InvalidBlankName');
end

% Empty products
for i = 1:nForms
    a.verifyError(@()AddReaction(m, 'r1', 'x1', emptySpeciesForms{i}, 'k1'), 'KroneckerBio:fixReactionSpecies:InvalidBlankName');
end
end

function testsymbolic2PseudoKroneckerMM(a)

m = michaelis_menten_model();

a.verifyEqual(m.nv, 1)
a.verifyEqual(m.nk, 2)
a.verifyEqual(m.ns, 1)
a.verifyEqual(m.nu, 1)
a.verifyEqual(m.nx, 2)
a.verifyEqual(m.ny, 3)

end

function testRemoveComponent(a)
m = InitializeModelAnalytic();
m = AddCompartment(m, 'v1', 3, 1);
m = AddParameter(m, 'k1', 4);
m = AddSeed(m, 's1', 1);
m = AddState(m, 'x1', 'v1', 1);
m = AddInput(m, 'u1', 'v1', 1);
m = AddOutput(m, 'y1', 'x1');
m = AddReaction(m, 'r1', 'x1', {}, 'k1*x1');
m = AddRule(m, 'z1', 'x1*x1');
m = FinalizeModel(m);

test = RemoveCompartment(m, 'v1');
a.verifyEqual(test.nv, 0);
a.verifyEqual(length(test.Compartments), 0);
a.verifyError(@()RemoveCompartment(m, 'v2'), 'KroneckerBio:RemoveCompartment:CompartmentNotFound');

test = RemoveParameter(m, 'k1');
a.verifyEqual(test.nk, 0);
a.verifyEqual(length(test.Parameters), 0);
a.verifyError(@()RemoveParameter(m, 'k2'), 'KroneckerBio:RemoveParameter:ParameterNotFound');

test = RemoveSeed(m, 's1');
a.verifyEqual(test.ns, 0);
a.verifyEqual(length(test.Seeds), 0);
a.verifyError(@()RemoveSeed(m, 's2'), 'KroneckerBio:RemoveSeed:SeedNotFound');

test = RemoveState(m, 'x1');
a.verifyEqual(test.nx, 0);
a.verifyEqual(length(test.States), 0);
a.verifyError(@()RemoveState(m, 'x2'), 'KroneckerBio:RemoveState:StateNotFound');

test = RemoveInput(m, 'u1');
a.verifyEqual(test.nu, 0);
a.verifyEqual(length(test.Inputs), 0);
a.verifyError(@()RemoveInput(m, 'u2'), 'KroneckerBio:RemoveInput:InputNotFound');

test = RemoveOutput(m, 'y1');
a.verifyEqual(test.ny, 0);
a.verifyEqual(length(test.Outputs), 0);
a.verifyError(@()RemoveOutput(m, 'y2'), 'KroneckerBio:RemoveOutput:OutputNotFound');

test = RemoveReaction(m, 'r1');
a.verifyEqual(test.nr, 0);
a.verifyEqual(length(test.Reactions), 0);
a.verifyError(@()RemoveReaction(m, 'r2'), 'KroneckerBio:RemoveReaction:ReactionNotFound');

test = RemoveRule(m, 'z1');
a.verifyEqual(test.nz, 0);
a.verifyEqual(length(test.Rules), 0);
a.verifyError(@()RemoveRule(m, 'z2'), 'KroneckerBio:RemoveRule:RuleNotFound');
end

function testAddComponentsAsOutputs(a)
[m, x0, u0] = model_with_some_species();
m = AddInput(m, 'u1', 'v1');
m = AddInput(m, 'u2', 'v1');
m = AddRule(m, 'z1', 'x1*u1');
m = AddRule(m, 'z2', 'x2*u2');
m = FinalizeModel(m);

test = addStatesAsOutputs(m);
a.verifyEqual(strcat('v1.', {m.States.Name}), {test.Outputs.Name})
a.verifyEqual(strcat('v1.', {m.States.Name}), {test.Outputs.Expression})

test = addStatesAsOutputs(m, false);
a.verifyEqual({m.States.Name}, {test.Outputs.Name})

test = addInputsAsOutputs(m);
a.verifyEqual(strcat('v1.', {m.Inputs.Name}), {test.Outputs.Name})
a.verifyEqual(strcat('v1.', {m.Inputs.Name}), {test.Outputs.Expression})

test = addInputsAsOutputs(m, false);
a.verifyEqual({m.Inputs.Name}, {test.Outputs.Name})

test = addRulesAsOutputs(m);
a.verifyEqual({m.Rules.Name}, {test.Outputs.Name})
a.verifyEqual({m.Rules.Name}, {test.Outputs.Expression})
end
