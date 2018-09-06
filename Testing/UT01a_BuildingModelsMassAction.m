function tests = UT01a_BuildingModelsMassAction()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
end

function testInitializingModel(a)
model_name = 'TestName';
m = InitializeModelMassActionAmount('TestName');
a.verifyEqual(m.Name, model_name);
end

function testInitializingModelWithEmptyName(a)
m = InitializeModelMassActionAmount();
a.verifyEqual(m.Name, '');
end

function testAddCompartment0(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'test', 0, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 0);
end

function testNonUnityZeroCompartment(a)
m = InitializeModelMassActionAmount();
a.verifyError(@()AddCompartment(m, 'test', 0, 2), 'KroneckerBio:Compartment:NonUnityZeroCompartment');
end

function testAddCompartment1(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'test', 1, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 1);
end

function testAddCompartment2(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'test', 2, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 2);
end

function testAddCompartment3(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'test', 3, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 3);
end

function testAddCompartment4(a)
m = InitializeModelMassActionAmount();
a.verifyError(@()AddCompartment(m, 'test', 4, 1), 'KroneckerBio:Compartment:Dimension');
end

function testAddCompartmentTwoWays(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'v1', 3, 2);
m = AddCompartment(m, 'v2', 3, {'', 2});
m = FinalizeModel(m);

v = m.v(0, m.x0(m.s), m.u);

a.verifyEqual(v(1), v(2))
end

function [m, x0, u0] = model_with_some_species()
m = InitializeModelMassActionAmount();
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

function testDefaultSpecies(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'v', 3, 1);
m = AddInput(m, 'u', 'v');
m = AddState(m, 'x', 'v');
m = FinalizeModel(m);
a.verifyEqual(m.u(), 0)
a.verifyEqual(m.x0(m.s), 0)
end

function testAddReactionFwd(a)
[m, x0, u0] = model_with_some_species();
m = AddReaction(m, 'test', {'x1', 'x2'}, 'x3', 'k1');
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [-24;-24;24;0])
end

function testAddReactionRev(a)
[m, x0, u0] = model_with_some_species();
m = AddReaction(m, 'test', {'x1', 'x2'}, 'x3', '', 'k2');
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [15;15;-15;0])
end

function testRepeatedComponents(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'v1', 3, 1);
m = AddParameter(m, 'k', 4);
m = AddSeed(m, 's1', 1);
m = AddState(m, 'x1', 'v1', 1);
m = AddInput(m, 'u1', 'v1', 1);
m = AddOutput(m, 'o1', 'x1');
m = AddReaction(m, 'r1', 'x1', {}, 'k');

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

test = AddReaction(m, 'r1', 'x1', {}, 'k');
a.verifyWarning(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:RepeatReactions'); % note: throws warning and ignores
end

function testRepeatedReactionsInputs(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'v1', 3, 1);
m = AddParameter(m, 'k', 4);
m = AddState(m, 'x1', 'v1', 1);
m = AddInput(m, 'u1', 'v1', 1);
m = AddInput(m, 'u2', 'v1', 1);
m = AddReaction(m, 'r1', {'x1', 'u1'}, {}, 'k');
m = AddReaction(m, 'r1', {'x1', 'u2'}, {}, 'k');

a.verifyWarningFree(@()FinalizeModel(m))
end

function testDuplicateSpeciesNames(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'v1', 3, 1);
m = AddCompartment(m, 'v2', 3, 1);
m = AddParameter(m, 'k', 4);
m = AddState(m, 'x1', 'v1');
m = AddState(m, 'x1', 'v2');
m = AddState(m, 'x2', 'v2');

a.verifyError(@()AddReaction(m, '', {'', 'x1'}, {}, 'k'), 'KroneckerBio:fixReactionSpecies:InvalidBlankName')

a.verifyError(@()AddReaction(m, '', {1, 'x1'}, {}, 'k'), 'KroneckerBio:fixReactionSpecies:InvalidName')

a.verifyError(@()AddReaction(m, '', 'v1.x1.y1', {}, 'k'), 'KroneckerBio:fixReactionSpecies:InvalidNameDots')

a.verifyError(@()AddReaction(m, '', 'v1.x1"y1', {}, 'k'), 'KroneckerBio:fixReactionSpecies:InvalidNameDoubleQuotes')

test = AddReaction(m, '', 'v1.x1', {}, 'k');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1) % no error on finalization

test = AddReaction(m, '', 'x1', {}, 'k', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', {}, 'x1', 'k', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'x1', {}, 'k', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'v1.x2', {}, 'k');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingQualifiedSpeciesName')

test = AddReaction(m, '', 'x3', {}, 'k');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingUnqualifiedSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:AmbiguousSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:AmbiguousSpeciesName')

test = AddReaction(m, '', 'x2', {}, 'k', '', 'v1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingSpeciesInReactionCompartment')

test = AddReaction(m, '', 'x1', {}, 'k', '', 'v3');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingReactionCompartment')
end

function testAddReactionsWithEmptyNames(a)
% Test various ways that reactions with empty reaction and product species lists
% are added. Does not test finalization.
m = InitializeModelMassActionAmount();
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

function testOutputsWithConstants(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'v1', 3, 1);
m = AddState(m, 'x1', 'v1');
m = AddState(m, 'x2', 'v1');
m = AddInput(m, 'u1', 'v1');
m = AddOutput(m, 'y1', {'', 10});
m = AddOutput(m, 'y2', {'x1', 2; '', 1});
m = AddOutput(m, 'y3', {'x2', 2; '', 1; 'u1', 3});
m = FinalizeModel(m);

a.verifyEqual(m.y(0, [3;4], 5), [10; 7; 24])
end

function testRemoveComponent(a)
m = InitializeModelMassActionAmount();
m = AddCompartment(m, 'v1', 3, 1);
m = AddParameter(m, 'k1', 4);
m = AddSeed(m, 's1', 1);
m = AddState(m, 'x1', 'v1', 1);
m = AddInput(m, 'u1', 'v1', 1);
m = AddOutput(m, 'y1', 'x1');
m = AddReaction(m, 'r1', 'x1', {}, 'k1');
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
end

function testAddOutputAsRegex(a)
m = model_with_some_species();
m = addOutputAsRegex(m, 'y', {'x'});
m = FinalizeModel(m);
a.verifyEqual(size(m.Outputs(1).Expression), [4,2])
end
