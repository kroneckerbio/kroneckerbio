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

test = AddReaction(m, '', 'v1.x1', {}, 'k*"v1.x1"');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1) % no error on finalization

test = AddReaction(m, '', 'x1', {}, 'k*x1', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', {}, 'x1', 'k*x1', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'x1', {}, 'k*"v1.x1"', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'x1', {}, 'k*"v2.x2"', '', 'v1');
test = FinalizeModel(test);
a.verifyEqual(test.nr, 1)

test = AddReaction(m, '', 'v1.x2', {}, 'k*v1.x2');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingQualifiedSpeciesName')

test = AddReaction(m, '', 'x3', {}, 'k*x3');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingUnqualifiedSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k*x1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:AmbiguousSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k*"v1.x1"');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:AmbiguousSpeciesName')

test = AddReaction(m, '', 'x1', {}, 'k*x2', '', 'v1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingSpeciesInReactionCompartment')

test = AddReaction(m, '', 'x2', {}, 'k*"x2"', '', 'v1');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingSpeciesInReactionCompartment')

test = AddReaction(m, '', 'x1', {}, 'k*"x1"', '', 'v3');
a.verifyError(@()FinalizeModel(test), 'KroneckerBio:FinalizeModel:MissingReactionCompartment')
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

% Empty reactants
for i = 1:nForms
    m = AddReaction(m, 'r1', emptySpeciesForms{i}, 'x2', 'k1');
    a.verifyEqual(m.add.Reactions(i).Reactants, cell(1,0))
end

% Empty products
for i = 1:nForms
    m = AddReaction(m, 'r1', 'x1', emptySpeciesForms{i}, 'k1');
    a.verifyEqual(m.add.Reactions(nForms+i).Products, cell(1,0))
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