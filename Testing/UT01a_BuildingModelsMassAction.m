function tests = UT01a_BuildingModelsMassAction()
tests = functiontests(localfunctions);
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

function testAddReactionFwd(a)
[m, x0, u0] = model_with_some_species();
m = AddReaction(m, 'test', 'x1', 'x2', 'x3', '', 'k1');
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [-24;-24;24;0])
end

function testAddReactionRev(a)
[m, x0, u0] = model_with_some_species();
m = AddReaction(m, 'test', 'x1', 'x2', 'x3', '', '', 'k2');
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [15;15;-15;0])
end
