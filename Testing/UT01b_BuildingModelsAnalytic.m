function tests = UT01b_BuildingModelsAnalytic()
tests = functiontests(localfunctions);
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
m = AddReaction(m, 'test', 'x1', 'x2', 'x3', '', 'k1*x1*x2'); % must be specified explicitly in analytic model
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [-24;-24;24;0])
end

function testAddReactionRev(a)
[m, x0, u0] = model_with_some_species();
m = AddReaction(m, 'test', 'x1', 'x2', 'x3', '', '', 'k2*x3');
m = FinalizeModel(m);
a.verifyEqual(m.f(0,x0,u0), [15;15;-15;0])
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