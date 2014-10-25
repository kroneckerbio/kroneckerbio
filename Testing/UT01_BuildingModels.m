function tests = UT01_BuildingModels()
tests = functiontests(localfunctions);
end

function testInitializingModel(a)
model_name = 'TestName';
m = InitializeModel('TestName');
a.verifyEqual(m.Name, model_name);
end

function testInitializingModelWithEmptyName(a)
m = InitializeModel();
a.verifyEqual(m.Name, '');
end

function testAddCompartment0(a)
m = InitializeModel();
m = AddCompartment(m, 'test', 0, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 0);
end

function testAddCompartment1(a)
m = InitializeModel();
m = AddCompartment(m, 'test', 1, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 1);
end

function testAddCompartment2(a)
m = InitializeModel();
m = AddCompartment(m, 'test', 2, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 2);
end

function testAddCompartment3(a)
m = InitializeModel();
m = AddCompartment(m, 'test', 3, 1);
m = FinalizeModel(m);
a.verifyEqual(m.Compartments(1).Dimension, 3);
end

function testAddCompartment4(a)
m = InitializeModel();
a.verifyError(@()AddCompartment(m, 'test', 4, 1), 'KroneckerBio:Compartment:Dimension');
end
