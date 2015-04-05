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

function testNonUnityZeroCompartment(a)
m = InitializeModel();
a.verifyError(@()AddCompartment(m, 'test', 0, 2), 'KroneckerBio:Compartment:NonUnityZeroCompartment');
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

function testsymbolic2PseudoKroneckerMM(a)
syms Km kcat S0 E S P

m.kNames = {'Km';'kcat'};
m.kSyms = [kcat;Km];
m.k = [10;2];

m.sNames = {'S0'};
m.sSyms = S0;
m.s = 30;

m.uNames = {'E'};
m.uSyms = E;
m.u = 1;

m.xNames = {'S';'P'};
m.xSyms = [S;P];
m.x0 = [S0;0];

r = kcat*E*S/(Km+S);
m.f = [-r;r];

m.yNames = {'S','P'};
m.y = [S;P];

m_kron = symbolic2PseudoKronecker(m);

end
