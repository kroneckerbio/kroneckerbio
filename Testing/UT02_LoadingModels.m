function tests = UT02_LoadingModels()
tests = functiontests(localfunctions);
end

function testEquilibriumLoading(a)
m = LoadModelMassAction('Equilibrium.txt');
a.verifyEqual(m.Name, 'Equilibrium');
end

function testSimpleLoading(a)
m = LoadModelMassAction('Simple.txt');
a.verifyEqual(m.Name, 'Simple');
end