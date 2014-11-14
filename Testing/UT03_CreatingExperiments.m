function tests = UT03_CreatingExperiments()
tests = functiontests(localfunctions);
end

function testExperiment1(a)
m = LoadModelMassAction('Equilibrium.txt');
con = Experiment(m, 10);
a.verifyEqual(con.tF, 10);
a.verifyEqual(con.Periodic, false);
a.verifyEqual(con.SteadyState, false);
end

function testExperimentSeed(a)
m = LoadModelMassAction('Equilibrium.txt');
s = vec(1:m.ns);
con = Experiment(m, 10, s);
a.verifyEqual(con.s, s);
end

function testExperimentInputConstant(a)
m = LoadModelMassAction('Equilibrium.txt');
u = 1:m.nu;
con = Experiment(m, 10, [], false, false, u);
a.verifyEqual(con.u(-1), zeros(m.nu,1));
a.verifyEqual(con.u(5), vec(1:m.nu));
end

function testExperimentDoseConstant(a)
m = LoadModelMassAction('Equilibrium.txt');
d = 1:m.nu;
con = Experiment(m, 10, [], false, false, d);
a.verifyEqual(con.u(-1), zeros(m.nu,1));
a.verifyEqual(con.u(5), vec(1:m.nu));
end
