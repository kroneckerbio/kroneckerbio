function tests = UT05_Sensitivity()
tests = functiontests(localfunctions);
end

function testEquilibrium(a)
m = LoadModel('Equilibrium.txt');
con = InitialValueExperiment(m, 10);
tGet = 1:10;
opts.Verbose = false;

verifySensitivity(a, m, con, tGet, opts)
end

function testDoseModel(a)
[m, con, unused, opts] = dose_model();

sim = SimulateSelect(m, con, 6, opts);

a.verifyEqual(sim(1).y(:,end), [0;8;0])
a.verifyEqual(sim(2).y(:,end), [0;4;0])
end

function testSimulateSensitivitySimple(a)
[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
tGet = 1:6;

verifySensitivity(a, m, con, tGet, opts)
end

function verifySensitivity(a, m, con, tGet, opts)
sim1 = SimulateSensitivity(m, con, opts);

sim2 = SimulateSensitivitySelect(m, con, tGet, opts);

sim3 = FiniteSimulateSensitivitySelect(m, con, tGet, opts);

a.verifyEqual(sim1.dydT(tGet,1:m.ny), sim3.dydT, 'RelTol', 0.01, 'AbsTol', 1e-4)
a.verifyEqual(sim2.dydT, sim3.dydT, 'RelTol', 0.01, 'AbsTol', 1e-4)
end
