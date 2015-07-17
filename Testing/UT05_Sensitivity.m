function tests = UT05_Sensitivity()
tests = functiontests(localfunctions);
end

function testEquilibrium(a)
m = LoadModel('Equilibrium.txt');
con = experimentInitialValue(m);
tGet = 1:10;
opts.Verbose = false;

verifySensitivity(a, m, con, tGet, opts)
end

function testDoseModel(a)
[m, con, unused, opts] = dose_model();
tGet = 1:6;

verifySensitivity(a, m, con(1), tGet, opts)
verifySensitivity(a, m, con(2), tGet, opts)
verifySensitivity(a, m, con(3), tGet, opts)
end

function testSimulateSensitivitySimple(a)
[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
tGet = 1:6;

verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivitySimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, unused, opts] = simple_model(simpleopts);
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
tGet = 1:6;

verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivitySimpleEvent(a)
[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

verifySensitivityEvent(a, m, con, obs, opts)
end

function testSimulateSensitivityMichaelisMenten(a)
[m, con, ~, opts] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);
tGet = 1:10;

verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivityMichaelisMentenEvent(a)
[m, con, ~, opts, eve] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);

obs = observationEvents(10, eve);

verifySensitivityEvent(a, m, con, obs, opts)
end

function verifySensitivity(a, m, con, tGet, opts)
obsSelect = observationSelect(tGet);

sim1 = SimulateSensitivity(m, con, max(tGet), opts);

sim2 = SimulateSensitivity(m, con, obsSelect, opts);

sim3 = FiniteSimulateSensitivity(m, con, obsSelect, opts);

a.verifyEqual(sim1.dydT(tGet,1:m.ny), sim3.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(sim2.dydT, sim3.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4)
end

function verifySensitivityEvent(a, m, con, obs, opts)
sim1 = SimulateSensitivity(m, con, obs, opts);

sim2 = FiniteSimulateSensitivity(m, con, obs, opts);

a.verifyEqual(sim1.dyedT, sim2.dyedT, 'RelTol', 0.001, 'AbsTol', 1e-4)
end
