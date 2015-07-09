function tests = UT06_Curvature()
tests = functiontests(localfunctions);
end

function testEquilibrium(a)
m = LoadModel('Equilibrium.txt');
con = experimentInitialValue(m);
tGet = 1:10;
opts.Verbose = false;
opts.UseParams = 1:m.nk;
opts.UseSeeds = [];
opts.UseInputControls = [];
opts.UseDoseControls = [];

verifyCurvature(a, m, con, tGet, opts)
end

function testDoseModel(a)
[m, con, unused, opts] = dose_model();
tGet = 1:6;

verifyCurvature(a, m, con(1), tGet, opts)
verifyCurvature(a, m, con(2), tGet, opts)
verifyCurvature(a, m, con(3), tGet, opts)
end

function testSimulateCurvatureSimple(a)
[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
tGet = 1:6;

verifyCurvature(a, m, con, tGet, opts)
end

function testSimulateCurvatureSimpleEvent(a)
[m, con, unused, opts] = simple_model();
m = m.Update(rand(m.nk,1)+1);
con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);

eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

verifyCurvatureEvent(a, m, con, obs, opts)
end

function testSimulateCurvatureSimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, unused, opts] = simple_model(simpleopts);
% m = m.Update(rand(m.nk,1)+1);
% con = con.Update(rand(con.ns,1)+1, rand(con.nq,1)+1, rand(con.nh,1)+1);
tGet = 1:6;

verifyCurvature(a, m, con, tGet, opts)
end

function testSimulateCurvatureMichaelisMenten(a)
[m, con, ~, opts] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);
tGet = 1:10;

verifyCurvature(a, m, con, tGet, opts)
end

function testSimulateCurvatureMichaelisMentenEvent(a)
[m, con, ~, opts, eve] = michaelis_menten_model();
m = m.Update(rand(m.nk,1)+m.k+1);
con = con.Update(rand(con.ns,1)+con.s+1, rand(con.nq,1)+con.q+1, rand(con.nh,1)+con.h+1);

obs = observationEvents(10, eve);

verifyCurvatureEvent(a, m, con, obs, opts)
end

function verifyCurvature(a, m, con, tGet, opts)
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);
obsSelect = observationSelect(tGet);

sim1 = SimulateCurvature(m, con, max(tGet), opts);

sim2 = SimulateCurvature(m, con, obsSelect, opts);

sim3 = FiniteSimulateCurvature(m, con, obsSelect, opts);

a.verifyEqual(size(sim2.d2xdT2), [m.nx*nT*nT, numel(tGet)])
a.verifyEqual(size(sim2.d2udT2), [m.nu*nT*nT, numel(tGet)])
a.verifyEqual(size(sim2.d2ydT2), [m.ny*nT*nT, numel(tGet)])

a.verifyEqual(sim1.d2ydT2(tGet,1:m.ny), sim3.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(sim2.d2ydT2, sim3.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4)
end

function verifyCurvatureEvent(a, m, con, obs, opts)
sim1 = SimulateCurvature(m, con, obs, opts);

sim2 = FiniteSimulateCurvature(m, con, obs, opts);

a.verifyEqual(sim1.d2yedT2, sim2.d2yedT2, 'RelTol', 0.001, 'AbsTol', 1e-4)
end
