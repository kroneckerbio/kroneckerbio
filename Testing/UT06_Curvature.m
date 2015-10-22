function tests = UT06_Curvature()
tests = functiontests(localfunctions);
if nargout < 1
    tests.run;
end
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
[m, con, ~, opts] = dose_model();
tGet = 1:6;

verifyCurvature(a, m, con(1), tGet, opts)
verifyCurvature(a, m, con(2), tGet, opts)
verifyCurvature(a, m, con(3), tGet, opts)
end

function testHigherOrderDose(a)
[m, con, ~, opts] = higher_order_dose_model();
tGet = 0:10;
verifyCurvature(a, m, con, tGet, opts)
end

function testSimulateCurvatureSimple(a)
[m, con, ~, opts] = simple_model();
tGet = 1:6;

verifyCurvature(a, m, con, tGet, opts)
end

function testSimulateCurvatureSimpleEvent(a)
[m, con, ~, opts] = simple_model();

eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

verifyCurvatureEvent(a, m, con, obs, opts)
end

function testSimulateCurvatureSimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, ~, opts] = simple_model(simpleopts);
tGet = 1:6;

verifyCurvature(a, m, con, tGet, opts)
end

function testSimulateCurvatureMichaelisMenten(a)
[m, con, ~, opts] = michaelis_menten_model();
tGet = 1:10;

verifyCurvature(a, m, con, tGet, opts)
end

function testSimulateCurvatureMichaelisMentenEvent(a)
[m, con, ~, opts, eve] = michaelis_menten_model();

obs = observationEvents(10, eve);

verifyCurvatureEvent(a, m, con, obs, opts)
end

function testSimulateCurvatureSimpleAnalytic(a)
[m, con, ~, opts] = simple_analytic_model();
tGet = 0:5:40;

verifyCurvature(a, m, con, tGet, opts)
end

function verifyCurvature(a, m, con, tGet, opts)
nT = nnz(opts.UseParams)+nnz(opts.UseSeeds)+nnz(opts.UseInputControls)+nnz(opts.UseDoseControls);
obsSelect = observationSelect(tGet);
opts.ComplexStep = true;

opts.Normalized = false;
sim1 = SimulateCurvature(m, con, max(tGet), opts);
sim2 = SimulateCurvature(m, con, obsSelect, opts);
sim3 = FiniteSimulateCurvature(m, con, obsSelect, opts);

a.verifyEqual(size(sim2.d2xdT2), [m.nx*nT*nT, numel(tGet)])
a.verifyEqual(size(sim2.d2udT2), [m.nu*nT*nT, numel(tGet)])
a.verifyEqual(size(sim2.d2ydT2), [m.ny*nT*nT, numel(tGet)])

a.verifyEqual(sim1.d2ydT2(tGet,1:m.ny), sim3.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(sim2.d2ydT2, sim3.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4)

opts.Normalized = true;
sim4 = SimulateCurvature(m, con, max(tGet), opts);
sim5 = SimulateCurvature(m, con, obsSelect, opts);
sim6 = FiniteSimulateCurvature(m, con, obsSelect, opts);

a.verifyEqual(sim4.d2ydT2(tGet,1:m.ny), sim6.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(sim5.d2ydT2, sim6.d2ydT2, 'RelTol', 0.001, 'AbsTol', 1e-4)
end

function verifyCurvatureEvent(a, m, con, obs, opts)
opts.ComplexStep = true;

opts.Normalized = false;
sim1 = SimulateCurvature(m, con, obs, opts);
sim2 = FiniteSimulateCurvature(m, con, obs, opts);

a.verifyEqual(sim1.d2yedT2, sim2.d2yedT2, 'RelTol', 0.001, 'AbsTol', 1e-4)

opts.Normalized = true;
sim3 = SimulateCurvature(m, con, obs, opts);
sim4 = FiniteSimulateCurvature(m, con, obs, opts);

a.verifyEqual(sim3.d2yedT2, sim4.d2yedT2, 'RelTol', 0.001, 'AbsTol', 1e-4)
end
