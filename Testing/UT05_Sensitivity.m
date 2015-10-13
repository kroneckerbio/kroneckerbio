function tests = UT05_Sensitivity()
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

verifySensitivity(a, m, con, tGet, opts)
end

function testDoseModel(a)
[m, con, ~, opts] = dose_model();
tGet = 1:6;

verifySensitivity(a, m, con(1), tGet, opts)
verifySensitivity(a, m, con(2), tGet, opts)
verifySensitivity(a, m, con(3), tGet, opts)
end

function testHigherOrderDose(a)
[m, con, ~, opts] = higher_order_dose_model();
tGet = 0:10;
verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivitySimple(a)
[m, con, ~, opts] = simple_model();
tGet = 1:6;

verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivitySimpleSteadyState(a)
simpleopts.steadyState = true;
[m, con, ~, opts] = simple_model(simpleopts);
tGet = 1:6;

verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivitySimpleEvent(a)
[m, con, ~, opts] = simple_model();

eve1 = eventDropsBelow(m, 10, 15);
eve2 = eventDropsBelow(m, 1, 2);
obs = observationEvents(6, [eve1;eve2]);

verifySensitivityEvent(a, m, con, obs, opts)
end

function testSimulateSensitivitySimpleAnalytic(a)
[m, con, ~, opts] = simple_analytic_model();
tGet = 0:10:40;

verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivityMichaelisMenten(a)
[m, con, ~, opts] = michaelis_menten_model();
tGet = 1:10;

verifySensitivity(a, m, con, tGet, opts)
end

function testSimulateSensitivityMichaelisMentenEvent(a)
[m, con, ~, opts, eve] = michaelis_menten_model();

obs = observationEvents(10, eve);

verifySensitivityEvent(a, m, con, obs, opts)
end

function verifySensitivity(a, m, con, tGet, opts)
opts.ImaginaryStep = true;
obsSelect = observationSelect(tGet);

opts.Normalized = false;
sim1 = SimulateSensitivity(m, con, max(tGet), opts);
sim2 = SimulateSensitivity(m, con, obsSelect, opts);
sim3 = FiniteSimulateSensitivity(m, con, obsSelect, opts);

a.verifyEqual(sim1.dydT(tGet,1:m.ny), sim3.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(sim2.dydT, sim3.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4)

opts.Normalized = true;
sim4 = SimulateSensitivity(m, con, max(tGet), opts);
sim5 = SimulateSensitivity(m, con, obsSelect, opts);
sim6 = FiniteSimulateSensitivity(m, con, obsSelect, opts);

a.verifyEqual(sim4.dydT(tGet,1:m.ny), sim6.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4)
a.verifyEqual(sim5.dydT, sim6.dydT, 'RelTol', 0.001, 'AbsTol', 1e-4)
end

function verifySensitivityEvent(a, m, con, obs, opts)
opts.ImaginaryStep = true;

opts.Normalized = false;
sim1 = SimulateSensitivity(m, con, obs, opts);
sim2 = FiniteSimulateSensitivity(m, con, obs, opts);

a.verifyEqual(sim1.dyedT, sim2.dyedT, 'RelTol', 0.001, 'AbsTol', 1e-4)

opts.Normalized = true;
sim3 = SimulateSensitivity(m, con, obs, opts);
sim4 = FiniteSimulateSensitivity(m, con, obs, opts);

a.verifyEqual(sim3.dyedT, sim4.dyedT, 'RelTol', 0.001, 'AbsTol', 1e-4)
end
